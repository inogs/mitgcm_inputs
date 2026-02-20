import argparse
import json
import logging
from pathlib import Path

import copernicusmarine as cm
import numpy as np
import pandas as pd
import xarray as xr
from bitsea.commons.mask import Mask
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import path_inside_an_existing_dir

LOGGER = logging.getLogger()


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Reads sewage info on excel file and for all the domains
    Takes salinity data from MDS service
    build the json files with the coordinate
    of the pointsource for E.Coli tracers
    PointSource_NAD.json
    PointSource_SAD.json
    ...

    """
    )

    parser.add_argument(
        "--inputfile",
        "-i",
        type=existing_file_path,
        default="Allegato_1_Capitolato_Tecnico_B32_B35_scarichi.xlsx",
        required=False,
        help="Scarichi, sewage discharge locations ",
    )
    parser.add_argument(
        "--domain",
        "-d",
        type=str,
        required=True,
        help="Domain name (e.g. NAD, SAD, etc.)",
    )
    parser.add_argument(
        "--maskfile",
        "-m",
        type=existing_file_path,
        required=True,
        help="Path to the mask file",
    )
    parser.add_argument(
        "--outfile",
        "-o",
        type=path_inside_an_existing_dir,
        required=True,
    )

    return parser.parse_args()


def load_salinity_dataset(
    *,
    minimum_longitude: float,
    maximum_longitude: float,
    minimum_latitude: float,
    maximum_latitude: float,
    maximum_depth: float,
) -> xr.Dataset:
    """
    Get average salinity data from Copernicus Marine DataStore service.

    Loads a salinity dataset for a specified geographical and depth range.
    The dataset provides monthly averaged salinity data for the Mediterranean
    Sea over a defined time period (2012-01-01 to 2021-12-31).
    The dataset ID used is 'cmems_mod_med_phy-sal_my_4.2km_P1M-m'.
    The function computes the mean salinity over the time dimension and
    returns the resulting dataset (as a 3D xarray.DataArray).

    Args:
        minimum_longitude: The minimum longitude of the geographical area.
        maximum_longitude: The maximum longitude of the geographical area.
        minimum_latitude: The minimum latitude of the geographical area.
        maximum_latitude: The maximum latitude of the geographical area.
        maximum_depth: The maximum depth (in meters) for retrieving salinity
            data.

    Returns:
        The salinity dataset averaged over the time dimension within the
        specified geographic and depth constraints.
    """
    t0 = "2012-01-01T00:00:00"
    t1 = "2021-12-31T00:00:00"
    dataset = "cmems_mod_med_phy-sal_my_4.2km_P1M-m"

    sal = cm.open_dataset(
        dataset_id=dataset,
        # dataset_version = version,
        variables=["so"],
        minimum_longitude=minimum_longitude,
        maximum_longitude=maximum_longitude,
        minimum_latitude=minimum_latitude,
        maximum_latitude=maximum_latitude,
        start_datetime=t0,
        end_datetime=t1,
        maximum_depth=maximum_depth,
    ).mean(dim="time")

    return sal.load()


def read_sewage_positions(excel_file: Path, meshmask_file: Path):
    LOGGER.info("Reading sewage positions from excel file %s", excel_file)
    df = pd.read_excel(excel_file)

    # remove rows with no domain assigned, here we filter the 50 discharges
    filtered_df = df[df["Dominio"] > 0]
    LOGGER.info("Excel file read")

    LOGGER.debug("Reading meshmask file %s", meshmask_file)
    domain_mask = Mask.from_mer_file(meshmask_file)

    x0, x1 = domain_mask.lon.min() - 0.5, domain_mask.lon.max() + 0.5
    y0, y1 = domain_mask.lat.min() - 0.5, domain_mask.lat.max() + 0.5
    z1 = domain_mask.zlevels.max() + 10.0

    # To avoid a warning when we download the data, we ensure that summing 0.5
    # does not move us outside the Mediterranean domain. So we add the maximum
    # number between 0.5 and the one that brings us to the edge of the domain
    MAX_MED_LATITUDE = 45.97916793823242
    y1 = max(min(y1, MAX_MED_LATITUDE), domain_mask.lat.max())

    LOGGER.info("Downloading salinity data from MDS")
    sal = load_salinity_dataset(
        minimum_longitude=x0,
        maximum_longitude=x1,
        minimum_latitude=y0,
        maximum_latitude=y1,
        maximum_depth=z1,
    )
    LOGGER.info("Download completed")

    sal["tmask"] = ("depth", "latitude", "longitude"), ~np.isnan(sal.so.values)
    cms_mask = Mask.from_xarray(sal)
    LOGGER.debug("CMS Mask built")

    bathy = cms_mask.bathymetry_in_cells()

    # Select only the interesting columns
    filtered_df = filtered_df[
        [
            "Codice_Scarico",
            "Lat",
            "Long",
            "Carico_Ingresso_AE",
            "Nome_scarico",
            "Nome_impianto",
        ]
    ]
    LOGGER.debug("Columns selected")

    # Here we define the function that we will apply to each row of
    # filtered_df dataset.

    def get_salinity_stats(row: pd.Series):
        """
        Get salinity average for a specific location.

        Retrieves salinity average for a given data row containing
        geographical and discharge-related information from a dataset.
        The function checks if the given coordinates fall within the
        computational model's domain, retrieves bathymetric and
        salinity data, and validates the position. If the position is valid,
        the function extracts and returns depth and average salinity;
        otherwise, it returns a discarded status.

        Args:
            row: A pandas Series containing the following information:
                - "Long" (float): Longitude of the sampled location.
                - "Lat" (float): Latitude of the sampled location.
                - "Nome_scarico" (str): Name of the discharge point.

        Returns:
            pd.Series: A pandas Series containing:
                - "keep" (bool): Whether the data row is valid based on
                        domain and position checks.
                - "CMS_depth" (float): Depth value from the computational
                        model system (CMS), or 0 if discarded.
                - "CMS_avgS" (float): Average salinity value from the CMS,
                        or 0 if discarded.
        """
        longitude = row["Long"]
        latitude = row["Lat"]
        name_scol = row["Nome_scarico"]

        # If this point is discarded, return the following discarded value
        discarded_point_value = pd.Series(
            {"keep": False, "CMS_depth": 0, "CMS_avgS": 0}
        )

        # If the point is outside the domain, return the discarded value
        # (but we still have to check if it is on land or water)
        if not domain_mask.is_inside_domain(lon=longitude, lat=latitude):
            return discarded_point_value

        # Find the coordinate of the closest water point on the
        # CopernicusMarine grid
        i, j = cms_mask.convert_lon_lat_wetpoint_indices(
            lon=longitude, lat=latitude, max_radius=2
        )

        # k is the number of water cells that there are on the position (j, i).
        # If k == -1, this means that there are no water cells in the radius
        # specified by the convert_lon_lat_wetpoint_indices and therefore we
        # discard the point.
        k = bathy[j, i] - 1  # bathymetry from CMS
        if k == -1:
            LOGGER.warning(
                f"Bad position for {name_scol} (lon = {longitude}, lat = "
                f"{latitude}) accordingly to CMS bathymetry, skipping"
            )
            return discarded_point_value

        # Here instead we check the indices for our domain
        j_d, i_d = domain_mask.convert_lat_lon_to_indices(
            lon=longitude, lat=latitude
        )
        if not domain_mask.mask[0, j_d, i_d]:
            LOGGER.warning(
                f"Position of {name_scol} is on land according to the mask"
            )
            return discarded_point_value

        depth = sal.depth[k].values.item()
        salinity = sal.so[k, j, i].values.item()
        return pd.Series(
            {"keep": True, "CMS_depth": depth, "CMS_avgS": salinity}
        )

    # We apply the function to each row of filtered_df and store the results in
    # three new columns: keep, CMS_depth, and CMS_avgS.
    filtered_df[["keep", "CMS_depth", "CMS_avgS"]] = filtered_df.apply(
        get_salinity_stats, axis=1
    )

    # Now we remove all the points where keep is False, and then we drop the
    # "keep" column
    filtered_df = filtered_df[filtered_df["keep"]].drop(columns="keep")

    # We add a new column with the dilution factor
    filtered_df[["Dilution_factor"]] = 1.0

    return filtered_df


def main():
    args = argument()

    excel_file = args.inputfile
    meshmask_file = args.maskfile
    domain_name = args.domain
    outfile = args.outfile

    sewage_df = read_sewage_positions(excel_file, meshmask_file)

    # Convert the selected rows to a list of dictionaries
    sewage_df_rows = sewage_df.to_dict(orient="records")

    # Create the final JSON structure
    output_data = {
        "file_name_origin": excel_file.name,
        "domain_name": domain_name,
        "n_points": len(sewage_df_rows),
        "discharge_points": sewage_df_rows,
    }

    # Write the final JSON structure to a JSON file
    with outfile.open("w") as json_file:
        json.dump(output_data, json_file, indent=4)

    LOGGER.info(f"JSON file created: {outfile}")


if __name__ == "__main__":
    main()
