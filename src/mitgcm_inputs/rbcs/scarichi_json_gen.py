import argparse
import json
import logging

import copernicusmarine as cm
import numpy as np
import pandas as pd
from bitsea.commons.grid import RegularGrid
from bitsea.commons.mask import Mask
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import path_inside_an_existing_dir

LOGGER = logging.getLogger()


def configure_logger():
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    LOGGER.setLevel(logging.INFO)

    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)

    LOGGER.addHandler(handler)


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Reads seawage info on excel file and for all the domains
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


args = argument()
configure_logger()

excel_file = args.inputfile


def load_datasets(x0, x1, y0, y1, z1):
    t0 = "2012-01-01T00:00:00"
    t1 = "2021-12-31T00:00:00"
    dataset = "cmems_mod_med_phy-sal_my_4.2km_P1M-m"
    var = "so"
    z0 = 1.0182366371154785  # 1 meter depth

    sal = cm.open_dataset(
        dataset_id=dataset,
        # dataset_version = version,
        username="",
        password="",
        variables=[var],
        minimum_longitude=x0,
        maximum_longitude=x1,
        minimum_latitude=y0,
        maximum_latitude=y1,
        start_datetime=t0,
        end_datetime=t1,
        minimum_depth=z0,
        maximum_depth=z1,
    ).mean(dim="time")

    # import matplotlib.pyplot as plt
    # import cmocean as cmo
    # sal.so[0,0,:,:].plot.imshow(); plt.show()

    print("Loaded datasets!")
    return sal


LOGGER.info("Start")
# Read the Excel file into a DataFrame
df = pd.read_excel(excel_file)
filtered_df = df[
    df["Dominio"] > 0
]  # remove rows with no domain assigned, here we filter the 50 discharges

LOGGER.info("Excel read")
namedomain = args.domain
M = Mask.from_mer_file(args.maskfile)
x0, x1 = M.lon.min() - 0.5, M.lon.max() + 0.5
y0, y1 = M.lat.min() - 0.5, M.lat.max() + 0.5
z1 = M.zlevels.max() + 10.0
sal = load_datasets(x0, x1, y0, y1, z1)
sal = sal.load()

LOGGER.info("Datasets loaded")
grid = RegularGrid(lat=sal.latitude.values, lon=sal.longitude.values)
LOGGER.info("Grid built")
CMSmask = Mask(
    grid=grid, zlevels=sal.depth.values, mask_array=~np.isnan(sal.so.values)
)
bathy = CMSmask.bathymetry_in_cells()
LOGGER.info("CMS Mask built")
# Select only the interesting columns
selected_columns = filtered_df[
    [
        "Codice_Scarico",
        "Lat",
        "Long",
        "Carico_Ingresso_AE",
        "Nome_scarico",
        "Nome_impianto",
    ]
]
LOGGER.info("Columns selected")
n_scol = selected_columns.shape[0]
good = np.zeros(n_scol, dtype=bool)
scol_list = []
depths = []
Savg = []
Sstd = []
dilut_fac = 1.0
for iscol in range(n_scol):
    lat_scol = selected_columns.Lat.to_numpy()[iscol]
    lon_scol = selected_columns.Long.to_numpy()[iscol]
    name_scol = selected_columns.Nome_impianto.to_numpy()[iscol]
    if M.is_inside_domain(lon=lon_scol, lat=lat_scol):
        LOGGER.info(
            f"Processing point {iscol}: lon={lon_scol}, lat={lat_scol}"
        )
        i, j = CMSmask.convert_lon_lat_wetpoint_indices(
            lon=lon_scol, lat=lat_scol
        )

        k = bathy[j, i] - 1  # bathymetry from CMS
        if k == -1:
            LOGGER.warning(
                f"Bad position for {name_scol} for CMS bathymetry, skipping"
            )
            continue
        jd, id = M.convert_lat_lon_to_indices(lon=lon_scol, lat=lat_scol)
        if not M.mask[0, jd, id]:
            LOGGER.warning(
                f"Position of {name_scol} is on land according to the mask"
            )
            continue
        depths.append(sal.depth[k].values.item())
        Savg.append(sal.so[k, j, i].values.item())
        good[iscol] = True

selected_columns = selected_columns[good]
selected_columns.insert(5, "CMS_depth", depths, True)
selected_columns.insert(6, "CMS_avgS", Savg, True)
selected_columns.insert(8, "Dilution_factor", dilut_fac, True)

# Convert the selected rows to a list of dictionaries
filtered_rows = selected_columns.to_dict(orient="records")

# Create the final JSON structure
output_data = {
    "file_name_origin": excel_file.name,
    "domain_name": namedomain,
    "n_points": len(filtered_rows),
    "discharge_points": filtered_rows,
}

# Write the final JSON structure to a JSON file
with open(args.outfile, "w") as json_file:
    json.dump(output_data, json_file, indent=4)

print(f"JSON file created: {args.outfile}")
