import argparse
import json
import logging
from collections import OrderedDict
from collections.abc import Iterable
from collections.abc import Sequence
from itertools import product as cart_prod
from pathlib import Path
from typing import Literal

import dask.array as da
import numpy as np
import xarray as xr
from bitsea.commons.mask import Mask
from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import existing_file_path
from numpy.typing import DTypeLike
from ogs_riverger.read_config import RiverConfig

from mitgcm_inputs.tools.read_mesh_mask import read_mesh_mask


LOGGER = logging.getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Generates RBCs files
    3D Files:
    - bottom_sources_S_mask.bin
           mask with 1.0 at sewage point sources locations, 0.0 elsewhere
    - bottom_sources_S_relaxation.bin MDS
           clim salinity at sewage point sources locations, 0.0 elsewhere
    - SST_mask.bin   1.0 at surface, 0.0 elsewhere

    2D Files with zero values except at point sources locations:
    - conc01_bottom_fluxes.bin
    - conc02_bottom_fluxes.bin
    - ...

    A check file in netcdf format:
    - check_fluxes.nc

    A text file RBCS_names.txt to append on mer config GSN_domain.json:
    ----------------------------------------
    ,
    "tracers": [
        {"id": 1, "name": "Scarico_Trieste_Servola"},
        {"id": 2, "name": "Scarico_Grado"},
        {"id": 3, "name": "Scarico_Depuratore_Di_Staran"},
        {"id": 4, "name": "Timavo"},
        {"id": 5, "name": "Isonzo"}
   ]
   ----------------------------------------

    Args:
        --sewage/-s : Json file output of scarichi_json_gen.py
        --river/-r  : Json file got from bathytools github repository
        --domdir    : Directory of the domain where
                        - MIT_static.nc is expected
                        - Rivers_positions.json is expected
                        - output files are dumped

    """
    )
    parser.add_argument(
        "--sewage",
        "-s",
        type=existing_file_path,
        required=True,
        help="Json file output of scarichi_json_gen.py",
    )
    parser.add_argument(
        "--river",
        "-r",
        type=existing_file_path,
        required=True,
        help="Json file output got from bathytools github repository",
    )

    parser.add_argument(
        "--domdir",
        type=existing_dir_path,
        required=True,
    )
    return parser.parse_args()


def get_spatial_description_from_mit_static(input_file: Path) -> xr.Dataset:
    """
    Read domain information from MITgcm_static.nc

    A "spatial description" is a dataset containing the following variables:
    - depth: 1D variable with the depth of the cells
    - latitude: 1D variable with the latitude of the cells
    - longitude: 1D variable with the longitude of the cells
    - tmask: 3D variable of booleans: True for water, False for land
    - bathymetry: 2D variable (latitude, longitude) with the bathymetry
        of the center of the cells (in meters, downwards positive, it uses
        1e20 as a fill value on the land cells)

    This function constructs the spatial description starting from the
    MITgcm_static.nc file.

    Args:
        input_file: The path of the MITgcm_static.nc file.

    Returns:
        A spatial description of our domain
    """
    with xr.open_dataset(input_file) as ds:
        spatial_description = xr.Dataset(
            data_vars={
                "bathymetry": (("latitude", "longitude"), ds.Depth.values),
            },
            coords={
                "latitude": ds.YC.values,
                "longitude": ds.XC.values,
                "depth": ds.Z.values,
            },
        )

    # Now we compute the tmask (True for water, False for land). First, we
    # find the index of the cell that is closest to the current bathymetry
    bathy_diff = np.abs(
        spatial_description.bathymetry - spatial_description.depth
    )
    # noinspection PyArgumentList
    bottom_index = bathy_diff.argmin(dim="depth")

    # We keep only the indices that are smaller than bottom_index. First, we
    # create a dummy array that simply stores the level number
    level_number = xr.DataArray(
        np.arange(spatial_description.sizes["depth"]),
        dims=("depth",),
    )

    # We also explicitly define the order of the dimensions, otherwise depth
    # will be the last one
    # noinspection PyTypeChecker
    tmask = (level_number <= bottom_index).transpose(
        "depth", "latitude", "longitude"
    )

    # bottom_index is always >= 0. This means that, at the moment, the first
    # level is made only of water cells. We fix this with the following line
    tmask.values[:, spatial_description.bathymetry.values == 0.0] = False

    spatial_description["tmask"] = tmask

    # We fix 1e20 as the fill value for the cells that are on the land
    spatial_description["bathymetry"] = xr.where(
        spatial_description.tmask.isel(depth=0),
        spatial_description.bathymetry,
        1e20,
    )

    return spatial_description


def get_spatial_description_from_meshmask(input_file: Path):
    """
    Read domain information from meshmask.nc

    A "spatial description" is a dataset containing the following variables:
    - depth: 1D variable with the depth of the cells
    - latitude: 1D variable with the latitude of the cells
    - longitude: 1D variable with the longitude of the cells
    - tmask: 3D variable of booleans: True for water, False for land
    - bathymetry: 2D variable (latitude, longitude) with the bathymetry
        of the center of the cells (in meters, downwards positive, it uses
        1e20 as a fill value on the land cells)

    This function constructs the spatial description starting from the
    meshmask.nc file.

    Args:
        input_file: The path of the meshmask.nc file.

    Returns:
        A spatial description of our domain
    """
    mask = read_mesh_mask(input_file)

    spatial_description = mask.to_xarray()

    bathy = xr.DataArray(
        mask.bathymetry(),
        dims=("latitude", "longitude"),
    )
    spatial_description["bathymetry"] = bathy

    return spatial_description


def open_point_sources(
    json_file,
):
    """
    Given the json of a specific domain's point sources,
    reads the data within.
    """

    with open(json_file) as jfile:
        jdata = json.load(jfile)
    return jdata["discharge_points"]


def _allocate_dataarray(
    spatial_description: xr.Dataset,
    dtype: DTypeLike = np.float32,
    dim: Literal[2, 3] = 3,
) -> xr.DataArray:
    """
    Allocate a 2D or 3D DataArray that is compatible with the current domain.

    This function reads the depth, latitude, and longitude information from the
    spatial_description and allocates a DataArray with the appropriate
    dimensions and dtype.

    Args:
        spatial_description: xr.Dataset containing depth, latitude, and
            longitude information. Usually produced by
            `get_spatial_description_from_mit_static` or
            `get_spatial_description_from_meshmask`.
        dtype: The data type of the DataArray.
        dim: 2 if the array is 2D (and its dimensions will be latitude and
            longitude), 3 if the array is 3D (and in this case it will also
            include the depth dimension).

    Returns: A DataArray with the appropriate dimensions and dtype filled with
        zeros.

    """
    if dim < 2 or dim > 3:
        raise ValueError("Only 2D and 3D data are supported.")

    var_3d_sizes = (
        spatial_description.sizes["depth"],
        spatial_description.sizes["latitude"],
        spatial_description.sizes["longitude"],
    )
    var_3d_dims = ("depth", "latitude", "longitude")

    if dim == 2:
        var_sizes = var_3d_sizes[1:]
        var_dims = var_3d_dims[1:]
    else:
        var_sizes = var_3d_sizes
        var_dims = var_3d_dims

    coords = {
        "latitude": spatial_description.latitude,
        "longitude": spatial_description.longitude,
    }

    if dim == 3:
        coords["depth"] = spatial_description.depth

    return xr.DataArray(
        np.zeros(var_sizes, dtype=dtype), dims=var_dims, coords=coords
    )


def get_opensea_swg_buoyant_plume(
    spatial_description: xr.Dataset,
    sewage_points: Iterable[dict],
    fixed_conc: float | None,
    water_freshener: float = 0.5,
) -> xr.Dataset:
    """

    Args:
        spatial_description: The spatial dataset providing depth, latitude,
            and longitude information necessary for constructing the salinity
            relaxation fields.
        sewage_points: Collection of sewage data points where each point is a
            dictionary containing information about longitude, latitude,
            sewage parameters, and other relevant data.
        fixed_conc: If this entry is a float, all the sewage concentrations are
            fixed to this value. Otherwise, if it is `None`, the concentrations
            are computed from the data inside the sewage_points dict.
        water_freshener (float, optional): A scalar used to reduce the salinity
            relaxation value by a constant amount to model the effect of fresh
            water introduction. Defaults to 0.5.

    Returns:
        A dataset containing data variables for relaxed salinity, salinity
        mask, sewage concentration indices, and concentration values.
    """
    relax_salt = _allocate_dataarray(spatial_description, dtype=np.float32)
    salinity_mask = xr.zeros_like(relax_salt, dtype=np.float32)
    salinity_mask.values[:] = 1.0

    sst_mask = xr.zeros_like(relax_salt, dtype=np.float32)
    sst_mask[{"depth": 0}] = 1.0

    conc_indices = xr.zeros_like(relax_salt.isel(depth=0), dtype=int)
    conc_values = xr.zeros_like(relax_salt.isel(depth=0), dtype=np.float32)

    mask = Mask.from_xarray(spatial_description)
    bottom_level_index = mask.bathymetry_in_cells() - 1

    for ns, jd in enumerate(sewage_points, start=1):
        lon = jd["Long"]
        lat = jd["Lat"]

        i, j = mask.convert_lat_lon_to_indices(lon=lon, lat=lat)
        k = bottom_level_index[i, j]

        if k == -1:
            raise ValueError(
                f"The following conc is on the land: {jd['Nome_impianto']}"
            )

        relaxed_salinity_value = jd["CMS_avgS"] - water_freshener

        if fixed_conc is None:
            c = jd["Carico_Ingresso_AE"] * jd["Dilution_factor"]
        elif isinstance(fixed_conc, float):
            c = fixed_conc
        else:
            raise ValueError(
                f"fixed_conc must be 'None' or a float, got {fixed_conc}"
            )

        relax_salt[k, i, j] = relaxed_salinity_value
        salinity_mask[k, i, j] = 0.0
        conc_indices[i, j] = ns
        conc_values[i, j] = c

    return xr.Dataset(
        data_vars={
            "relaxed_salinity": relax_salt,
            "salinity_mask": salinity_mask,
            "sst_mask": sst_mask,
            "conc_values": conc_values,
            "conc_indices": conc_indices,
        }
    )


def get_river_swg_plume(
    *,
    spatial_description: xr.Dataset,
    rivers: Iterable[dict],
    uniform_concentration=1000.0,
    fixed_conc: float | None = None,
) -> xr.Dataset:
    conc_indices = _allocate_dataarray(spatial_description, dtype=int, dim=2)
    conc_values = xr.zeros_like(conc_indices, dtype=np.float32)

    for r_num, river in enumerate(rivers, start=1):
        source_cells = cart_prod(
            river["latitude_indices"], river["longitude_indices"]
        )

        for j, i in source_cells:
            # shift indices of one cell downstream, according to the side of
            # the river, to avoid overwriting of open boundary conditions
            # (which are all zeros by construction for the E. coli
            # concentrations)
            if river["side"] == "E":
                i -= 1
            elif river["side"] == "W":
                i += 1
            elif river["side"] == "S":
                j += 1
            elif river["side"] == "N":
                j -= 1

            if fixed_conc is None:
                c = uniform_concentration
            elif isinstance(fixed_conc, float):
                c = fixed_conc
            else:
                raise ValueError(
                    f"fixed_conc must be 'None' or a float, got {fixed_conc}"
                )
            conc_indices[j, i] = r_num
            conc_values[j, i] = c

    return xr.Dataset(
        {
            "conc_values": conc_values,
            "conc_indices": conc_indices,
        }
    )


def write_binary_files(conc_and_relax: xr.Dataset, out_dir: Path):
    """
    Args:
        relax_sal : ndarray [depth,lat,lon] values to relax to
        S_mask    : ndarray [depth,lat,lon] of 0.0 1.0 (float)
        sst_mask  : ndarray [depth,lat,lon] of 0.0 1.0 (float)
        concs     : list of xr DataArrays [lat,lon] with sources concentrations

    Given the relaxation salinity, salinity mask and concentration
    arrays, writes them to MITgcm-appropriate binary files:
            - bottom_sources_S_mask.bin
            - bottom_sources_S_relaxation.bin
            - SST_mask.bin
            - conc01_bottom_fluxes.bin
            - conc02_bottom_fluxes.bin
            - ...
            - check_fluxes.nc

    """
    filename = out_dir / "bottom_sources_S_relaxation.bin"
    conc_and_relax.relaxed_salinity.values.astype("f4").tofile(filename)
    LOGGER.info("File %s has been written", filename)

    filename = out_dir / "bottom_sources_S_mask.bin"
    conc_and_relax.salinity_mask.values.astype("f4").tofile(filename)
    LOGGER.info("File %s has been written", filename)

    filename = out_dir / "SST_mask.bin"
    conc_and_relax.sst_mask.values.astype("f4").tofile(filename)
    LOGGER.info("File %s has been written", filename)

    for data_var in conc_and_relax.data_vars:
        if not data_var.startswith("conc"):
            continue

        filename = out_dir / f"{data_var}_bottom_fluxes.bin"
        conc_and_relax[data_var].values.astype("f4").tofile(filename)
        LOGGER.info("File %s has been written", filename)

    filename = out_dir / "check_fluxes.nc"
    conc_and_relax.to_netcdf(filename)


def build_domain_config_string(sources: list[tuple[int, str, str]]) -> str:
    # text to insert in GSN_domain.json file,
    # just after
    #     "model_parameters": {
    #     "OBCSsponge_N": ".FALSE.",
    #     "OBCSsponge_S": ".TRUE.",
    #     "OBCSsponge_E": ".FALSE.",
    #     "OBCSsponge_W": ".TRUE.",

    #     "useOBCSbalance": ".TRUE.",
    #     "OBCS_balanceFacN": "0.",
    #     "OBCS_balanceFacS": "1.",
    #     "OBCS_balanceFacE": "0.",
    #     "OBCS_balanceFacW": "1."
    # }
    output = ',\n    "tracers": [\n'
    for i, s in enumerate(sources):
        comma = "," if i < len(sources) - 1 else ""
        output += "        " + json.dumps(s) + comma + "\n"
    output += "   ]"
    return output


def build_conc_and_relax_variables(
    spatial_description: xr.Dataset,
    sewage_points: Sequence[dict],
    rivers: Sequence[dict] | None = None,
) -> tuple[list, xr.Dataset]:
    if rivers is None:
        rivers = []

    swg_buoyant_plume = get_opensea_swg_buoyant_plume(
        spatial_description=spatial_description,
        sewage_points=sewage_points,
        fixed_conc=1.0,
    )
    n_concs = int(swg_buoyant_plume.conc_indices.max())
    LOGGER.info("Produced %i concentrations from the sewage sources", n_concs)

    rivers_conc = get_river_swg_plume(
        spatial_description=spatial_description,
        rivers=rivers,
        fixed_conc=1.0,
    )

    rivers_conc_indices = xr.where(
        rivers_conc.conc_indices != 0, rivers_conc.conc_indices + n_concs, 0
    )
    LOGGER.info(
        "Produced %i concentrations from the rivers",
        rivers_conc.conc_indices.max(),
    )

    # Create a list of the sources, together with their name and their kind
    sources = []
    for s in sewage_points:
        sources.append(
            OrderedDict(
                [
                    ("id", len(sources) + 1),
                    ("name", s["Nome_impianto"]),
                    ("kind", "sewage"),
                ]
            )
        )
    for r in rivers:
        sources.append(
            OrderedDict(
                [
                    ("id", len(sources) + 1),
                    ("name", r["name"]),
                    ("kind", "river"),
                ]
            )
        )

    total_conc_indices = swg_buoyant_plume.conc_indices + rivers_conc_indices
    total_conc_values = swg_buoyant_plume.conc_values + rivers_conc.conc_values

    salinity_and_concs = xr.Dataset(
        data_vars={
            "relaxed_salinity": swg_buoyant_plume.relaxed_salinity,
            "salinity_mask": swg_buoyant_plume.salinity_mask,
            "sst_mask": swg_buoyant_plume.sst_mask,
        },
        coords={
            "latitude": spatial_description.latitude,
            "longitude": spatial_description.longitude,
            "depth": spatial_description.depth,
        },
    )

    # We transform conc_indices and conc_values into different concs variables
    n_concs = int(total_conc_indices.max())

    for i in range(1, n_concs + 1):
        current_conc_values = da.where(
            total_conc_indices == i, total_conc_values, 0.0
        )
        current_conc = xr.DataArray(
            current_conc_values,
            dims=("latitude", "longitude"),
            coords={
                "latitude": salinity_and_concs.latitude,
                "longitude": salinity_and_concs.longitude,
            },
        )
        salinity_and_concs["conc" + f"{i:02}"] = current_conc

    return sources, salinity_and_concs


def main():
    args = argument()

    input_file = args.domdir / "MIT_static.nc"
    rivers_positions_file = args.domdir / "rivers_positions.json"

    river_positions = json.loads(
        rivers_positions_file.read_text(encoding="utf-8")
    )

    if args.river is not None:
        rivers_config = RiverConfig.from_json(args.river)
    else:
        rivers_config = RiverConfig(root=OrderedDict())

    def has_tracer(river):
        river_id = river["id"]
        return len(rivers_config.root[river_id].concentrations) > 0

    # Keep only rivers with tracers
    rivers = tuple(r for r in river_positions if has_tracer(r))

    # Read the description of the domain from MITgcm static file
    spatial_description = get_spatial_description_from_mit_static(input_file)

    # Read the JSON with the sewage sources
    sewage_points = open_point_sources(args.sewage)

    sources, conc_and_relax = build_conc_and_relax_variables(
        spatial_description=spatial_description,
        sewage_points=sewage_points,
        rivers=rivers,
    )

    config_string = build_domain_config_string(sources)

    point_source_output = args.domdir / "RBCS_names.txt"
    point_source_output.write_text(config_string)

    write_binary_files(conc_and_relax, args.domdir)


if __name__ == "__main__":
    main()
