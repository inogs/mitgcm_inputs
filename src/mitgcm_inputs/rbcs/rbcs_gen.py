import argparse
import json
import logging
from collections import OrderedDict
from collections.abc import Iterable
from collections.abc import Sequence
from itertools import product as cart_prod
from pathlib import Path
from typing import Any
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
    Generates a Dataset of buoyant plume data for sewage discharge points.

    This function creates a dataset based on spatial descriptions and a set of
    sewage points. The output dataset has 5 variables:
    - `relaxed_salinity`
    - `salinity_mask`
    - `sst_mask`
    - `conc_values`
    - `conc_indices`

    The first two, `relaxed_salinity` and `salinity_mask`, are 3D arrays
    that represent where and how the salinity must be relaxed. `salinity_mask`
    is 1 everywhere except at the point where the salinity is relaxed.
    `relaxed_salinity` contains the values of the salinity at the relaxed
    points (and it is zero everywhere else).

    The third, `sst_mask`, is a 3D array that is 0 everywhere except at the
    surface.

    Finally, `conc_values` and `conc_indices` are 2D arrays that contain the
    information about the position and the values of the tracers (and their
    concentration). If there are n tracers, `conc_indices` contains integer
    values between 0 and n. The concentration values of the `i`-th tracers are
    0 everywhere besides on the cells where `conc_indices` is equal to `i`.
    On these cells, the concentration values are the values of `conc_values`.
    In this way, we can store all the information about the concentration of
    many tracers using just two variables.

    Args:
        spatial_description (xr.Dataset): The spatial data describing the
            area of interest, which includes grid coordinates and bathymetry
            information. Look at the documentation of the functions
            `get_spatial_description_from_mit_static` or
            `get_spatial_description_from_meshmask` to have an accurate
            description of a spatial description.
        sewage_points: An iterable of dictionaries, where each dictionary
            describes properties of individual sewage discharge points,
            including geographical coordinates and sewage properties.
        fixed_conc: An optional fixed value for sewage point concentration.
            If not provided, concentrations are computed dynamically based on
            input properties.
        water_freshener (float): The adjustment value for fresh water in
            salinity calculations, defaulting to 0.5.

    Returns:
        A dataset containing the computed buoyant plume data for the sewage
        discharge points.

    Raises:
        ValueError: If any sewage discharge point is on land or if the
        `fixed_conc` value is invalid.
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

        j, i = mask.convert_lat_lon_to_indices(lon=lon, lat=lat)
        k = bottom_level_index[j, i]

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

        relax_salt[k, j, i] = relaxed_salinity_value
        salinity_mask[k, j, i] = 0.0
        conc_indices[j, i] = ns
        conc_values[j, i] = c

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
    """
    Generate a dataset with the position of the pollutant from the rivers.

    This function does for the rivers mouths what
    `get_opensea_swg_buoyant_plume` does for the sewage points. It iterates
    over the rivers and generates a dataset containing the position of the
    E.coli concentration sources and their concentration values.
    The output dataset has two variables:
      - "conc_values"
      - "conc_indices"

    These two variables have the same meaning that they have for the
    `get_opensea_swg_buoyant_plume` function.

    Args:
        spatial_description: A dataset describing the spatial grid in
            which the river plume will be projected.
        rivers: A collection of dictionaries where each dictionary
            describes an individual river discharge source. Each dictionary
            must include at least the keys "latitude_indices",
            "longitude_indices", and "side", detailing the source cell indices
            and river outflow orientation.
        uniform_concentration: Default concentration value to be used for all
            river discharge sources when `fixed_conc` is None. Defaults to
            1000.0.
        fixed_conc (optional): Fixed concentration value to override
            `uniform_concentration` for all river discharge sources. Must be a
            float or None. Defaults to None.

    Returns:
        A dataset containing two data arrays:
            - "conc_values": A 2D data array containing the pollutant
              concentration values across the spatial grid.
            - "conc_indices": A 2D data array containing numerical indices
              marking the source regions of pollutant discharge.
    """
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
) -> tuple[list[OrderedDict[str, Any]], xr.Dataset]:
    """
    Builds concentration and relaxation vars based on sewage and river data.

    This function processes sewage points and optional river data to generate
    concentration indices and values for both, combining them into a
    consolidated dataset. It also produces a list of source metadata for the
    inputs.

    This function calls `get_opensea_swg_buoyant_plume` and
    `get_river_swg_plume` and then merges the results into a single dataset.

    Args:
        spatial_description: Dataset containing spatial details, including
            latitude, longitude, and depth.
        sewage_points: List of dictionaries representing sewage points.
            Each dictionary should include metadata about a sewage source such
            as its name.
        rivers: Optional list of dictionaries representing river sources. Each
            dictionary should include metadata about a river source, such as
            its name. Defaults to None.

    Returns:
        A tuple where the first element is a list of sources with their
        metadata and the second element is an xarray dataset with relaxation
        and concentration variables.

        The list of sources contains a dictionary for each source. This
        dictionary contains 4 keys:
        - "id": A unique identifier for the source.
        - "name": The name of the source. If the sources are sewage points,
          this is the "Name_impianto" key from the sewage point dictionary.
          Otherwise, it is the name of the river source.
        - "kind": The type of the source. Can be "sewage" or "river".
        - "frequency": The frequency of the source. This describes if and how
          the source changes during time. At the moment, it is always "fixed".

        The second element is an xarray dataset containing `3 + n` variables,
        where `n` is the number of unique concentration indices produced.
        There is one concentration variable for each tracer and three other
        variables:
        - `relaxed_salinity`
        - `salinity_mask`
        - `sst_mask`
        These 3 variables are generated by the `get_opensea_swg_buoyant_plume`
        function.

    """
    if rivers is None:
        rivers = []

    swg_buoyant_plume = get_opensea_swg_buoyant_plume(
        spatial_description=spatial_description,
        sewage_points=sewage_points,
        fixed_conc=1.0,
    )
    n_concs = int(swg_buoyant_plume.conc_indices.max())
    LOGGER.info("Produced %i concentrations from the sewage sources", n_concs)

    # Now we read the concentrations from the rivers
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

    # Create a list of the sources, together with their name and kind
    sources = []
    for s in sewage_points:
        sources.append(
            OrderedDict(
                [
                    ("id", len(sources) + 1),
                    ("name", s["Nome_impianto"]),
                    ("kind", "sewage"),
                    ("frequency", "fixed"),
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
                    ("frequency", "fixed"),
                ]
            )
        )

    total_conc_indices = swg_buoyant_plume.conc_indices + rivers_conc_indices
    total_conc_values = swg_buoyant_plume.conc_values + rivers_conc.conc_values

    # Now we prepare the dataset that is the second part of the output.
    # We introduce the three masks about the salinity and SST and then we copy
    # the concentrations of the tracers from the conc_values and conc_indices
    # variables.
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
