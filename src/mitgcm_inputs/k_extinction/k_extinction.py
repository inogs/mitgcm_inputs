import logging
from os import PathLike
from pathlib import Path

import numpy as np
import xarray as xr
from bitsea.commons.mask import Mask
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay


DATA_FILE = Path(__file__).resolve().parent / "data" / "kExt_climatology.nc"


LOGGER = logging.getLogger(__name__)


def compute_k_extinction(
    meshmask: Mask, data_file: PathLike = DATA_FILE
) -> xr.Dataset:
    """
    Compute the K extinction coefficients by interpolating climatological
    values precomputed on a coarse but larger grid. The function produces a
    map of K extinction coefficients for each day of the year.

    Parameters:
    meshmask (Mask): The meshmask defining the domain.
    data_file (PathLike): Path to the climatological K extinction data file.

    Returns:
    xr.Dataset: Dataset containing K extinction coefficients
    """
    LOGGER.info("Computing K extinction coefficients")

    LOGGER.debug("Meshmask shape: %s", meshmask.shape)

    LOGGER.debug("Reading data file: %s", data_file)
    climatological_data = xr.load_dataset(data_file)

    data_lon_min = np.min(climatological_data.longitude.values)
    data_lon_max = np.max(climatological_data.longitude.values)
    data_lat_min = np.min(climatological_data.latitude.values)
    data_lat_max = np.max(climatological_data.latitude.values)
    LOGGER.debug(
        "Data file loaded; grid lat: (%s, %s), grid lon: (%s, %s)",
        data_lat_min,
        data_lat_max,
        data_lon_min,
        data_lon_max,
    )

    mesh_lon = meshmask.xlevels[~meshmask[0]]
    mesh_lat = meshmask.ylevels[~meshmask[0]]

    if (min_mesh_lon := np.min(mesh_lon)) < data_lon_min:
        raise ValueError(
            f"The minimum longitude of the meshmask ({min_mesh_lon:.2f}) is "
            f"outside of the range of the climatological data "
            f"(from {data_lon_min:.2f} to {data_lon_max:.2f})"
        )
    if (max_mesh_lon := np.max(mesh_lon)) > data_lon_max:
        raise ValueError(
            f"The maximum longitude of the meshmask ({max_mesh_lon:.2f}) is "
            f"outside of the range of the climatological data "
            f"(from {data_lon_min:.2f} to {data_lon_max:.2f})"
        )
    if (min_mesh_lat := np.min(mesh_lat)) < data_lat_min:
        raise ValueError(
            f"The minimum latitude of the meshmask ({min_mesh_lat:.2f}) is "
            f"outside of the range of the climatological data "
            f"(from {data_lat_min:.2f} to {data_lat_max:.2f})"
        )
    if (max_mesh_lat := np.max(mesh_lat)) > data_lat_max:
        raise ValueError(
            f"The maximum longitude of the meshmask ({max_mesh_lat:.2f}) is "
            f"outside of the range of the climatological data "
            f"(from {data_lat_min:.2f} to {data_lat_max:.2f})"
        )

    climatological_kext = climatological_data["kExt"].to_masked_array()
    water_data_cells = ~np.ma.getmaskarray(climatological_kext[0, :, :])

    grid_point_lon = climatological_data.longitude.values[water_data_cells]
    grid_point_lat = climatological_data.latitude.values[water_data_cells]

    n_days = climatological_kext.shape[0]
    output_shape = (n_days,) + meshmask.shape[1:]
    LOGGER.debug("Output shape: %s", output_shape)
    output = np.full(output_shape, np.nan, dtype=climatological_kext.dtype)

    triangulation = Delaunay(np.column_stack((grid_point_lon, grid_point_lat)))

    for d in range(n_days):
        LOGGER.debug("Processing day %s of %s", d + 1, n_days)
        # It would be really nice if we could compute the interpolation once
        # and then just change its values using the `set_values` method.
        # Unfortunately, at the time being, the pull request
        # https://github.com/scipy/scipy/pull/18376/commits/3b0bf8bc2b010d5d4f7bde5da75df6de8403d3c8
        # is not released yet. So let's rebuild the same interpolator for 365
        # times!
        interpolator = LinearNDInterpolator(
            triangulation, climatological_kext.data[d, water_data_cells]
        )
        day_data = interpolator(mesh_lon, mesh_lat)

        output[d, ~meshmask[0]] = day_data

    LOGGER.debug("All days have been computed")

    return xr.Dataset(
        data_vars={
            "kExt": (["day", "x", "y"], output),
        },
        coords={
            "day": np.arange(1, n_days + 1),
            "latitude": (("x", "y"), meshmask.ylevels),
            "longitude": (("x", "y"), meshmask.xlevels),
        },
    )
