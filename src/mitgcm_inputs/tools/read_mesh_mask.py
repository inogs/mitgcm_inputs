from pathlib import Path

import xarray as xr
from bitsea.commons.grid import RegularGrid
from bitsea.commons.mask import Mask


def read_mesh_mask(mask_file_path: Path):
    with xr.open_dataset(mask_file_path) as ds:
        if "mer_mesh_mask_version" not in ds.attrs:
            return Mask.from_file(mask_file_path)
        latitude = ds.latitude.values
        longitude = ds.longitude.values
        depth = ds.depth.values
        mesh_mask = ds.tmask.values == 1
    grid = RegularGrid(lat=latitude, lon=longitude)
    mask = Mask(grid, depth, mesh_mask)
    return mask
