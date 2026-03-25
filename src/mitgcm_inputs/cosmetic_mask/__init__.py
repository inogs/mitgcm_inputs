import argparse
import logging

import numpy as np
import xarray as xr
from bitsea.commons.mask import MaskWithRivers
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import path_inside_an_existing_dir

from mitgcm_inputs.tools.read_mesh_mask import read_mesh_mask

if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


COMMAND_NAME = "cosmetic_mask"


def sub_arguments(subparser):
    parser = subparser.add_parser(
        COMMAND_NAME,
        help="Generate a new meshmask file that has no rivers and that can be "
        "used to produce plots without the rivers stems",
    )

    parser.add_argument(
        "-m",
        "--mask",
        required=True,
        type=existing_file_path,
        help="The meshmask file that defines the domain",
    )

    parser.add_argument(
        "-r",
        "--river-file",
        required=False,
        default=None,
        type=existing_file_path,
        help="An auxiliary file that contains a map that associates each "
        "river cell with its river id. This file is useful only if this "
        "information is not already stored in the meshmask file and if there "
        "are rivers on the domain. If the rivers-positions file is not "
        "provided, we will assume that the information about the rivers is "
        'stored inside the meshmask file (in a variable named "rivers")',
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=path_inside_an_existing_dir,
        help="Path of the NetCDF output file to write",
    )

    parser.add_argument(
        "--mer",
        action="store_true",
        help="""
        If this flag is set, this script will generate the meshmask in the
        format used by the MER project. This format is CF-1.4 compliant.
        """,
    )


def main(args: argparse.Namespace) -> int:
    if args.cmd != COMMAND_NAME:
        LOGGER.error(
            "%s command has been invoked with command: %s",
            COMMAND_NAME,
            args.cmd,
        )
        return 1

    domain_mask = read_mesh_mask(args.mask)

    # On which file do we have the information about the cells occupied by
    # the rivers?
    if args.river_file is not None:
        file_with_river = args.river_file
    else:
        file_with_river = args.mask

    with xr.open_dataset(file_with_river) as river_ds:
        river_map = river_ds["rivers"].load()

    river_mask = MaskWithRivers(
        grid=domain_mask.grid,
        zlevels=domain_mask.zlevels,
        mask_array=domain_mask[:],
        river_positions=river_map,
    )

    if args.mer:
        mesh_mask_coords = {
            "depth": (
                ("depth",),
                river_mask.zlevels,
                {
                    "units": "m",
                    "positive": "down",
                    "standard_name": "depth",
                    "axis": "Z",
                },
            ),
            "latitude": (
                ("latitude",),
                river_mask.lat,
                {
                    "units": "degrees_north",
                    "standard_name": "latitude",
                    "axis": "Y",
                },
            ),
            "longitude": (
                ("longitude",),
                river_mask.lon,
                {
                    "units": "degrees_east",
                    "standard_name": "longitude",
                    "axis": "X",
                },
            ),
        }
        e1t = xr.DataArray(data=river_mask.e1t, dims=("latitude", "longitude"))
        e2t = xr.DataArray(data=river_mask.e2t, dims=("latitude", "longitude"))
        e3t = xr.DataArray(
            data=river_mask.e3t, dims=("depth", "latitude", "longitude")
        )
        del_z = xr.DataArray(
            data=river_mask.zlevels,
            dims=("depth",),
        )
        tmask = xr.DataArray(
            data=np.asarray(river_mask, dtype=np.int8),
            dims=("depth", "latitude", "longitude"),
        )

        mask_array = xr.Dataset(
            data_vars={
                "e1t": e1t,
                "e2t": e2t,
                "e3t": e3t,
                "delZ": del_z,
                "tmask": tmask,
            },
            coords=mesh_mask_coords,
            attrs={"mer_mesh_mask_version": "1.0"},
        )
        mask_array.to_netcdf(args.output)
    else:
        river_mask.save_as_netcdf(args.output)

    return 0
