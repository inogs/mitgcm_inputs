import argparse
import logging

import xarray as xr
from bitsea.commons.mask import MaskWithRivers
from bitsea.utilities.argparse_types import existing_file_path

from mitgcm_inputs.tools.read_mesh_mask import read_mesh_mask

if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


COMMAND_NAME = "ob_indices"


def sub_arguments(subparser):
    parser = subparser.add_parser(
        COMMAND_NAME, help="Generate OB indices for rivers and sponge layers"
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
        help="An auxiliary file that contains the a map that associates each "
        "river cell with its river id. This file is useful only if this "
        "information is not already stored in the meshmask file",
    )


def main(args: argparse.Namespace) -> int:
    if args.cmd != COMMAND_NAME:
        LOGGER.error(
            "obindices command has been invoked with command: %s", args.cmd
        )
        return 1

    meshmask = read_mesh_mask(args.mask)

    if args.river_file is not None:
        file_with_river = args.river_file
    else:
        file_with_river = args.mask
    with xr.open_dataset(file_with_river) as river_ds:
        river_map = river_ds["rivers"].load()

    river_mask = MaskWithRivers(
        grid=meshmask.grid,
        zlevels=meshmask.zlevels,
        mask_array=meshmask[:],
        river_positions=river_map,
    )
    print(river_mask)
