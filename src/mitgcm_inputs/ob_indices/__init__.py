import argparse
import json
import logging

import xarray as xr
from bitsea.commons.mask import MaskWithRivers
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import path_inside_an_existing_dir

from mitgcm_inputs.ob_indices.ob_indices import generate_ob_indices
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
        "information is not already stored in the meshmask file and if there "
        "are rivers on the domain. If the rivers-positions file is not "
        "provided, we will assume that there are no rivers on the domain "
        "and this file is simply ignored",
    )

    parser.add_argument(
        "-p",
        "--rivers-positions",
        required=False,
        type=existing_file_path,
        default=None,
        help="The json file that describe the position of the mouths of the"
        "rivers. If this file is not provided, we will assume that there are "
        "no rivers on the domain",
    )

    parser.add_argument(
        "-s",
        "--sponge-extent",
        required=False,
        default=0,
        type=int,
        help="The extent of the sponge layer in number of cells",
    )

    parser.add_argument(
        "-o",
        "--ob-indices",
        required=True,
        type=path_inside_an_existing_dir,
        help="The path of the output file where the OB indices will be "
        "stored. If this is not submitted, this output file will not be "
        "written.",
    )

    parser.add_argument(
        "-n",
        "--nudge-indices",
        required=False,
        type=path_inside_an_existing_dir,
        default=None,
        help="The path of the output file where the nudge indices will be "
        "stored",
    )


def main(args: argparse.Namespace) -> int:
    if args.cmd != COMMAND_NAME:
        LOGGER.error(
            "obindices command has been invoked with command: %s", args.cmd
        )
        return 1

    domain_mask = read_mesh_mask(args.mask)

    # On which file do we have the information about the cells occupied by
    # the rivers?
    if args.river_file is not None:
        file_with_river = args.river_file
    else:
        file_with_river = args.mask

    # If there are rivers in this domain
    if args.rivers_positions is not None:
        with xr.open_dataset(file_with_river) as river_ds:
            river_map = river_ds["rivers"].load()

        river_mask = MaskWithRivers(
            grid=domain_mask.grid,
            zlevels=domain_mask.zlevels,
            mask_array=domain_mask[:],
            river_positions=river_map,
        )
        rivers_positions = json.loads(args.rivers_positions.read_text())
        mesh_mask = river_mask
    else:
        rivers_positions = []
        mesh_mask = domain_mask

    ob_indices, ob_sponge = generate_ob_indices(
        river_mask=mesh_mask,
        sponge_extent=args.sponge_extent,
        rivers_positions=rivers_positions,
    )

    args.ob_indices.write_text(ob_indices)
    if args.nudge_indices is not None:
        args.nudge_indices.write_text(ob_sponge)

    return 0
