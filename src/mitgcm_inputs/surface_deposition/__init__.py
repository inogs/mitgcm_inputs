import argparse
import logging

import numpy as np
from bitsea.utilities.argparse_types import dir_to_be_created_if_not_exists
from bitsea.utilities.argparse_types import existing_file_path

from mitgcm_inputs.surface_deposition.surface_deposition import (
    compute_surface_deposition,
)
from mitgcm_inputs.tools.read_mesh_mask import read_mesh_mask
from mitgcm_inputs.tools.save_dataset import save_dataset


if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


COMMAND_NAME = "surface_deposition"


def sub_arguments(subparser):
    parser = subparser.add_parser(
        COMMAND_NAME,
        help="Generate the surface deposition climatological data",
    )

    parser.add_argument(
        "-m",
        "--mask",
        required=True,
        type=existing_file_path,
        help="The meshmask file that defines the domain",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=dir_to_be_created_if_not_exists,
        help="Path of the directory where the output files will be written",
    )

    parser.add_argument(
        "-z",
        "--output-name",
        required=False,
        default="${VAR}_surface_fluxes.bin",
        type=str,
        help="Name of the output files",
    )


def main(args: argparse.Namespace) -> int:
    if args.cmd != COMMAND_NAME:
        LOGGER.error(
            "Surface deposition command has been invoked with command: %s",
            args.cmd,
        )
        return 1

    mask = read_mesh_mask(args.mask)

    output_data = compute_surface_deposition(mask)
    output_dtype = np.float32

    save_dataset(
        output_data,
        args.output,
        args.output_name,
        output_dtype,
    )

    LOGGER.info("Execution completed!")

    return 0
