import argparse
import logging

import numpy as np
from bitsea.commons.mask import Mask
from bitsea.utilities.argparse_types import dir_to_be_created_if_not_exists
from bitsea.utilities.argparse_types import existing_file_path

from mitgcm_inputs.bottom_fluxes.bottom_fluxes import compute_bottom_fluxes
from mitgcm_inputs.tools.save_dataset import save_dataset


if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


COMMAND_NAME = "bottom_fluxes"


def sub_arguments(subparser):
    parser = subparser.add_parser(
        "bottom_fluxes", help="Generate vertical fluxes"
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
        default="${VAR}_bottom_fluxes.dat",
        type=str,
        help="Name of the output files",
    )


def main(args: argparse.Namespace) -> int:
    if args.cmd != COMMAND_NAME:
        LOGGER.error(
            "Bottom fluxes command has been invoked with command: %s", args.cmd
        )
        return 1

    mask = Mask.from_file(args.mask)

    bottom_fluxes_ds = compute_bottom_fluxes(mask)
    output_dtype = np.float32

    save_dataset(
        bottom_fluxes_ds,
        args.output,
        args.output_name,
        output_dtype,
    )
    LOGGER.info("Execution completed!")

    return 0
