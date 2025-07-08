import argparse
import logging

import numpy as np
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import path_inside_an_existing_dir

from mitgcm_inputs.k_extinction.k_extinction import compute_k_extinction
from mitgcm_inputs.tools.read_mesh_mask import read_mesh_mask


if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


COMMAND_NAME = "k_extinction"


def sub_arguments(subparser):
    parser = subparser.add_parser(
        COMMAND_NAME, help="Generate K extinction coefficients"
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
        type=path_inside_an_existing_dir,
        help="Path of the output file",
    )


def main(args: argparse.Namespace) -> int:
    if args.cmd != COMMAND_NAME:
        LOGGER.error(
            "K extinction command has been invoked with command: %s", args.cmd
        )
        return 1

    mask = read_mesh_mask(args.mask)

    output_data = compute_k_extinction(mask)
    output_dtype = np.float32

    LOGGER.info("Writing output file: %s", args.output)
    with open(args.output, "wb") as fw:
        fw.write(output_data.kExt.values.astype(output_dtype).tobytes())
    LOGGER.info("Execution completed!")

    return 0
