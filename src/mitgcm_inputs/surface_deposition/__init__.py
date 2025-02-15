import argparse
import logging

import numpy as np
from bitsea.commons.mask import Mask
from bitsea.utilities.argparse_types import dir_to_be_created_if_not_exists
from bitsea.utilities.argparse_types import existing_file_path

from mitgcm_inputs.surface_deposition.surface_deposition import (
    compute_surface_deposition,
)


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
        help="Path of the directyory where the output files will be written",
    )

    parser.add_argument(
        "-z",
        "--output-name",
        required=False,
        default="${VAR}_surface_fluxes.dat",
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

    mask = Mask.from_file(args.mask)

    output_data = compute_surface_deposition(mask)
    output_dtype = np.float32

    for var_name in output_data.data_vars:
        file_output_name = args.output_name.replace("${VAR}", var_name)
        file_output_path = args.output / file_output_name

        LOGGER.info("Writing output file: %s", file_output_path)
        with open(file_output_path, "wb") as fw:
            fw.write(
                output_data[var_name].values.astype(output_dtype).tobytes()
            )

    LOGGER.info("Execution completed!")

    return 0
