import argparse
import logging

import xarray as xr
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import path_inside_an_existing_dir

if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


COMMAND_NAME = "exf_albedo"


def sub_arguments(subparser):
    parser = subparser.add_parser(
        COMMAND_NAME, help="Compute the EXF albedo for a specific domain"
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


def compute_albedo(average_lat: float) -> float:
    return 0.06 + 0.0022 * (average_lat - 37)


def main(args: argparse.Namespace) -> int:
    if args.cmd != COMMAND_NAME:
        LOGGER.error(
            "EXF albedo command has been invoked with command: %s", args.cmd
        )
        return 1

    LOGGER.debug("Opening file %s", args.mask)
    with xr.open_dataset(args.mask) as ds:
        min_lat = float(ds.attrs["domain_minimum_latitude"])
        max_lat = float(ds.attrs["domain_maximum_latitude"])
    LOGGER.debug("Domain latitude range: %s - %s", min_lat, max_lat)

    average_lat = (min_lat + max_lat) / 2
    albedo = compute_albedo(average_lat)

    LOGGER.debug("Writing output to file %s", args.output)
    with open(args.output, "w") as f:
        f.write(f"exf_albedo={albedo:.3f}\n")

    LOGGER.info("Execution completed!")

    return 0
