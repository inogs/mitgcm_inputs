import argparse
import logging


if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


def sub_arguments(subparser):
    parser = subparser.add_parser(
        "vertical_fluxes", help="Generate vertical fluxes"
    )

    # Just an example for the future! This is just a draft
    parser.add_argument("-o", "--output", help="Output directory")


def main(args: argparse.Namespace) -> int:
    if args.cmd != "vertical_fluxes":
        LOGGER.error(
            "Vertical flux command has been invoked with command: %s", args.cmd
        )
        return 1
    return 0
