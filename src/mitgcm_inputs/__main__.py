import argparse
import logging
from sys import exit as sys_exit

from mitgcm_inputs.vertical_fluxes import main as vflux_main
from mitgcm_inputs.vertical_fluxes import sub_arguments as vflux_sub_arguments


if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


def configure_logger():
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    LOGGER.setLevel(logging.DEBUG)

    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)

    LOGGER.addHandler(handler)


def argument():
    """
    Generate a subparser for each type of file that can be produced by this
    script and delegate the task of configuring the subparser to the
    corresponding implementation
    """
    parser = argparse.ArgumentParser(description="Generate MITgcm input files")
    subparsers = parser.add_subparsers(
        title="static_type",
        dest="cmd",
        required=True,
    )

    vflux_sub_arguments(subparsers)
    return parser.parse_args()


def main() -> int:
    args = argument()

    cmd_map = {
        "vertical_fluxes": vflux_main,
    }

    def unknown_command(_args: argparse.Namespace):
        LOGGER.error("Invalid command name: %s", _args.cmd)
        return 1

    return cmd_map.get(args.cmd, unknown_command)(args)


if __name__ == "__main__":
    sys_exit(main())
