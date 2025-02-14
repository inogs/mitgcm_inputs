import argparse
import logging
from sys import exit as sys_exit

from mitgcm_inputs.k_extinction import COMMAND_NAME as kext_command_name
from mitgcm_inputs.k_extinction import main as kext_main
from mitgcm_inputs.k_extinction import sub_arguments as kext_sub_arguments
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

    kext_sub_arguments(subparsers)
    vflux_sub_arguments(subparsers)
    return parser.parse_args()


def main() -> int:
    args = argument()
    configure_logger()

    cmd_map = {
        "vertical_fluxes": vflux_main,
        kext_command_name: kext_main,
    }

    def unknown_command(_args: argparse.Namespace):
        LOGGER.error("Invalid command name: %s", _args.cmd)
        return 1

    return cmd_map.get(args.cmd, unknown_command)(args)


if __name__ == "__main__":
    sys_exit(main())
