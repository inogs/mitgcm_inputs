import argparse
import logging
import tarfile
from os import fspath
from pathlib import Path
from sys import exit as sys_exit
from tempfile import TemporaryDirectory
from time import time

from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import path_inside_an_existing_dir

from mitgcm_inputs.bottom_fluxes import COMMAND_NAME as BFLUX_COMMAND_NAME
from mitgcm_inputs.bottom_fluxes import main as bflux_main
from mitgcm_inputs.bottom_fluxes import sub_arguments as bflux_sub_arguments
from mitgcm_inputs.k_extinction import COMMAND_NAME as KEXT_COMMAND_NAME
from mitgcm_inputs.k_extinction import main as kext_main
from mitgcm_inputs.k_extinction import sub_arguments as kext_sub_arguments
from mitgcm_inputs.surface_deposition import COMMAND_NAME as SD_COMMAND_NAME
from mitgcm_inputs.surface_deposition import main as sd_main
from mitgcm_inputs.surface_deposition import sub_arguments as sd_sub_arguments


if __name__ == "__main__" or __name__ == "mitgcm_inputs.__main__":
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


def argument(sys_argv=None):
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
    sd_sub_arguments(subparsers)
    bflux_sub_arguments(subparsers)

    subparser = subparsers.add_parser(
        "ALL", help="Execute all the available commands"
    )
    subparser.add_argument(
        "-m",
        "--mask",
        required=True,
        type=existing_file_path,
        help="The meshmask file that defines the domain",
    )

    subparser.add_argument(
        "-o",
        "--output",
        required=True,
        type=path_inside_an_existing_dir,
        help="""
        Path of the compressed tar.gz file that will store all the output
        files
        """,
    )

    if sys_argv is not None:
        return parser.parse_args(sys_argv)
    else:
        return parser.parse_args()


def execute_all_commands(args: argparse.Namespace):
    mask_path = args.mask
    output_file = args.output

    with TemporaryDirectory() as tmp_dir:
        tmp_dir_path = Path(tmp_dir)

        command_args = {
            BFLUX_COMMAND_NAME: [
                BFLUX_COMMAND_NAME,
                "-m",
                fspath(mask_path),
                "-o",
                fspath(tmp_dir_path),
            ],
            SD_COMMAND_NAME: [
                SD_COMMAND_NAME,
                "-m",
                fspath(mask_path),
                "-o",
                fspath(tmp_dir_path),
            ],
            KEXT_COMMAND_NAME: [
                KEXT_COMMAND_NAME,
                "-m",
                fspath(mask_path),
                "-o",
                fspath(tmp_dir_path / "Kext.bin"),
            ],
        }

        for command_name, command in CMD_MAP.items():
            if command_name.lower() == "all":
                continue
            current_args = command_args[command_name]
            LOGGER.info("Executing %s", command_name)
            LOGGER.debug(
                "Executing %s with the following args: %s",
                command_name,
                current_args,
            )
            command(argument(current_args))

        output_stem = output_file.stem
        if output_stem.lower().endswith(".tar"):
            output_stem = Path(output_stem).stem

        LOGGER.info("Compressing data into file %s", output_file)
        with tarfile.open(output_file, "w:gz") as tar:
            for file_path in tmp_dir_path.iterdir():
                tar_relative_path = Path(output_stem) / file_path.name
                tar.add(
                    file_path,
                    arcname=tar_relative_path,
                )

            # Get one random file inside the tar
            reference_file = next(iter(tar.getmembers()))

            # Add the directory where the files are stored as a tar object
            # (otherwise MER produces an error)
            data_dir = tarfile.TarInfo(name=output_stem)
            data_dir.type = tarfile.DIRTYPE
            data_dir.mtime = int(time())
            data_dir.uid = reference_file.uid
            data_dir.gid = reference_file.gid
            data_dir.uname = reference_file.uname
            data_dir.gname = reference_file.gname

            tar.addfile(data_dir)

        LOGGER.info("Everything done!")
    return 0


CMD_MAP = {
    BFLUX_COMMAND_NAME: bflux_main,
    SD_COMMAND_NAME: sd_main,
    KEXT_COMMAND_NAME: kext_main,
    "ALL": execute_all_commands,
}


def main() -> int:
    args = argument()
    configure_logger()

    def unknown_command(_args: argparse.Namespace):
        LOGGER.error("Invalid command name: %s", _args.cmd)
        return 1

    return CMD_MAP.get(args.cmd, unknown_command)(args)


if __name__ == "__main__":
    sys_exit(main())
