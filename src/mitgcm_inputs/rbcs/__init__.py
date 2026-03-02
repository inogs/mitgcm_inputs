import argparse
import io
import json
import logging
import tarfile
from collections import OrderedDict
from pathlib import Path

from bitsea.utilities.argparse_types import dir_to_be_created_if_not_exists
from bitsea.utilities.argparse_types import existing_file_path
from ogs_riverger.read_config import RiverConfig

from mitgcm_inputs.rbcs.rbcs_gen import build_conc_and_relax_variables
from mitgcm_inputs.rbcs.rbcs_gen import build_domain_config_string
from mitgcm_inputs.rbcs.rbcs_gen import get_spatial_description_from_meshmask
from mitgcm_inputs.rbcs.scarichi_json_gen import read_sewage_positions
from mitgcm_inputs.tools.tar_utils import set_tar_file_ownerships


if __name__ == "__main__":
    LOGGER = logging.getLogger()
else:
    LOGGER = logging.getLogger(__name__)


COMMAND_NAME = "rbcs"

# The current directory of this file
MAIN_DIR = Path(__file__).resolve().parent
DEFAULT_SEWAGE_FILE = (
    MAIN_DIR / "Allegato_1_Capitolato_Tecnico_B32_B35_scarichi.xlsx"
)


def sub_arguments(subparser):
    parser = subparser.add_parser(
        COMMAND_NAME,
        help="Generate the files conc.tar.gz and relax.tar.gz that store the "
        "relaxation and concentration values for the tracers",
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
        "-s",
        "--sewage",
        required=False,
        type=existing_file_path,
        default=DEFAULT_SEWAGE_FILE,
        help="The xlsx file with the sewage locations",
    )

    parser.add_argument(
        "-p",
        "--river-positions",
        required=False,
        type=existing_file_path,
        default=None,
        help="The json file that describe the position of the mouth of the"
        "rivers",
    )

    parser.add_argument(
        "-r",
        "--river-config",
        required=False,
        type=existing_file_path,
        default=None,
        help="The main configuration file for the rivers; it can be ignored "
        "if there are no rivers in the domain but it is mandatory if a "
        "river positions file is provided",
    )

    parser.add_argument(
        "-d",
        "--river-domain-file",
        required=False,
        type=existing_file_path,
        default=None,
        help="A further configuration for the rivers, specific to the current "
        "domain; if it is `None`, it will be ignored ",
    )


def main(args: argparse.Namespace) -> int:
    if args.cmd != COMMAND_NAME:
        LOGGER.error(
            "RBCS command has been invoked with command: %s",
            args.cmd,
        )
        return 1

    # Check if the arguments from the command line are coherent
    if args.river_config is None:
        if args.river_positions is not None:
            LOGGER.error(
                "If river positions are provided, a river config must be "
                "provided"
            )
            return 101
        if args.river_domain_file is not None:
            LOGGER.error(
                "If river domain file is provided, a river config must be "
                "provided"
            )
            return 102
    if args.river_config is not None:
        if args.river_positions is None:
            LOGGER.error(
                "If river config is provided, river positions must be provided"
            )
            return 103

    # Read the meshmask to build the spatial description of the domain
    # and then read the position of the sewage points
    spatial_description = get_spatial_description_from_meshmask(args.mask)
    sewage_positions = read_sewage_positions(args.sewage, args.mask)

    # If the river config file is available, read it using ogs_riverger
    if args.river_config is not None:
        rivers_config = RiverConfig.from_json(
            args.river_config, args.river_domain_file
        )
    else:
        rivers_config = RiverConfig(root=OrderedDict())

    # This function checks if a river has the E. coli tracer enabled
    def has_tracer(river):
        river_id = river["id"]
        if river_id not in rivers_config.root:
            raise IndexError(
                f"River id {river_id} not found in the main river config"
            )
        return len(rivers_config.root[river_id].concentrations) > 0

    # Read the rivers_positions.json file
    if args.river_positions is not None:
        river_positions = json.loads(args.river_positions.read_text())
    else:
        river_positions = tuple()

    # Keep only rivers with tracers
    rivers = tuple(r for r in river_positions if has_tracer(r))

    sources, conc_and_relax = build_conc_and_relax_variables(
        spatial_description=spatial_description,
        sewage_points=sewage_positions.to_dict(orient="records"),
        rivers=rivers,
    )

    if len(sources) == 0:
        LOGGER.info(
            "No tracers found for this domain, exiting without writing any "
            "file"
        )
        return 0

    conc_file_path = args.output / "conc.tar.gz"
    LOGGER.info("Compressing conc files into %s", conc_file_path)
    with tarfile.open(conc_file_path, "w:gz") as tar:
        # Create the directory conc inside the file
        relax_dir = tarfile.TarInfo(name="conc")
        relax_dir.type = tarfile.DIRTYPE
        set_tar_file_ownerships(relax_dir)
        relax_dir.mode = 0o755
        tar.addfile(relax_dir)

        for data_var in conc_and_relax.data_vars:
            if not data_var.startswith("conc"):
                continue
            current_conc_data = (
                conc_and_relax[data_var].values.astype("f4", copy=False).tobytes()
            )
            current_conc = tarfile.TarInfo(f"conc/{data_var}.bin")
            set_tar_file_ownerships(current_conc)
            current_conc.mode = 0o644
            current_conc.size = len(current_conc_data)
            tar.addfile(current_conc, io.BytesIO(current_conc_data))

    config_string = build_domain_config_string(sources)
    point_source_output = args.output / "RBCS_names.txt"
    LOGGER.info("Writing file %s", point_source_output)
    point_source_output.write_text(config_string)

    if not any(s["kind"] == "sewage" for s in sources):
        LOGGER.info("No sewage inside this domain, skipping relaxation files")
        LOGGER.info("Execution completed!")
        return 0

    relax_files = {
        "relaxed_salinity": "bottom_sources_S_relaxation.bin",
        "salinity_mask": "bottom_sources_S_mask.bin",
    }
    relax_file_path = args.output / "relax.tar.gz"
    LOGGER.info("Compressing relax files into %s", relax_file_path)
    with tarfile.open(relax_file_path, "w:gz") as tar:
        # Create the directory relax inside the file
        relax_dir = tarfile.TarInfo(name="relax")
        relax_dir.type = tarfile.DIRTYPE
        set_tar_file_ownerships(relax_dir)
        relax_dir.mode = 0o755
        tar.addfile(relax_dir)

        for var_name, file_name in relax_files.items():
            current_data_array = conc_and_relax[var_name].values.astype(
                "f4", copy=False
            )
            current_data = current_data_array.tobytes()
            tar_file_pointer = tarfile.TarInfo(f"relax/{file_name}")
            set_tar_file_ownerships(tar_file_pointer)
            tar_file_pointer.mode = 0o644
            tar_file_pointer.size = len(current_data)
            tar.addfile(tar_file_pointer, io.BytesIO(current_data))

    LOGGER.info("Execution completed!")

    return 0
