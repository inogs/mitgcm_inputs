import logging
from pathlib import Path

import numpy as np
import xarray as xr

LOGGER = logging.getLogger(__name__)


def save_dataset(
    dataset: xr.Dataset,
    output_dir: Path,
    output_name_mask: str,
    output_dtype: np.typing.DTypeLike = np.float32,
):
    for var_name in dataset.data_vars:
        file_output_name = output_name_mask.replace("${VAR}", var_name)
        file_output_path = output_dir / file_output_name

        LOGGER.info("Writing output file: %s", file_output_path)
        with open(file_output_path, "wb") as fw:
            fw.write(dataset[var_name].values.astype(output_dtype).tobytes())
