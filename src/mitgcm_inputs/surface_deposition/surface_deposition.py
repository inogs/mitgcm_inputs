import logging
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np
import xarray as xr
from bitsea.commons.mask import Mask


LOGGER = logging.getLogger(__name__)


VAR_VALUES = MappingProxyType(
    {
        "N1p": 1.1055e-8,
        "N3n": 6.5938e-7,
    }
)


def compute_surface_deposition(
    meshmask: Mask, var_values: Mapping[str, float] = VAR_VALUES
):
    LOGGER.info("Computing surface deposition")

    LOGGER.debug("Meshmask shape: %s", meshmask.shape)
    surface_decomposition_dataset = xr.Dataset(
        coords={
            "latitude": (("x", "y"), meshmask.ylevels),
            "longitude": (("x", "y"), meshmask.xlevels),
        },
    )

    for var_name, var_value in var_values.items():
        LOGGER.debug("Processing variable: %s", var_name)
        var_data = xr.DataArray(
            data=np.broadcast_to(var_value, meshmask.shape[1:]),
            dims=("x", "y"),
        )
        surface_decomposition_dataset[var_name] = var_data

    return surface_decomposition_dataset
