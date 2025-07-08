import logging
from enum import StrEnum
from numbers import Real
from types import MappingProxyType

import numpy as np
import xarray as xr
from bitsea.commons.mask import Mask


LOGGER = logging.getLogger(__name__)


# This is the value that will be used for the cells that do not contain water
# FILL_VALUE = np.nan
FILL_VALUE = np.float32(0.0)


# Mean (±standard deviation) annual benthic fluxes of N, P, Si and
# O2 (mmol m−2 d−1).
# Reference              NO−3       NH+4       PO3−4       Si(OH)4    O2
# Bertuzzi et al.(1997)  0.17±0.73  0.8±0.7    0.029±0.05  2.59±2.3   −20.4±8.9
# BFM–POM 1D             0.27±0.16  0.63±0.38  0.048±0.03  0.47±0.41  −5.14±3.5


class BenthicVar(StrEnum):
    """
    Enumeration representing benthic variable types.

    This class provides an enumeration for specific benthic variables.
    Each value corresponds to a particular type of benthic variable identified
    by its associated character code.

    Attributes:
        N: Represents benthic variable Nitrogen
        P: Represents benthic variable Phosphorus
        B: Represents benthic variable Biological Oxygen Demand
        C: Represents benthic variable Chemical Oxygen Demand
        O: Represents benthic variable Oxygen
    """

    N = "N"
    P = "P"
    B = "B"
    C = "C"
    O = "O"  # noqa: E741


ATOMIC_WEIGHTS = MappingProxyType(
    {
        BenthicVar.N: 14,
        BenthicVar.P: 31,
        BenthicVar.B: 1,
        BenthicVar.C: 1,
        BenthicVar.O: 16,
    }
)

# Mean (±standard deviation) annual benthic fluxes of N, P, Si and
# O2 (mmol m−2 s−1).
BENTHIC_FLUXES = MappingProxyType(
    {
        BenthicVar.N: (0.17 + 0.8) / 86400,
        BenthicVar.P: 0.029 / 86400,
        BenthicVar.B: 0.0,
        BenthicVar.C: 0.0,
        BenthicVar.O: -20.4 / 86400,
    }
)


FLUX_PARAMETERS = {
    BenthicVar.N: {
        "depth_start_decay": 25.0,
        "depth_zero_flux": 75.0,
        "exp_coeff": 4,
    },
    BenthicVar.P: {
        "depth_start_decay": 25.0,
        "depth_zero_flux": 75.0,
        "exp_coeff": 4,
        "multiplication_factor": 2.0,
    },
    BenthicVar.O: {
        "depth_start_decay": 15.0,
        "depth_zero_flux": 50.0,
        "exp_coeff": 4,
    },
}


def time_variability_factor(t: np.ndarray | Real, _benthic_var: BenthicVar):
    """
    Generate a coefficient that represents the yearly variability of benthic
    fluxes.

    Specifically, the original flux value will be scaled by the output of this
    function, which depends on the day of the year.

    This function is vectorized, so if `t` is a 1D array, the output will also
    be a 1D array of the same shape.

    Args:
        t: The day of the year, where 1 corresponds to January 1st, and 365
            corresponds to December 31st.
        _benthic_var: The benthic variable for which the coefficient is
            calculated. Currently, this parameter is not used, as the same
            coefficient applies to all benthic variables.

    Returns:
        A coefficient that reflects the variability throughout the year.
    """
    return 0.5 * (1.0 - np.cos(2.0 * np.pi / 365.0 * (t + 10)))


def depth_variability_factor(
    depth: np.ndarray | Real,
    benthic_var: BenthicVar,
    depth_start_decay: float = 25.0,
    depth_zero_flux: float = 75.0,
    exp_coeff: int = 4,
    multiplication_factor: float = 1.0,
):
    """
    Compute a coefficient that represents the vertical variability of benthic
    fluxes.

    This function returns the value of the benthic flux at the given depth.

    Args:
        depth: The depth at which the coefficient is calculated.
        benthic_var: The variable for which the coefficient is calculated.
        depth_start_decay: The depth at which the flux starts to decrease. At
            this depth, the flux is exactly equal to the reference values
            stored in the `BENTHIC_FLUXES` map.
        depth_zero_flux: The depth at which the flux is exactly equal to zero.
            For depths that are between `depth_start_decay` and
            `depth_zero_flux`, the value is interpolated. Below this depth, the
            flux is always equal to zero.
        exp_coeff: The exponent used to compute the decay. If it is 1, the
            flux will decrease linearly from `depth_start_decay` to
            `depth_zero_flux`. If it is 2, the flux will decrease quadratically
            and so on.
        multiplication_factor: If this value is different from 1, the output
            of this function will be multiplied by this value.

    Returns:
        The value of the benthic flux at the given depth.
    """
    # bounded_depth = np.clip(
    #     depth, a_max=depth_zero_flux, a_min=depth_start_decay
    # )
    bounded_depth = np.minimum(depth, depth_zero_flux)

    # Value between 0 and 1; it is 1 when we are on depth_start_decay, it is
    # 0 when we are on depth_zero_flux
    depth_fraction = (depth_zero_flux - bounded_depth) / (
        depth_zero_flux - depth_start_decay
    )

    benthic_flux = BENTHIC_FLUXES[benthic_var]
    main_coeff = depth_fraction**exp_coeff

    return benthic_flux * main_coeff * multiplication_factor


def compute_bottom_fluxes(meshmask: Mask):
    n_days = 365

    depth = meshmask.bathymetry()
    water_cells = meshmask[0]

    dataset = xr.Dataset(
        coords={
            "day": np.arange(1, n_days + 1),
            "latitude": (("x", "y"), meshmask.ylevels),
            "longitude": (("x", "y"), meshmask.xlevels),
        },
    )

    for variable in BenthicVar:
        LOGGER.debug("Computing bottom fluxes for %s", variable)
        data = np.full(depth.shape, FILL_VALUE, dtype=np.float32)
        data[water_cells] = depth_variability_factor(
            depth[water_cells], variable, **FLUX_PARAMETERS.get(variable, {})
        )
        time_coefficients = time_variability_factor(
            np.arange(1, 366), variable
        )
        dataset[variable] = xr.DataArray(
            data=time_coefficients[:, np.newaxis, np.newaxis] * data,
            dims=("day", "x", "y"),
        )
    return dataset
