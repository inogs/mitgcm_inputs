import logging
from enum import StrEnum
from numbers import Real

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
    by its associated character code. It is a subclass of StrEnum to ensure the
    values are treated as strings.

    Attributes:
        N: Represents benthic variable Nitrogen
        P: Represents benthic variable Phosphorus
        B: Represents benthic variable Boron
        C: Represents benthic variable Carbon
        O: Represents benthic variable Oxygen
    """

    N = "N"
    P = "P"
    B = "B"
    C = "C"
    O = "O"  # noqa: E741


p_atom = {
    BenthicVar.N: 14,
    BenthicVar.P: 31,
    BenthicVar.B: 1,
    BenthicVar.C: 1,
    BenthicVar.O: 16,
}

# Mean (±standard deviation) annual benthic fluxes of N, P, Si and
# O2 (mmol m−2 s−1).
benthic_fluxes = {
    BenthicVar.N: (0.17 + 0.8) / 86400,
    BenthicVar.P: 0.029 / 86400,
    BenthicVar.B: 0.0,
    BenthicVar.C: 0.0,
    BenthicVar.O: -20.4 / 86400,
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
    depth: np.ndarray | Real, benthic_var: BenthicVar
):
    """
    Compute a coefficient that represents the vertical variability of benthic
    fluxes.

    This function returns the value of the benthic flux at the given depth.

    Args:
        depth: The depth at which the coefficient is calculated.
        benthic_var: The variable for which the coefficient is calculated.

    Returns:
        The value of the benthic flux at the given depth.
    """
    # Parameters
    exp_coeff = 3
    ref_depth = 35.0  # meters
    min_depth = 3.0  # meters

    benthic_flux = benthic_fluxes[benthic_var]
    main_coeff = np.maximum(
        ((ref_depth - depth) / (ref_depth - min_depth)) ** exp_coeff, 0
    )
    return benthic_flux * main_coeff


def compute_bottom_fluxes(meshmask: Mask):
    n_days = 365

    depth = meshmask.bathymetry()

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
        data[meshmask[0]] = depth_variability_factor(
            depth[meshmask[0]], variable
        )
        time_coefficients = time_variability_factor(
            np.arange(1, 366), variable
        )
        dataset[variable] = xr.DataArray(
            data=time_coefficients[:, np.newaxis, np.newaxis] * data,
            dims=("day", "x", "y"),
        )
    return dataset
