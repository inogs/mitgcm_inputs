import numpy as np

is_plot = False

# Grid
nx = 450
ny = 300
nz = 37
nt = 365

# Global Area (Integrated horizontal Area [m^2])
# GA=1.179E+09

# Mean (±standard deviation) annual benthic fluxes of N, P, Si and
# O2 (mmol m−2 d−1).
# Reference              NO−3       NH+4       PO3−4       Si(OH)4    O2
# Bertuzzi et al.(1997)  0.17±0.73  0.8±0.7    0.029±0.05  2.59±2.3   −20.4±8.9
# BFM–POM 1D             0.27±0.16  0.63±0.38  0.048±0.03  0.47±0.41  −5.14±3.5

# Mean (±standard deviation) annual benthic fluxes of N, P, Si and
# O2 (mmol m−2 d−1).
Nb = (0.17 + 0.8) / 86400  # [mmol/m^2/s]
Pb = 0.029 / 86400  # [mmol/m^2/s]
Ob = -20.4 / 86400  # [mmol/m^2/s]

namef = np.array(["N", "P", "B", "C", "O"])
patom = np.array([14, 31, 1, 1, 16])
bentflux = np.array([Nb, Pb, 0.0, 0.0, Ob])
