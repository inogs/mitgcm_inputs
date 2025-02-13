# from legacy code legacy_code/make_surface_deposition_GoT_iNEST.m
# create surface fluxes (atmospheric deposition)
import matplotlib.pyplot as plt
import numpy as np

import mitgcm_inputs.pathes_to_data as ptd
from mitgcm_inputs.parameters import is_plot
from mitgcm_inputs.parameters import nt
from mitgcm_inputs.parameters import nx
from mitgcm_inputs.parameters import ny

print("Start make_surface_deposition_GoT_iNEST")
# Parametri
namef = ["N1p", "N3n"]
# options: print the contouur plot of surface deposition
cont = True

# fields for BFMcoupler_N1pSurfForcFile and BFMcoupler_N3nSurfForcFile
# N1p=1.1055 x 10^(-8) mmol/m^2/s
# N3n=6.5938 x 10^(-7) mmol/m^2/s
# Valori dei flussi superficiali (in mmol/m^2/s)
dep = [1.1055e-8, 6.5938e-7]

# Ciclo per creare i flussi superficiali e salvarli nei file
for var in range(2):
    with open(
        ptd.fluxes_path + f"{namef[var]}_surface_fluxes_GoT_iNEST_V2.dat", "wb"
    ) as fidout:
        for t in range(nt):
            sflux = np.ones((nx, ny)) * dep[var]
            fidout.write(sflux)
        fidout.close()
        if cont:
            # Creazione del grafico
            if is_plot:
                plt.figure()
                plt.pcolor(sflux.T, shading="flat")
                plt.colorbar()
                plt.title(f"Superficie di flusso {namef[var]}")
                plt.show()
print("Done make_surface_deposition_GoT_iNEST")
