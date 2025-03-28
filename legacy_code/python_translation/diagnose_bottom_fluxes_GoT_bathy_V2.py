# from legacy code diagnose_bottom_fluxes_GoT_bathy_V2.py
import os

import matplotlib.pyplot as plt
import numpy as np

import legacy_code.python_translation.parameters as p
import legacy_code.python_translation.pathes_to_data as ptd

# diagnose bottom fluxes (no point sources)
print("Start diagnose_bottom_fluxes_GoT_bathy_V2")
if not p.is_plot:
    print("No figures printed")
curdir = os.getcwd()

# Load surface mask
surf = np.fromfile(
    ptd.mask_path + "mask_surface_GoT_iNEST_V2", dtype="float32"
).reshape(p.ny, p.nx)
surf[surf == 0] = np.nan

# Load bathymetry data
os.chdir(ptd.bathymetry_path)

dxg = np.fromfile("DXG.data", dtype="float32").reshape(p.ny, p.nx)
dyg = np.fromfile("DYG.data", dtype="float32").reshape(p.ny, p.nx)

# FLUXES
os.chdir(ptd.bottom_fluxes_path)
for var in range(len(p.namef)):
    if p.namef[var] == "N":
        fidin = open(f"{p.namef[var]}_bottom_fluxes_GoT_iNEST_V2_x2.dat", "rb")
    elif p.namef[var] == "P":
        fidin = open(
            f"{p.namef[var]}_bottom_fluxes_GoT_iNEST_V2_x2.dat", "rb"
        )  ## io ho solo x2
    else:
        fidin = open(
            f"{p.namef[var]}_bottom_fluxes_GoT_iNEST_V2_flat.dat", "rb"
        )

    tbflux = np.zeros((p.ny, p.nx, p.nt), dtype="float32")
    tbf = np.zeros((p.ny, p.nx, p.nt), dtype="float32")
    for t in range(p.nt):
        bflux = np.fromfile(fidin, dtype="float32", count=p.nx * p.ny).reshape(
            (p.ny, p.nx)
        )  # [mmol/m^2/s]
        tbflux[:, :, t] = bflux  # [mmol/m^2/s]
        tbf[:, :, t] = bflux * dxg * dyg  # [mmol/s]
    tbfint = np.sum(tbf, axis=(0, 1))

    tonnyr = np.sum(tbfint) * 86400 / 1e9 * p.patom[var]  # [tonn/yr]
    print(f"{p.namef[var]}: {tonnyr:.2e} [tonn/yr]")

    # plt.figure()
    if p.is_plot:
        plt.plot(tbfint)  # [mmol/s]
        plt.title(
            f"instant flux of {p.namef[var]} - spatially integrated [mmol/s]"
        )
        plt.savefig(f"instant_flux_{p.namef[var]}.png", dpi=600)

    mbflux = np.mean(tbflux, axis=2)
    if p.is_plot:
        plt.figure()
        plt.pcolor(mbflux * surf)
        plt.colorbar()
        plt.title(f"mean flux of {p.namef[var]} [mmol/m^2/s]")
        plt.show()
        plt.savefig(f"mean_flux_{p.namef[var]}.png", dpi=600)

    # plt.figure()
    if p.is_plot:
        for stat in range(20, 241, 20):
            plt.plot(tbflux[200, stat, :], "b")
            plt.plot(
                np.repeat(np.mean(tbflux[165, stat, :], axis=0), 365), "--r"
            )
        plt.title(
            f"mean bottom flux of {p.namef[var]} = {np.mean(mbflux):.2f} [mmol/m^2/s]"
        )
        plt.show()
        plt.savefig(f"time_series_flux_transect_{p.namef[var]}.png", dpi=600)

    fidin.close()
print("Done diagnose_bottom_fluxes_GoT_bathy_V2")
