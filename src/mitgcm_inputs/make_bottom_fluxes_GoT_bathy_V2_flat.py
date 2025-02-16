# from legacy code make_bottom_fluxes_GoT_bathy_V2_flat.m
import os

import matplotlib.pyplot as plt
import numpy as np

import mitgcm_inputs.parameters as p
import mitgcm_inputs.pathes_to_data as ptd


def check_z_levels(lev, nz):
    if np.sum(lev) != nz:
        print("Wrong vertical discretization")
        check = 0
    else:
        check = 1
    return check


# create bottom fluxes (no point sources)
print("Start make_bottom_fluxes_GoT_bathy_V2_flat")
curdir = os.getcwd()

# Global Area (Integrated horizontal Area [m^2])
# GA=1.179E+09
levels = [6, 29, 2]
check_z_levels(levels, p.nz)
delZ = np.array(
    [0.5] * levels[0] + [1.0] * levels[1] + [2.0] * levels[2]
).reshape(-1, 1)
# print("DetaZ",delZ)


mask_bottom = np.fromfile(
    ptd.mask_path + "mask_bottom_GoT_iNEST_V2", dtype=int
).reshape(p.ny, p.nx)
mask_surface = np.fromfile(
    ptd.mask_path + "mask_surface_GoT_iNEST_V2", dtype="float32"
).reshape(p.ny, p.nx)

# Load bathymetry data
os.chdir(ptd.bathymetry_path)

dxg = np.fromfile("DXG.data", dtype="float32").reshape(p.ny, p.nx)
dyg = np.fromfile("DYG.data", dtype="float32").reshape(p.ny, p.nx)
btm = np.fromfile("Depth.data", dtype="float32").reshape(p.ny, p.nx)

# FLUXES
# minlev=6; % 6(AZAL)=~15 m; 3.0010 (profondità minima)
# minlev=19; % 19(AZAL_HR)=~15 m; 3.00 (profondità minima)
minlev = 15  # 15(GoT_iNEST)=~11 m; 3.00 (profondità minima)
mindepth = np.sum(delZ[:minlev])
# levlim=15; % -> ~70 m; levlim=12(AZAL) -> W-Interf. bottom 48.8900 (limite dei ~50 m)
# levlim=45; % -> ~70 m; levlim=42(AZAL_HR) -> limite dei ~50 m
levlim = 37  # -> ~35 m; levlim=37(GoT_iNEST) -> limite dei ~35 m
exp = 3
lim = np.sum(delZ[:levlim])
btml = btm.copy()
btml[btml >= lim] = lim  # taglio il fondo a ~35 m
sp = ((mask_surface * lim - btml) / (lim - mindepth)) ** exp


btml = mask_bottom.copy()
btml[btml >= levlim] = levlim  # taglio il fondo al livello levlim
sp2 = ((mask_surface * levlim - btml) / (levlim - minlev)) ** exp
if p.is_plot:
    plt.figure()
    plt.pcolor(sp.T)
    plt.colorbar()
    plt.title(f"minlev={minlev} levlim={levlim} exp={exp}")
    plt.show()

    plt.figure()
    plt.pcolor(sp2.T)
    plt.colorbar()
    plt.title(f"minlev={minlev} levlim={levlim} exp={exp}")
    plt.show()

# bottom fluxes
os.chdir(ptd.bottom_fluxes_path)
for var in range(len(p.namef)):
    fidout = open(f"{p.namef[var]}_bottom_fluxes_GoT_iNEST_V2_flat.dat", "wb")
    tflux = np.zeros(p.nt)
    for t in range(p.nt):
        time = t + 1
        cy = 0.5 * (
            1 - np.cos(((2 * np.pi / 365) * time) + 2 * np.pi / 365 * 10)
        )  # ok<SAGROW
        bflux = sp * cy * p.bentflux[var]
        fidout.write(bflux.astype("float32").tobytes())

        bf = bflux * dxg * dyg
        bf[bf == 0] = np.nan
        if np.any(~np.isnan(bf)):
            tflux[t] = np.nansum(bf) * 86400 / 1e9 * p.patom[var]  # [tonn/gg]
        # else:
        # print(tflux[t])
    yrflux = np.sum(tflux)
    fidout.write(bflux.astype("float32").tobytes())
    fidout.close()
    print(f"{p.namef[var]} integral flux", yrflux)
    if p.is_plot:
        plt.figure()
        plt.pcolor(bf.T)
        # plt.shading('flat')
        plt.colorbar()
        plt.title(f"{p.namef[var]}")
        plt.show()
print("Done make_bottom_fluxes_GoT_bathy_V2_flat")
