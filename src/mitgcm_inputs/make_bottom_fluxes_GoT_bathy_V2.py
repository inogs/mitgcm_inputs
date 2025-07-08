# from legacy code make_bottom_fluxes_GoT_bathy_V2.m
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


# Create bottom fluxes (no point sources)
print("Start make_bottom_fluxes_GoT_bathy_V2")
curdir = os.getcwd()

# Global Area (Integrated horizontal Area [m^2])
# GA=1.179E+09

levels = [6, 29, 2]
check_z_levels(levels, p.nz)
delZ = np.array(
    [0.5] * levels[0] + [1.0] * levels[1] + [2.0] * levels[2]
).reshape(-1, 1)
# print("DetaZ",delZ)

# Load masks
mask_bottom = np.fromfile(
    ptd.mask_path + "mask_bottom_GoT_iNEST_V2", dtype=int
).reshape(p.ny, p.nx)
mask_qbottom = np.fromfile(
    ptd.mask_path + "mask_qbottom_GoT_iNEST_V2", dtype=int
).reshape(p.ny, p.nx)
mask_surface = np.fromfile(
    ptd.mask_path + "mask_surface_GoT_iNEST_V2", dtype="float32"
).reshape(p.ny, p.nx)


# Load bathymetry data
os.chdir(ptd.bathymetry_path)

dxg = np.fromfile("DXG.data", dtype="float32").reshape(p.ny, p.nx)
dyg = np.fromfile("DYG.data", dtype="float32").reshape(p.ny, p.nx)
btm = np.fromfile("Depth.data", dtype="float32").reshape(p.ny, p.nx)

print("dxg,dyg, btm shapes")
print(np.shape(dxg))
print(np.shape(dyg))
print(np.shape(btm))
print("DXG limits ", np.min(dxg), np.max(dxg))
print("DYG limits ", np.min(dyg), np.max(dyg))
print("Depth limits ", np.min(btm), np.max(btm))

# FLUXES
# minlev=6; % 6(AZAL)=~15 m; 3.0010 (profondità minima)
# minlev=8; % 8(AZAL)=~25 m;
# minlev=28; % 28(AZAL_HR)=~25 m; 3.00 (profondità minima)
minlev = 18  # 18(GoT_iNEST)=~15 m; 3.00 (profondità minima)
mindepth = np.sum(delZ[:minlev])
# levlim=15; % -> ~70 m; levlim=12(AZAL) -> W-Interf. bottom 48.8900 (limite dei ~50 m)
# levlim=45; % -> ~70 m; levlim=42(AZAL_HR) -> limite dei ~50 m
levlim = 37  # -> ~35 m; levlim=37(GoT_iNEST) -> limite dei ~35 m
exp = 4
mult = 2
lim = np.sum(delZ[0:levlim])
print("min depth", mindepth, "lev min", lim)
btml = btm.copy()
btml[btml >= lim] = lim  # Cut the bottom at ~35 m
sp = mult * ((mask_surface * lim - btml) / (lim - mindepth)) ** exp

btml = mask_bottom.copy()
btml[btml >= levlim] = levlim  # Cut the bottom at levlim
sp2 = mult * ((mask_surface * levlim - btml) / (levlim - minlev)) ** exp
# plot
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
    fidout = open(f"{p.namef[var]}_bottom_fluxes_GoT_iNEST_V2_x2.dat", "wb")
    tflux = np.zeros(p.nt)
    for t in range(p.nt):
        time = t + 1
        cy = 0.5 * (
            1 - np.cos(((2 * np.pi / 365) * time) + 2 * np.pi / 365 * (10 + 0))
        )  # ok<SAGROW>
        # c1[t]=cy;
        # if t>=210:
        # y(t)=cy(time).*1/(time-209).^0.5;
        # c1(t)=c1(time);
        bflux = sp * cy * p.bentflux[var]
        fidout.write(bflux.astype("float32").tobytes())

        bf = bflux * dxg * dyg
        bf[bf == 0] = np.nan
        if np.any(~np.isnan(bf)):
            tflux[t] = np.nansum(bf) * 86400 / 1e9 * p.patom[var]  # [tonn/gg]
        # else:
        # print(tflux[t])

    if p.nt > 0:
        yrflux = np.sum(tflux)
        print(f"{p.namef[var]} integral flux", yrflux)
        if p.is_plot:
            plt.figure()
            plt.pcolor(bf.T)
            plt.title(f"{p.namef[var]}")
            plt.colorbar()

            plt.show()
    fidout.close()
print("Done make_bottom_fluxes_GoT_bathy_V2")
