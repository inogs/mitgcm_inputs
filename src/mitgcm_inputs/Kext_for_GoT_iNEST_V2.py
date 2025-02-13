# from legacy code Kext_for_GoT_iNEST_V2.m
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

import mitgcm_inputs.pathes_to_data as ptd
from mitgcm_inputs.parameters import is_plot


def plot2D(is_plot, A):
    if is_plot:
        plt.figure()
        plt.pcolor(A)
        plt.colorbar()
        plt.show()
    else:
        print("Not plotting figure")
    return 0


print("Start Kext_for_GoT_iNEST_V2")

if not is_plot:
    print("No figures printed")
# Path to files
data_file = ptd.med_data_path
strdir_grid = ptd.med_grid_path
mask_file_got = ptd.bathymetry_path + "/hFacC.data"

# parameters mediterran sea
nx = 570
ny = 264
nz = 75
days = 1921

# parameters Gulf of Trieste
yc_origin = 45.4
xc_origin = 13.22
delX = 1 / 768
delY = 1 / 768
nx_got = 450
ny_got = 300
nz_got = 40
# load mask for GoT
fr2 = open(mask_file_got, "rb")
mask_got = np.fromfile(fr2, dtype=np.float32).reshape((nz_got, ny_got, nx_got))
fr2.close()
mask_got = mask_got[0, :, :]
print("SHAPE GoT", np.shape(mask_got))
plot2D(is_plot, mask_got)

fhfacC = open(f"{strdir_grid}Computational_hFacC", "rb")
mask_hfac = np.fromfile(fhfacC, dtype=">f4").reshape((nz, ny, nx))
fhfacC.close()

mask_hfac[mask_hfac != 1] = 0
mask_hfac = mask_hfac[0, :, :]
print("SHAPE Med", np.shape(mask_hfac))

plot2D(is_plot, mask_hfac)
kext = np.zeros((days, ny, nx))
fKext = open(data_file, "rb")
A = np.zeros((days, ny, nx))
A = np.fromfile(fKext, dtype=">f4", count=days * nx * ny).reshape(
    (days, ny, nx)
)
fKext.close()
for idays in range(days):
    kext[idays, :, :] = A[0, :, :] * mask_hfac
plot2D(is_plot, kext[0, :, :])

mask_nan = mask_hfac.copy()
mask_nan[mask_nan == 0] = np.nan
plot2D(is_plot, mask_nan)
annual = np.nanmean(kext[92 : 92 + 364, :, :], axis=0)
annual = np.squeeze(annual) * mask_nan
print("Shape annual", np.shape(annual))

plot2D(is_plot, annual)

# ora interpolo sulla griglia GoT
mlat = np.arange(yc_origin, yc_origin + ny_got * delY, delY)
mlon = np.arange(xc_origin, xc_origin + nx_got * delX, delX)
clat, clon = np.meshgrid(mlat, mlon)
print("Shape lat and lon", np.shape(mlat), np.shape(mlon))

fh_lon = open(f"{strdir_grid}Computational_XC", "rb")
lon = np.fromfile(fh_lon, dtype=">f4").reshape((ny, nx))
fh_lon.close()
fh_lat = open(f"{strdir_grid}Computational_YC", "rb")
lat = np.fromfile(fh_lat, dtype=">f4").reshape((ny, nx))
fh_lat.close()

lat_nan = lat * mask_nan
lon_nan = lon * mask_nan

varint = griddata(
    (lon[~np.isnan(mask_nan)], lat[~np.isnan(mask_nan)]),
    annual[~np.isnan(mask_nan)],
    (clon, clat),
    method="linear",
)
print("Shape varint", np.shape(varint))

# plot2D(is_plot,clat)
# plot2D(is_plot,clon)
# plot2D(is_plot,lat_nan)
# plot2D(is_plot,lon_nan)
plot2D(is_plot, varint.T)

if is_plot:
    plt.figure()
    plt.pcolor(annual)
    plt.colorbar()
    # Set color bar limits
    plt.clim(0.06, 0.22)
    plt.show()

    plt.figure()
    plt.pcolor((varint.T * mask_got))
    plt.colorbar()
    plt.clim(0.06, 0.22)
    plt.show()

fw = open(ptd.kext_path + "Kext_GoT_iNEST_365.dat", "wb")
for _ in range(365):
    fw.write(varint.astype(np.float32).tobytes())
fw.close()
print("Done Kext_for_GoT_iNEST_V2")
