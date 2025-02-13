# from legacy code make_masks_GoT.m
import os

import numpy as np

import mitgcm_inputs.parameters as p
import mitgcm_inputs.pathes_to_data as ptd

# create masks
# path change
cwd = os.getcwd()
os.chdir(ptd.bathymetry_path)

hFacC_3D = np.zeros((p.nx, p.ny, p.nz), dtype=np.float32)
with open("hFacC.data", "rb") as fidin:
    for iz in range(p.nz):
        C = np.fromfile(fidin, dtype=np.float32, count=p.nx * p.ny)
        C = C.reshape(p.nx, p.ny)
        hFacC_3D[:, :, iz] = C
fidin.close()
mask = hFacC_3D
surf = mask[:, :, 0]  # attention:ndarray is not C-contiguous

os.chdir(cwd)

qbottom = np.zeros(
    (p.nx, p.ny), dtype=int
)  # penultimate water cell (above the bottom water cell)
bottom = np.zeros((p.nx, p.ny), dtype=int)  # bottom water cell
datum = np.zeros(
    (p.nx, p.ny), dtype=np.float32
)  # hFacC of the bottom water cell

for iz in range(p.nz - 1):
    for i in range(p.nx):
        for j in range(p.ny):
            if mask[i, j, iz] == 1.0 and mask[i, j, iz + 1] != 0.0:
                qbottom[i, j] = iz

for i in range(p.nx):
    for j in range(p.ny):
        if qbottom[i, j] != 0:
            datum[i, j] = hFacC_3D[i, j, qbottom[i, j] + 1]
            bottom[i, j] = qbottom[i, j] + 1

# mask_values=[bottom,qbottom,surf,datum] #to loop over files?

# write mask values
os.chdir(ptd.mask_path)
filenames = [
    "mask_bottom_GoT_iNEST_V2",
    "mask_qbottom_GoT_iNEST_V2",
    "mask_surface_GoT_iNEST_V2",
    "mask_bhFacC_GoT_iNEST_V2",
]
surf = surf.copy(
    order="C"
)  # force to be C-contiguos, necessary rto save it with write, with np.save worked already
c = 0
for name in filenames:
    f = open(name, "wb")
    # f.write(mask_values[c,nx,:])
    if name == filenames[0]:
        f.write(bottom)
    if name == filenames[1]:
        f.write(qbottom)
    if name == filenames[2]:
        f.write(surf)
    if name == filenames[3]:
        f.write(datum)
    f.close()
print("done make_masks_GoT.m")
