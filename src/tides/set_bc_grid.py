import argparse
import json
from pathlib import Path

import numpy as np
import scipy.io as sio

# I put the fuinctions below
# from readbin import readbin
# create_bc_grid


def create_bc_grid(path, dir_in, grid, cbdry):
    """
    Input:  dir_in, directory that has MITgcm grid files
     (following files are expected,
      YC.data. XC.data, hFacC.data
      RF.data, AngleCS.data, AngleSN.data )
      nx, first dimension of MITgcm grid, (west-east)
      ny, second  dimension of MITgcm grid (south-north)
      nz, dimension in z (or r)  direction
      cbdry, character for bdry to generate: s, n, w, e,
      (south, north, west, east)
    Output: No output
    The script will generate a grid file named bdry_grid_[snwe].mat
    """
    path_grid = path
    nx = grid[0]
    ny = grid[1]
    nz = grid[2]
    debug = False
    debug2 = True
    # Default location of boundary
    ny_north = ny - 1
    ny_south = 0
    nx_west = 0
    nx_east = ny - 1

    # Read latitude and longitude data
    lat = readbin(f"{dir_in}/YC.data", (ny, nx))
    if debug:
        print("lat shape", lat.shape)
        print(lat)
    lon = readbin(f"{dir_in}/XC.data", (ny, nx))
    if debug:
        print("lon shape", lon.shape)
    # Read RF data and calculate thickness
    rf = -readbin(f"{dir_in}/RF.data", (nz + 1))
    thk = np.diff(rf)  # thickness cells
    if debug:
        print("rf shape", rf.shape, "thk", thk.shape)
    # Read hFacC data and calculate depth
    hfac = readbin(f"{dir_in}/hFacC.data", (nz, ny, nx))
    thk3d = np.zeros_like(hfac)

    for k in range(nz):
        temp1 = hfac[k, :, :]
        thk3d[k, :, :] = thk[k] * temp1
    depth = np.sum(thk3d, axis=0)

    # Initialize AngleCS and AngleSN
    AngleCS = np.ones((ny, nx))
    AngleSN = np.zeros((ny, nx))

    try:
        AngleCS = readbin(f"{dir_in}/AnglesCS.data", (ny, nx))
    except FileNotFoundError:
        print("AngleCS file not found")
    try:
        AngleSN = readbin(f"{dir_in}/AnglesSN.data", (ny, nx))
    except FileNotFoundError:
        print("AngleSN file not found")

    # Determine boundary based on cbdry
    if cbdry == "s":
        lat_bc = lat[ny_south, :]
        lon_bc = lon[ny_south, :]
        depth_bc = depth[ny_south, :]
        AngleCS_bc = AngleCS[ny_south, :]
        AngleSN_bc = AngleSN[ny_south, :]
        print("Determining south boundary")
    elif cbdry == "n":
        lat_bc = lat[ny_north, :]
        lon_bc = lon[ny_north, :]
        depth_bc = depth[ny_north, :]
        AngleCS_bc = AngleCS[ny_north, :]
        AngleSN_bc = AngleSN[ny_north, :]
        print("Determining north boundary")
    elif cbdry == "w":
        lat_bc = lat[:, nx_west]
        lon_bc = lon[:, nx_west]
        depth_bc = depth[:, nx_west]
        AngleCS_bc = AngleCS[:, nx_west]
        AngleSN_bc = AngleSN[:, nx_west]
        print("Determining west boundary")
    elif cbdry == "e":
        lat_bc = lat[:, nx_east]
        lon_bc = lon[:, nx_east]
        depth_bc = depth[:, nx_east]
        AngleCS_bc = AngleCS[:, nx_east]
        AngleSN_bc = AngleSN[:, nx_east]
        print("Determining east boundary")
    else:
        print(f"No such selection: {cbdry}")
        return

    # Save boundary data to a .mat or python file
    mat = False
    if mat:
        bdry_file = f"bdry_grid_{cbdry}.mat"
        sio.savemat(
            bdry_file,
            {
                "lat_bc": lat_bc,
                "lon_bc": lon_bc,
                "depth_bc": depth_bc,
                "AngleCS_bc": AngleCS_bc,
                "AngleSN_bc": AngleSN_bc,
            },
        )
    else:
        bdry_file = f"bdry_grid_{cbdry}.npy"
        data = {
            "lat_bc": lat_bc,
            "lon_bc": lon_bc,
            "depth_bc": depth_bc,
            "AngleCS_bc": AngleCS_bc,
            "AngleSN_bc": AngleSN_bc,
        }
        np.save(path_grid + bdry_file, data)
        print(bdry_file, "file saved in ", path_grid)
        if debug2:
            loaded_array = np.load(
                path_grid + "bdry_grid_" + cbdry + ".npy", allow_pickle=True
            ).item()
            for name, arr in loaded_array.items():
                print("loaded_array:\n", name, arr.shape)
            # or_mat_file = sio.loadmat(path_grid + "bdry_grid_"+cbdry+".mat")
    return


def readbin(
    fnam, siz=(360, 224, 46), typ=1, prec="float32", skip=0, mform="<"
):
    """
    Function to read N-D binary field from a file.

    Parameters:
    fnam  : str  - Input path and file name
    siz   : tuple - Grid dimension (default (360, 224, 46))
    typ   : int  - 0 for sequential FORTRAN, 1 for plain binary (default)
    prec  : str  - Numeric precision (default 'float32')
    skip  : int  - Records to skip before reading (default 0)
    mform : str  - Machine format (default '<' for little-endian)

    Returns:
    fld   : np.ndarray - Output array of dimension siz

    Example usage:
    fld = readbin('data.bin', (360, 224, 46), 1, 'float32', 0, '>')
    """

    dtype_map = {
        "int8": np.int8,
        "int16": np.int16,
        "int32": np.int32,
        "int64": np.int64,
        "uint8": np.uint8,
        "uint16": np.uint16,
        "uint32": np.uint32,
        "uint64": np.uint64,
        "float32": np.float32,
        "float64": np.float64,  # doen't exist in python?
    }

    dtype = dtype_map.get(prec, np.float32)

    with open(fnam, "rb") as fid:
        if skip > 0:
            if typ == 0:
                for _ in range(skip):
                    _ = np.fromfile(
                        fid, dtype=dtype
                    )  # read_record() matlab function
            else:
                reclength = np.prod(siz) * dtype().itemsize
                fid.seek(skip * reclength, 0)

        if typ == 0:
            tmp = np.fromfile(
                fid, dtype=dtype
            )  # read_record() matlab function
        else:
            tmp = np.fromfile(fid, dtype=dtype, count=np.prod(siz))

    fld = tmp.reshape(siz)
    return fld


def main(grid, boundary, path, path_grid):
    # the script to generate bc grid files
    # which boundary
    bc_list = boundary
    # Grid
    nx = grid[0]
    ny = grid[1]
    nz = grid[2]
    print("Generating grid BC files, \n the grid size: ", nx, ny, nz)
    dir_in = path_grid
    print("From: \n ", dir_in)

    for cbdry in bc_list:
        create_bc_grid(path, dir_in, grid, cbdry)
        print(cbdry, "boundary created")


if __name__ == "__main__":
    # load default values from JSON (if exists)
    default_config = {}
    json_path = "tides_config.json"
    if Path(json_path).exists():
        with open(json_path) as f:
            default_config = json.load(f)
    # if present set up argparse that overrides JSON
    parser = argparse.ArgumentParser(
        description=(
            "Script for preparing the bc grid for the tidal boundary\
        condition for MITgcm (*.obcs)\ninputs can be provided by\
        json dictionary and/or by command line,ls\nthese overide json inputs."
        )
    )

    parser.add_argument(
        "--path",
        default=default_config.get(
            "path", "/g100_work/OGS23_PRACE_IT/apetroni/tides/"
        ),
        help=(
            "path to tide MITgcm grid boundary\
            (default: /g100_work/OGS23_PRACE_IT/apetroni/tides/)"
        ),
    )
    parser.add_argument(
        "--path_grid",
        default=default_config.get(
            "path_grid",
            "/g100_work/OGS23_PRACE_IT/apetroni/tides/\
                    grid/run_160p_AZAL_HR_for_grid",
        ),
        help="path to mitGCM grid (default:\
             '/g100_work/OGS23_PRACE_IT/apetroni/tides/\
             grid/run_160p_AZAL_HR_for_grid')",
    )
    parser.add_argument(
        "--boundary",
        nargs="+",  # accetta una lista di valori separati da spazio
        default=default_config.get("boundary", ["s"]),
        help="Constituent list (default: ['w'])",
        # Complete list boundaries ['n', 's', 'e', 'w']
    )
    parser.add_argument(
        "--grid",
        nargs="+",  # accetta una lista di valori separati da spazio
        default=default_config.get("grid", [320, 12, 50]),
        help="grid points nx, ny (default: [320,128,50])",
    )
    args = parser.parse_args()
    main(args.path, args.boundary, args.path, args.path_grid)
