import argparse
import json
from pathlib import Path
import numpy as np


def create_bc_grid(path_results, path_grid, cbdry, grid_dims, endian_val):
    """
    Input:
     path_results, directory where results are saved
     path_grid, directory of MITgcm grid files
     (following files are expected,
      YC.data. XC.data, hFacC.data
      RF.data, AngleCS.data, AngleSN.data )
      cbdry, character for bdry to generate: s, n, w, e,
      (south, north, west, east)
      grid_dims:
      nx, first dimension of MITgcm grid, (west-east)
      ny, second  dimension of MITgcm grid (south-north)
      nz, dimension in z (or r)  direction
    Output:
    The script will generate a grid file named bdry_grid_[snwe].npy and will save
    it in path_results
    """

    # path_results= "/home/sgaliana/ShareMER_tides/Resulting_NWMed_tide_files/"
    # path_results = Path(path_results).expanduser().resolve()
    # path_grid= "/home/sgaliana/ShareMER_tides/MITgcm_NWMed_grid/"
    # path_grid = Path(path_grid).expanduser().resolve()
    # cbdry="s"
    # grid_dims= [850,664,61]
    # endian_val="little"

    nx = grid_dims[0]
    ny = grid_dims[1]
    nz = grid_dims[2]
    debug = False
    debug2 = True
    # Default location of boundary
    ny_north = ny-1
    ny_south = 0
    nx_west = 0
    nx_east = nx-1
    cbdry = cbdry.lower() #Enforce lowercase boundary labels

    # Read latitude and longitude data
    lat = readbin(str(path_grid/ f"YC.data"), (ny, nx), endian=endian_val)
    if lat.shape != (ny,nx):
        raise ValueError(f"YC.data shape {lat.shape} != ({ny},{nx})")
    if debug:
        print("lat shape", lat.shape)
        print(lat)

    lon = readbin(str(path_grid / f"XC.data"), (ny, nx), endian=endian_val)
    if lon.shape != (ny,nx):
        raise ValueError(f"XC.data shape {lon.shape} != ({ny},{nx})")
    if debug:
        print("lon shape", lon.shape)

    # Read RF data and calculate thickness
    rf = -readbin(str(path_grid / f"RF.data"), (nz + 1), endian=endian_val)
    if rf.shape != (nz+1,):
        raise ValueError(f"RF.data shape {rf.shape} != ({nz+1})")

    thk = np.diff(rf)  # thickness cells
    if debug:
        print("rf shape", rf.shape, "thk", thk.shape)

    # Read hFacC data and calculate depth
    hfac = readbin(str(path_grid / f"hFacC.data"), (nz, ny, nx), endian=endian_val)
    if hfac.shape != (nz,ny,nx):
        raise ValueError(f"hFacC.data shape {hfac.shape} != ({nz},{ny},{nx})")


    # thk3d = np.zeros_like(hfac)
    # for k in range(nz):
    #     temp1 = hfac[k, :, :]
    #     thk3d[k, :, :] = thk[k] * temp1
    depth = np.sum(hfac*thk[:, None, None], axis=0)

    # # Initialize AngleCS and AngleSN
    # AngleCS = np.ones((ny, nx))
    # AngleSN = np.zeros((ny, nx))

    try:
        AngleCS = readbin(str(path_grid / f"AngleCS.data"), (ny, nx), endian=endian_val)
    except FileNotFoundError:
        AngleCS = np.ones((ny, nx))
        print("AngleCS file not found. Generating a {ny,nx} ones array")

    try:
        AngleSN = readbin(str(path_grid / f"AngleSN.data"), (ny, nx), endian=endian_val)
    except FileNotFoundError:
        AngleSN = np.zeros((ny, nx))
        print("AngleSN file not found. Generating a {ny,nx} zeros array")


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

    # Save boundary data to a python file

    bdry_file = f"bdry_grid_{cbdry}.npy"
    data = {
    "lat_bc": lat_bc,
    "lon_bc": lon_bc,
    "depth_bc": depth_bc,
    "AngleCS_bc": AngleCS_bc,
    "AngleSN_bc": AngleSN_bc,
    }
    np.save(path_results / bdry_file, data)
    print(bdry_file, "file saved in ", str(path_results))
    if debug2:
        loaded_array = np.load(
            path_results / bdry_file, allow_pickle=True
        ).item()
    for name, arr in loaded_array.items():
        print("loaded_array:\n", name, arr.shape)
    # or_mat_file = sio.loadmat(path_grid + "bdry_grid_"+cbdry+".mat")
    return


def readbin(fnam, siz=(850, 664, 61), typ=1, prec="float32", skip=0, endian="big"):
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
    fld = readbin('data.bin', (360, 224, 46), 1, 'float32', 0, 'big')
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
        "float64": np.float64,
    }

    dtype = dtype_map.get(prec, np.float32)
    endian_map = {"big": ">", "little": "<"}

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


def main(path_results, path_grid, boundaries, grid_dims, endian):
    # the script to generate bc grid files

    # Grid
    grid_dims = list(map(int, grid_dims))
    nx = grid_dims[0]
    ny = grid_dims[1]
    nz = grid_dims[2]
    print("Generating grid BC files, \n the grid size: ", nx, ny, nz)

    print("From: \n ", path_grid)

    path_results = Path(path_results).expanduser().resolve()
    path_grid = Path(path_grid).expanduser().resolve()

    path_results.mkdir(parents=True, exist_ok=True)

    for cbdry in boundaries:
        create_bc_grid(path_results, path_grid, cbdry, grid_dims, endian)
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
        "--path_results",
        default=default_config.get("path_results"),
        help="path where results will be saved",
    )
    parser.add_argument(
        "--path_grid",
        default=default_config.get("path_grid"),
        help="path to MITgcm grid",
    )
    parser.add_argument(
        "--boundaries",
        nargs="+",  # accetta una lista di valori separati da spazio
        default=default_config.get("boundaries"),
        help="Boundaries where tides are implemented",
        # Complete list boundaries ['n', 's', 'e', 'w']
    )
    parser.add_argument(
        "--grid_dims",
        nargs="+", # accetta una lista di valori separati da spazio
        default=default_config.get("grid_dims"),
        help="grid dimensions [nx, ny, nz]",
    )
    parser.add_argument(
        "--binary_data_endianess",
        default=default_config.get("binary_data_endianess"),
        help="MITgcm binary data endianess",
    )

    args = parser.parse_args()
    main(args.path_results, args.path_grid, args.boundaries, args.grid_dims, args.binary_data_endianess)
