import argparse
import json
from pathlib import Path
import numpy as np
import xarray as xr


def create_bc_grid(path_results, path_grid, cbdry, grid_dims):
    """
    Input:
     path_results, directory where results are saved
     path_grid, path to MIT_static.nc file or directory containing it
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

    # Read netCDF file
    nc_file = Path(path_grid) / "MIT_static.nc" if Path(path_grid).is_dir() else Path(path_grid)
    ds = xr.open_dataset(nc_file)

    # Read latitude and longitude data
    YC = ds["YC"].values
    XC = ds["XC"].values
    lon, lat = np.meshgrid(XC, YC)

    if lat.shape != (ny, nx):
        raise ValueError(f"YC shape {lat.shape} != ({ny},{nx})")
    if debug:
        print("lat shape", lat.shape)
        print(lat)

    if lon.shape != (ny, nx):
        raise ValueError(f"XC shape {lon.shape} != ({ny},{nx})")
    if debug:
        print("lon shape", lon.shape)

    # Read drF data (vertical levels thickness)
    thk = ds["drF"].values

    print("drf shape", thk.shape)

    # Read HFacC data and calculate depth
    hfac = ds["hFacC"].values
    if hfac.shape != (nz, ny, nx):
        raise ValueError(f"hFacC shape {hfac.shape} != ({nz},{ny},{nx})")

    depth = np.sum(hfac*thk[:, None, None], axis=0)

    # Try to read angle data from netCDF, or use defaults
    try:
        AngleCS = ds["AngleCS"].values
    except KeyError:
        AngleCS = np.ones((ny, nx))
        print("AngleCS variable not found in netCDF. Generating a {ny,nx} ones array")

    try:
        AngleSN = ds["AngleSN"].values
    except KeyError:
        AngleSN = np.zeros((ny, nx))
        print("AngleSN variable not found in netCDF. Generating a {ny,nx} zeros array")

    # Close dataset
    ds.close()

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
    return


def main(path_results, path_grid, boundaries, grid_dims):
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
        create_bc_grid(path_results, path_grid, cbdry, grid_dims)
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
        help="path to MIT_static.nc file or directory containing it",
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

    args = parser.parse_args()
    main(args.path_results, args.path_grid, args.boundaries, args.grid_dims)