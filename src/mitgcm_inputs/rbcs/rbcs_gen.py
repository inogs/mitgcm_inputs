import argparse
import json
from pathlib import Path

import numpy as np
import xarray as xr
from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import existing_file_path


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Generates RBCs files
    3D Files:
    - bottom_sources_S_mask.bin
           mask with 1.0 at sewage point sources locations, 0.0 elsewhere
    - bottom_sources_S_relaxation.bin MDS
           clim salinity at sewage point sources locations, 0.0 elsewhere
    - SST_mask.bin   1.0 at surface, 0.0 elsewhere

    2D Files with zero values except at point sources locations:
    - conc01_bottom_fluxes.bin
    - conc02_bottom_fluxes.bin
    - ...

    A check file in netcdf format:
    - check_fluxes.nc

    A text file RBCS_names.txt to append on mer config GSN_domain.json:
    ----------------------------------------
    ,
    "tracers": [
        {"id": 1, "name": "Scarico_Trieste_Servola"},
        {"id": 2, "name": "Scarico_Grado"},
        {"id": 3, "name": "Scarico_Depuratore_Di_Staran"},
        {"id": 4, "name": "Timavo"},
        {"id": 5, "name": "Isonzo"}
   ]
   ----------------------------------------

    Args:
        --sewage/-s : Json file output of scarichi_json_gen.py
        --river/-r  : Json file got from bathytools github repository
        --domdir    : Directory of the domain where
                        - MIT_static.nc is expected
                        - Rivers_positions.json is expected
                        - output files are dumped


    """
    )
    parser.add_argument(
        "--sewage",
        "-s",
        type=existing_file_path,
        required=True,
        help="Json file output of scarichi_json_gen.py",
    )
    parser.add_argument(
        "--river",
        "-r",
        type=existing_file_path,
        required=True,
        help="Json file output got from bathytools github repository",
    )

    parser.add_argument(
        "--domdir",
        type=existing_dir_path,
        required=True,
    )
    return parser.parse_args()


args = argument()


def get_spatial_masks(
    inputfile,
):
    """ """
    with xr.open_dataset(inputfile) as ds:
        depth = ds.Depth
        hFac = ds.hFacC
        xc = ds.XC
        yc = ds.YC
        zc = ds.Z

    return {"xc": xc, "yc": yc, "zc": zc, "depth": depth, "hFac": hFac}


def open_point_sources(
    json_file,
):
    """
    Given the json of a specific domain's point sources,
    reads the data within.
    """

    with open(json_file) as jfile:
        jdata = json.load(jfile)
    return jdata["discharge_points"]


def open_river_sources(
    json_file,
):
    """
    Given the json of a specific domain's river sources,
    reads the data within.
    """

    with open(json_file) as jfile:
        jdata = json.load(jfile)
    return jdata


#


def open_sewage_rivers(
    json_file,
):
    print("reading", json_file)
    with open(json_file) as jfile:
        jdata = json.load(jfile)
    return jdata["rivers"]


#


def get_opensea_swg_buoyant_plume(
    relax_sal=None,
    conc=None,
    x=None,
    y=None,
    z=None,
    depth=None,
    json_data=None,
    water_freshener=0.5,
    only_one=False,
    fixed_conc=None,
):
    """

    Returns:
     relax_sal: xr DataArray [depth,lat,lon] with salinity to relax to
     S_mask   : xr DataArray [depth,lat,lon] of 0.0 1.0 (float)
     conc_list: a list of concentration arrays, one per river source.
    """

    n_sources = len(json_data)
    conc_list = []
    for ns in range(n_sources):
        conc_list.append(conc * 0)

    Lon, Lat = np.meshgrid(x, y)
    Lon *= np.where(depth == 0.0, np.nan, 1)
    Lat *= np.where(depth == 0.0, np.nan, 1)

    for ns, jd in enumerate(json_data):
        lon = jd["Long"]
        lat = jd["Lat"]
        # i = np.argmin(np.abs(x.values - lon))
        # j = np.argmin(np.abs(y.values - lat))

        idx = np.nanargmin((Lon - lon) ** 2 + (Lat - lat) ** 2)
        j, i = idx // len(x.values), idx % len(x.values)
        k = np.argmin(
            np.abs(z.values - depth[j, i].values)
        )  # can we stop changing sign to the bathymetry?
        rel_S = jd["CMS_avgS"] - water_freshener

        if fixed_conc == "None":
            c = jd["Carico_Ingresso_AE"] * jd["Dilution_factor"]
        elif isinstance(fixed_conc, float):
            c = fixed_conc
        relax_sal[k, j, i] = rel_S
        conc_list[ns][j, i] = c
    S_mask = xr.where(relax_sal == 0.0, 0.0, 1.0).astype("f4")
    return relax_sal, S_mask, conc_list


#


def get_river_swg_plume(
    *,
    conc:np,
    rivers_positions:list[dict],
    sewage_rivers:list[dict],
    uniform_concentration=1000.0,
    fixed_conc:float=None,
)-> tuple[list,list]:
    """
    Args:
        conc: xr DataArray [lat,lon] with zero values
        rivers_positions: list of dict with river positions in the domain
        sewage_rivers: list of dict of Italy rivers with E.coli information
        uniform_concentration: float, concentration to assign to each river
        fixed_conc: if float, use this concentration for all rivers,
    Returns:
     conc_list: a list of concentration arrays, one per river source.
     ecoli_list: a list of river names corresponding to the conc_list
    """

    conc_list = []
    ecoli_list = []

    for jd in rivers_positions:
        i = np.array(jd["longitude_indices"])
        j = np.array(jd["latitude_indices"])
        if isinstance(i, str):
            rr = range(int(i[2:4]), int(i[-4:-1]))
            i = np.array([ii for ii in rr])

        # shift indices of one cell downstream according to side
        # of the river to avoid overwriting of
        # open boundary conditions (which are all zeros
        # by construction for the E. coli concentrations)
        if jd["side"] == "E":
            i -= 1
        elif jd["side"] == "W":
            i += 1
        elif jd["side"] == "S":
            j += 1
        elif jd["side"] == "N":
            j -= 1

            for jr in sewage_rivers:
                if jr["name"] == jd["name"]:
                    if "concentrations" in jr.keys():
                        ecoli_list.append(jr["name"])
                        river_conc = conc.copy()
                        if fixed_conc is None:
                            c = uniform_concentration
                        elif isinstance(fixed_conc, float):
                            c = fixed_conc
                        river_conc[j, i] = c
                        conc_list.append(river_conc)


    return conc_list, ecoli_list


#


def write_binary_files(
    relax_sal,
    S_mask,
    sst_mask,
    concs,
    out_dir: Path,
):
    """
    Args:
        relax_sal : ndarray [depth,lat,lon] values to relax to
        S_mask    : ndarray [depth,lat,lon] of 0.0 1.0 (float)
        sst_mask  : ndarray [depth,lat,lon] of 0.0 1.0 (float)
        concs     : list of xr DataArrays [lat,lon] with sources concentrations

    Given the relaxation salinity, salinity mask and concentration
    arrays, writes them to MITgcm-appropriate binary files:
            - bottom_sources_S_mask.bin
            - bottom_sources_S_relaxation.bin
            - SST_mask.bin
            - conc01_bottom_fluxes.bin
            - conc02_bottom_fluxes.bin
            - ...
            - check_fluxes.nc

    """

    filename = out_dir / "bottom_sources_S_relaxation.bin"
    relax_sal.values.astype("f4").tofile(filename)
    print(filename)

    for i, conc in enumerate(concs):
        filename = out_dir / f"conc{i + 1:02}_bottom_fluxes.bin"
        print(filename)
        (conc.values).astype("f4").tofile(filename)

    filename = out_dir / "bottom_sources_S_mask.bin"
    S_mask.values.astype("f4").tofile(filename)
    print(filename)

    filename = out_dir / "SST_mask.bin"
    sst_mask.values.astype("f4").tofile(filename)
    print(filename)

    # check
    filename = out_dir / "check_fluxes.nc"
    print(filename)
    xr.merge(
        [
            S_mask.rename("salinity_mask"),
            relax_sal.rename("relax_salinity"),
            sst_mask.rename("sst_mask"),
        ]
        + [conc.rename(f"CONC{i + 1:02}") for i, conc in enumerate(concs)]
    ).to_netcdf(filename)


#


###


if True:  # def main():
    inputfile = args.domdir / ("MIT_static.nc")

    # sewage
    sewers_opensea = open_point_sources(args.sewage)
    coords = get_spatial_masks(inputfile)
    DataArray3D = coords["hFac"]
    relax_salt = xr.zeros_like(DataArray3D)
    tracer_conc = xr.zeros_like(DataArray3D[0, :, :])

    SST_mask = xr.zeros_like(DataArray3D)
    SST_mask[0, :, :] = 1.0
    relax_salt, mask_salt, opensea_sew_conc_list = (
        get_opensea_swg_buoyant_plume(
            relax_sal=relax_salt,
            conc=tracer_conc,
            x=coords["xc"],
            y=coords["yc"],
            z=coords["zc"],
            depth=coords["depth"],
            json_data=sewers_opensea,
            fixed_conc=1.0,
        )
    )

    print("len(opensea_sew_conc_list):", len(opensea_sew_conc_list))
    # rivers
    
    domain_rivers = open_river_sources(args.domdir / "rivers_positions.json")
    italy_rivers = open_sewage_rivers(args.river)
    tracer_conc = xr.zeros_like(DataArray3D[0, :, :])

    swg_river_conc_list, ecoli_rivers = get_river_swg_plume(
        conc=tracer_conc,
        rivers_positions=domain_rivers,
        sewage_rivers=italy_rivers,
        fixed_conc=1.0,
    )
    print("len(swg_river_conc_list):", len(swg_river_conc_list))
    write_binary_files(
        relax_salt,
        mask_salt,
        SST_mask,
        opensea_sew_conc_list + swg_river_conc_list,
        out_dir=args.domdir,
    )

    sources=[s["Nome_scarico"] for s in sewers_opensea] + ecoli_rivers

    # text to insert in GSN_domain.json file,
    # just after 
    #     "model_parameters": {
    #     "OBCSsponge_N": ".FALSE.",
    #     "OBCSsponge_S": ".TRUE.",
    #     "OBCSsponge_E": ".FALSE.",
    #     "OBCSsponge_W": ".TRUE.",

    #     "useOBCSbalance": ".TRUE.",
    #     "OBCS_balanceFacN": "0.",
    #     "OBCS_balanceFacS": "1.",
    #     "OBCS_balanceFacE": "0.",
    #     "OBCS_balanceFacW": "1."
    # }
    point_source_output =  args.domdir / "RBCS_names.txt"
    with open(point_source_output, "w", encoding="utf-8") as f:
        f.write(',\n    "tracers": [\n')
        for i, s in enumerate(sources):
            comma = "," if i < len(sources) - 1 else ""
            point_source_dict = {"id": i+1, "name": s}
            f.write("        " + json.dumps(point_source_dict) + comma + "\n")
        f.write("   ]")
    
