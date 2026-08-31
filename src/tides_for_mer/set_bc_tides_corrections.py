import argparse
import json
import os
from datetime import datetime as dt
from pathlib import Path

import julian
import numpy as np
import pyTMD
from scipy.interpolate import interp1d


def create_tide_bc(
    path_results,
    tide_num,
    cbdry_used,
    grid_dims,
    time_corrections,
    endian_val,
):

    """ "
    Input: grid[0] the size of first dimension of MITgcm domain
    grid[1] the size of second dimension of MITgcm domain
    cbdry_used, cell array specify which boundary to add tidal
    boundary condition,
    {'s', 'n', 'w', 'e'}, South, North, West, and East
    If there is no ocean points on one boundary, do not specify it.
    Otherwise it will generate an error in interpolation.
    tide_num: tidal components to be included in tidal BC file
    time_corrections: time for nodal and MITgcm time reference corrections
    Output: Tidal boundary condition files will be generated with names like
        OB[N,S,E,W]_[u,v][Am,Ph].tide_obcs
    Note:
    1) The script will use  CAT2008b to generate tidal bc files
        with names like
        OB[NSWE]_[u,v]Am.tide_obcs, for tidal current amplitude   (m/s)
        OB[NSWE]_[u,v]Ph.tide_obcs, phase in seconds  (s)
    2) The following ten  tidal constituents can be included
              m2  s2  n2  k2  k1  o1  p1  q1 mf mm
              1   2   3   4   5   6   7   8  9  10
    3) The phase generated here inludes the Greenwich Phase of tidal
        constituent, nodal correction term, and a time constant in order to use
        nodal correction formula of CAT2008b
        phase_in_second = ( tidal_period*phase_in_radian/2*pi - mitgcm_timec)
    """

    # Dimension of Boundary files
    NX = grid_dims[0]
    NY = grid_dims[1]

    # Set the Boundary to add tides
    # AngleCS_bc, AngleSN_bc are the Cos and Sin and angle
    # of U axis of MITgcm
    # from due east

    for c_bdry in cbdry_used:

        bdry_file = path_results / f"bdry_grid_{c_bdry}.npy"

        data = np.load(bdry_file, allow_pickle=True).item()

        lat_bc = data["lat_bc"]
        lon_bc = data["lon_bc"]
        depth_bc = data["depth_bc"]
        AngleCS_bc = data["AngleCS_bc"]
        AngleSN_bc = data["AngleSN_bc"]

        boundary = c_bdry.upper()

        OBlength = NX if boundary in ["N", "S"] else NY
        print("\n --------------------------------------------------\n")
        print(
            "Applying corrections to static tidal forcing on boundary",
            c_bdry,
            "\n over a grid of lon size",
            len(lon_bc),
            "and lat size",
            len(lat_bc),
        )

        ## Load static tidal signal calculated with set_bc_tides_static

        u_am_st_all_file = path_results / f"data_u_am_st_all_{c_bdry}.npy"
        u_ph_st_all_deg_file = path_results / f"data_u_ph_st_all_deg_{c_bdry}.npy"
        v_am_st_all_file = path_results / f"data_v_am_st_all_{c_bdry}.npy"
        v_ph_st_all_deg_file = path_results / f"data_v_ph_st_all_deg_{c_bdry}.npy"

        u_amp_st = np.load(u_am_st_all_file)
        u_pha_st_deg = np.load(u_ph_st_all_deg_file)
        v_amp_st = np.load(v_am_st_all_file)
        v_pha_st_deg = np.load(v_ph_st_all_deg_file)

        ## Get nodal corrections for each constituent at time_corrections

        MJD = julian.to_jd(time_corrections, fmt="mjd")

        pu, pf, G = nodal(MJD, tide_num)

        print("Nodal corrections in radians:\n", pu, "\n", pf)


        ## Get constituent parametres for each constituent

        ph      = np.zeros(len(tide_num))
        omega   = np.zeros(len(tide_num))

        for i, c in enumerate(tide_num):
            amp_c, ph_c, omega_c, alpha_c, species_c = pyTMD.arguments._constituent_parameters(c)

            omega[i] = omega_c # rad/s - consituent angular frequency
            ph[i]    = ph_c # rad - constituent astronomical argument (reference epoch = 1 jan 1992)

        period = 2*np.pi/omega #seconds
        print(period)

        ## Loop over constituents and apply corrections

        # create out arrays for all boundary points and all constituents data
        L = NX if boundary in ['N','S'] else NY

        u_amp_out = np.zeros((L, len(tide_num)))
        u_pha_out_deg = np.zeros((L, len(tide_num)))
        u_pha_out_sec = np.zeros((L, len(tide_num)))
        v_amp_out = np.zeros((L, len(tide_num)))
        v_pha_out_deg = np.zeros((L, len(tide_num)))
        v_pha_out_sec = np.zeros((L, len(tide_num)))


        print("Applying corrections for each constitutent:")
        print("Nodal (time dependent), astronomical argument and MITgcm reference time (time dependent)")
        for i, c in enumerate(tide_num):

            ua_c = u_amp_st[:,i]
            va_c = v_amp_st[:,i]
            up_c = u_pha_st_deg[:,i]
            vp_c = v_pha_st_deg[:,i]

            ## Apply nodal corrections and astronomical argument correction
            # Tidal prediction in MItgcm will be:
            # pf*amp*cos[omega(t-time0+time0-time1992) -(pha - ph - pu)]
            # so the phase that we will give to MItgcm will be:
            # phi_MITgcm = (omega*(time0-time1992)-(pha - ph - pu))/omega, because MITgcm wants it in seconds

            ua_c *= pf[0,i] #dims of pf are= 1,len(tide_num), 1 because we only calculate nodal corrections for one time
            va_c *= pf[0,i]
            up_c = np.deg2rad(up_c) - ph[i] - pu[0,i] #nodal correction and astronomical argument correction (pha - ph - pu)
            vp_c = np.deg2rad(vp_c) - ph[i] - pu[0,i]

            up_c = np.rad2deg(up_c)
            vp_c = np.rad2deg(vp_c)

            # Transform from transports to velocities by dividing by the depth, but set to 0 wherever depth_bc == 0
            ua_c = np.divide(ua_c, depth_bc, out=np.zeros_like(ua_c), where=(depth_bc != 0))
            va_c = np.divide(va_c, depth_bc, out=np.zeros_like(va_c), where=(depth_bc != 0))

            # Rotate into model grid
            ua_cr, up_cr, va_cr, vp_cr = tide_rot(ua_c, up_c, va_c, vp_c, AngleCS_bc, AngleSN_bc)

            ## Apply MITgcm reference time correction and transform phases to seconds for MITgcm
            epoch = dt(1992, 1, 1)
            t0_mitgcm = time_corrections
            deltat_seconds = (t0_mitgcm-epoch).total_seconds()

            u_amp_out[:, i] = ua_cr
            u_pha_out_deg[:,i] = ((up_cr - (np.rad2deg((2 * np.pi * deltat_seconds)/period[i]))% 360.0 +180) % 360.0)-180 #The % 360.0 is to wrap it to [0,360)
            u_pha_out_sec[:, i] = ((period[i] * np.deg2rad(up_cr)) / (2 * np.pi) - deltat_seconds) % period[i] #The % period[i] Wrap to modulo period[i]

            v_amp_out[:, i] = va_cr
            v_pha_out_deg[:,i] = ((vp_cr- (np.rad2deg((2 * np.pi * deltat_seconds)/period[i])) % 360.0 +180) % 360.0)-180
            v_pha_out_sec[:, i] = ((period[i] * np.deg2rad(vp_cr)) / (2 * np.pi) - deltat_seconds) % period[i]


        u_amp_out[np.isnan(u_amp_out)] = 0
        u_pha_out_deg[np.isnan(u_pha_out_deg)] = 0
        u_pha_out_sec[np.isnan(u_pha_out_sec)] = 0

        v_amp_out[np.isnan(v_amp_out)] = 0
        v_pha_out_deg[np.isnan(v_pha_out_deg)] = 0
        v_pha_out_sec[np.isnan(v_pha_out_sec)] = 0

        # Create a boolean mask where depth is exactly 0
        land_mask = (depth_bc == 0)

        # Zero out all tidal constituents for these land points
        u_amp_out[land_mask, :] = 0
        u_pha_out_deg[land_mask, :] = 0
        u_pha_out_sec[land_mask, :] = 0
        v_amp_out[land_mask, :] = 0
        v_pha_out_deg[land_mask, :] = 0
        v_pha_out_sec[land_mask, :] = 0

        print("Saving results for each constitutent")

        fname_u_am_all = "data_u_am_all_" + c_bdry + ".npy"
        fname_u_ph_all = "data_u_ph_all_" + c_bdry + ".npy"
        fname_u_ph_all_deg = "data_u_ph_all_deg_" + c_bdry + ".npy"
        fname_v_am_all = "data_v_am_all_" + c_bdry + ".npy"
        fname_v_ph_all = "data_v_ph_all_" + c_bdry + ".npy"
        fname_v_ph_all_deg = "data_v_ph_all_deg_" + c_bdry + ".npy"

        np.save(path_results / fname_u_am_all, u_amp_out)
        np.save(path_results / fname_u_ph_all, u_pha_out_sec)
        np.save(path_results / fname_u_ph_all_deg, u_pha_out_deg)
        np.save(path_results / fname_v_am_all, v_amp_out)
        np.save(path_results / fname_v_ph_all, v_pha_out_sec)
        np.save(path_results / fname_v_ph_all_deg, v_pha_out_deg)

        vel = ["u", "v"]
        fld = ["Am", "Ph"]

        cfld = fld[0]

        cvel = vel[0]
        fnm = path_results / f"OB{boundary}_{cvel}{cfld}.tide_obcs"
        writebin(str(fnm), u_amp_out.T, endian=endian_val)

        cvel = vel[1]
        fnm = path_results / f"OB{boundary}_{cvel}{cfld}.tide_obcs"
        writebin(str(fnm), v_amp_out.T, endian=endian_val)


        cfld = fld[1]

        cvel = vel[0]
        fnm = path_results / f"OB{boundary}_{cvel}{cfld}.tide_obcs"
        writebin(str(fnm), u_pha_out_sec.T, endian=endian_val)

        cvel = vel[1]
        fnm = path_results / f"OB{boundary}_{cvel}{cfld}.tide_obcs"
        writebin(str(fnm), v_pha_out_sec.T, endian=endian_val)

    return True

def nodal(MJD, conlist):
    """
    function to compute nodal corrections for tidal constituents
    USAGE:
    [pu,pf,G]=nodal(MJD,tide_num)
    PARAMETERS:
    MJD: Modified Julian Date
    tide_num: list of tidal constituents
    OUTPUT:
    pu, pf: nodal corrections for the constituents
    G: phase correction in degrees
    """
    pu, pf, G = pyTMD.arguments.arguments(
        MJD, conlist, corrections='OTIS')  # pu in radians and G in degrees
    return pu, pf, G

def tide_rot(ua_c, up_c, va_c, vp_c, angleCS, angleSN):
    """
    Rotate tidal current parameter from NS and EW to a new coordinate
    defined by angleCS and angleSN
    Input:
       ua_c, up_c(in deg)
       va_c, vp_c(in deg)
       AngleCS  angleSN:  cos(alpha), sin(alpha);
       alpha is the angle between new coordinate (ua_cr, va_cr)
       and due east
    Output:
        ua_cr, up_cr(in deg)
        va_cr, vp_cr(in deg);
    Both input and output are assumed (1 L) array "
    """
    # Convert degrees to radians
    up_c = np.deg2rad(up_c)
    vp_c = np.deg2rad(vp_c)

    # Calculate the new u-component amplitude and phase
    u1_p1 = ua_c * angleCS * np.cos(
        up_c
    ) + va_c * angleSN * np.cos(vp_c)
    u1_p2 = ua_c* angleCS * np.sin(
        up_c
    ) + va_c * angleSN * np.sin(vp_c)

    ua_cr = np.sqrt(u1_p1**2 + u1_p2**2)

    # sin1 = u1_p2 / ua_cr
    # cos1 = u1_p1 / ua_cr
    # up_cr = np.arctan2(sin1, cos1)
    # up_cr = np.rad2deg(up_cr)
    up_cr = np.zeros_like(ua_cr)
    mask = ua_cr > 0 # To ensure that we do not divide by 0
    up_cr[mask] = np.rad2deg(np.arctan2(u1_p2[mask], u1_p1[mask]))

    # Calculate the new v-component amplitude and phase
    v1_p1 = -ua_c * angleSN * np.cos(
        up_c
    ) + va_c * angleCS * np.cos(vp_c)
    v1_p2 = -ua_c * angleSN * np.sin(
        up_c
    ) + va_c * angleCS * np.sin(vp_c)
    va_cr = np.sqrt(v1_p1**2 + v1_p2**2)
    # sin1 = v1_p2 / va_cr
    # cos1 = v1_p1 / va_cr
    # vp_cr = np.arctan2(sin1, cos1)
    # vp_cr = np.rad2deg(vp_cr)
    vp_cr = np.zeros_like(va_cr)
    mask = va_cr > 0 # To ensure that we do not divide by 0
    vp_cr[mask] = np.rad2deg(np.arctan2(v1_p2[mask], v1_p1[mask]))

    # Return the new amplitudes and phases
    return ua_cr, up_cr, va_cr, vp_cr

def writebin(fnam, fld, prec="float32", endian="big", skip=0, typ=1):
    """
    Write N-D binary field to a file.

    Parameters
    ----------
    fnam : str
        Output file name
    fld : ndarray
        Data to write
    prec : str
        Precision ('float32', 'float64', ...)
    endian : str
        'big' or 'little'
    skip : int
        Number of records to skip
    typ : int
        1 = plain binary (MITgcm style)
        0 = sequential FORTRAN (NOT implemented)
    """

    if typ != 1:
        raise NotImplementedError("Sequential FORTRAN (typ=0) not supported")

    endian_map = {"big": ">", "little": "<"}
    prec_map = {
        "float32": "f4",
        "float64": "f8",
        "int16": "i2",
        "int32": "i4",
        "int64": "i8",
    }

    if endian not in endian_map:
        raise ValueError("endian must be 'big' or 'little'")
    if prec not in prec_map:
        raise ValueError(f"Unsupported precision: {prec}")

    dtype = np.dtype(endian_map[endian] + prec_map[prec])
    fld = np.asarray(fld, dtype=dtype)

    mode = "r+b" if os.path.exists(fnam) else "wb"

    with open(fnam, mode) as f:
        if skip > 0:
            f.seek(skip * fld.size * dtype.itemsize, 0)
        f.write(fld.tobytes())


def main(path_results, constituents, boundaries, grid_dims, corrections_datetime, endianess):
    print("\n-----------------------------------------------------")
    print(f"Results path:          {path_results}")
    print(f"Constituents:          {constituents}")
    print(f"Boundaries:            {boundaries}")
    print(f"Grid dimensions:       {grid_dims}")
    print(f"Datetime for time dependent corrections: {corrections_datetime}")
    print(f"MITgcm data endianess: {endianess}")
    print("\n-----------------------------------------------------")
    print("\nCreating MITgcm Tidal BC")

    # Time for time dependent corrections

    time_corrections = dt.strptime(corrections_datetime, "%Y%m%d-%H:%M:%S")
    print("Corrections time", time_corrections, julian.to_jd(time_corrections, fmt="jd"))

    grid_dims = list(map(int, grid_dims))

    path_results = Path(path_results).expanduser().resolve()


    print("Domain sides", boundaries, "with constituents", constituents)
    create_tide_bc(
        path_results,
        constituents,
        boundaries,
        grid_dims,
        time_corrections,
        endianess,
    )


if __name__ == "__main__":
    # load default values from JSON (if exists)
    default_config = {}
    json_path = "tides_config_new.json"
    if Path(json_path).exists():
        with open(json_path) as f:
            default_config = json.load(f)
    # if present set up argparse that overrides JSON
    parser = argparse.ArgumentParser(
        description="Script for creating tidal boundary condition\
        for MITgcm (*.obcs)\ninputs can be provided by json dictionary\
        and/or by command line,ls\nthese overide json inputs."
    )
    subparsers = parser.add_subparsers(title="subcommand", dest="subcommand")
    parser.add_argument(
        "--path_results",
        default=default_config.get("path_results"),
        help=(
            "path to the folder where results are saved and where the results of\
            set_bc_grid and set_bc_tides_static have been saved"
        ),
    )
    parser.add_argument(
        "--constituents",
        nargs="+",  # accetta una lista di valori separati da spazio
        default=default_config.get("constituents"),
        help=("Constituent list"),
    )
    parser.add_argument(
        "--boundaries",
        nargs="+",  # accetta una lista di valori separati da spazio
        default=default_config.get("boundaries"),
        help=("Boundary list where tides have to be added"),
    )
    parser.add_argument(
        "--grid_dims",
        nargs="+",  # accetta una lista di valori separati da spazio
        default=default_config.get("grid_dims"),
        help="number of grid points in each dimension (longitude, latitude, depth): nx, ny, nz",
    )
    parser.add_argument(
        "--corrections_datetime",
        default=default_config.get("corrections_datetime"),
        help="Datetime for time dependent corrections calculations (yymmdd-hh:mm:ss)",
    )
    parser.add_argument(
        "--binary_data_endianess",
        default=default_config.get("binary_data_endianess"),
        help="MITgcm binary data endianess",
    )


    args = parser.parse_args()

    main(
        args.path_results,
        args.constituents,
        args.boundaries,
        args.grid_dims,
        args.corrections_datetime,
        args.binary_data_endianess,
    )