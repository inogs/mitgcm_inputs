import argparse
from bitsea.utilities.argparse_types import existing_dir_path
from bitsea.utilities.argparse_types import existing_file_path
from pathlib import Path
import numpy as np
import xarray as xr
import sys
import json

def argument():
	parser = argparse.ArgumentParser(description = """
	Generates RBCs files
    3D Files:
	- bottom_sources_S_mask.bin       mask with 1.0 at sewage point sources locations, 0.0 elsewhere
	- bottom_sources_S_relaxation.bin MDS clim salinity at sewage point sources locations, 0.0 elsewhere
    - SST_mask.bin   1.0 at surface, 0.0 elsewhere

	2D Files with zero values except at point sources locations:
	- conc01_bottom_fluxes.bin
	- conc02_bottom_fluxes.bin
	- ...

	A check file in netcdf format:
	- check_fluxes.nc

	A json file mapping tracers indexes to point sources:
	- RBCS_tracer_names.json
	{
    "TRAC52": "Scarico_Trieste_Servola",
    "TRAC53": "Scarico_Grado",
    "TRAC54": "Scarico_Depuratore_Di_Staran",
    "TRAC55": "Timavo",
    "TRAC56": "Isonzo"
	}

    Args:
        --sewage/-s : Json file output of scarichi_json_gen.py
		--river/-r  : Json file got from bathytools github repository
		--domdir    : Directory of the domain where
						- MIT_static.nc is expected
						- Rivers_positions.json is expected
						- output files are dumped


    """)
	parser.add_argument(
        '--sewage', '-s',
        type=existing_file_path,
        required=True,
        help="Json file output of scarichi_json_gen.py"
    )
	parser.add_argument(
        '--river', '-r',
        type=existing_file_path,
        required=True,
        help="Json file output got from bathytools github repository"
    )

	parser.add_argument(
        '--domdir',
        type=existing_dir_path,
        required=True,
    )
	return parser.parse_args()

args = argument()



def initialise_sal(
		hFac,
):
	"""
	Given the hFac mask from the static files generation,
	returns a zero-fill 3D array to be filled with
	salinity relaxation values for the bottom concentration
	point sources.
	"""

	return xr.zeros_like(hFac)

def initialise_conc_fluxes(
		hFac,
):
	"""
	Given the hFac mask from the static files generation,
	returns a zero-fill 2D array to be filled with
	concentration values for the bottom concentration
	point sources or surface river sources.
	"""

	return xr.zeros_like(hFac[0,:,:])

#

def get_spatial_masks(
		ds,
):
	"""
	Given the dataset from the static files generation,
	returns a zero-fill 2D array to be filled with
	concentration values for the bottom concentration
	point sources or surface river sources.
	"""

	depth = ds.Depth
	hFac = ds.hFacC
	xc = ds.XC
	yc = ds.YC
	zc = ds.Z

	return {'xc': xc, 'yc': yc, 'zc': zc, 'depth': depth, 'hFac': hFac}

def open_point_sources(
		json_file,
):
	"""
	Given the json of a specific domain's point sources,
	reads the data within.
	"""
	
	with open(json_file, 'r') as jfile:
		jdata = json.load(jfile)
	return jdata['discharge_points']

def open_river_sources(
		json_file,
):
	"""
	Given the json of a specific domain's river sources,
	reads the data within.
	"""

	with open(json_file, 'r') as jfile:
		jdata = json.load(jfile)
	return jdata
#

def open_sewage_rivers(
		json_file,
):
	print("reading" , json_file)
	with open(json_file, 'r') as jfile:
		jdata = json.load(jfile)
	return jdata['rivers']
#

def fill_sst_mask(
		hFac,
):
	"""
	Given a 3D field (hFacC from the static data),
	returns the mask for SST relaxation as a 3D array
	[depth, lat, lon] filled with 1s in the topmost level
	and 0s elsewhere
	"""

	SST_mask = xr.zeros_like(hFac)
	SST_mask[0,:,:] = 1.
	return SST_mask.astype('f4')
#

def get_opensea_swg_buoyant_plume(
		relax_sal = None,
		conc = None,
		x = None,
		y = None,
		z = None,
		depth = None,
		json_data = None,
		water_freshener = 0.5,
		only_one = False,
		fixed_conc = None,
):
	"""
	Given the info for the point sources and the spatial static data,
	fills the relaxation salinity array, creates the salinity mask and
	fills the concentration array.
	"""

	n_sources = len(json_data)
	conc_list = []
	for ns in range(n_sources):
		conc_list.append(conc * 0)
	
	Lon, Lat = np.meshgrid(x, y)
	Lon *= np.where(depth == 0., np.nan, 1)
	Lat *= np.where(depth == 0., np.nan, 1)

	for ns, jd in enumerate(json_data):
		lon = jd['Long']
		lat = jd['Lat']
		#i = np.argmin(np.abs(x.values - lon))
		#j = np.argmin(np.abs(y.values - lat))
		
		idx = np.nanargmin((Lon - lon)**2 + (Lat - lat)**2)
		j, i = idx//len(x.values), idx%len(x.values)
		k = np.argmin(np.abs(z.values - depth[j,i].values)) ### can we stop changing sign to the bathymetry?
		rel_S = jd['CMS_avgS'] - water_freshener

		if fixed_conc == 'None':
			c = jd['Carico_Ingresso_AE'] * jd['Dilution_factor']
		elif isinstance(fixed_conc, float):
			c = fixed_conc
		relax_sal[k,j,i] = rel_S
		conc_list[ns][j,i] = c
	return relax_sal, xr.where(relax_sal == 0., 0., 1.).astype('f4'), conc_list

#

def get_river_swg_plume(
		*,
		conc = None,
		x = None,
		y = None,
		z = None,
		depth = None,
		all_rivers = None,
		sewage_rivers = None,
		uniform_concentration = 1000.,
		fixed_conc = None,
):
	"""
	Given the info for the point sources and the spatial static data,
	fills the relaxation salinity array, creates the salinity mask and
	fills the concentration array.
	Returns:
	 conc_list: a list of concentration arrays, one per river source.
	"""

	#n_sources = len(json_data)
	n_sources = 0
	for dr in sewage_rivers:
		if 'concentrations' in dr.keys():
			n_sources += 1
	conc_list = []
	for ns in range(n_sources):
		conc_list.append(conc * 0)

	riv_cntr = 0
	for ns, jd in enumerate(all_rivers): #[direc]
		i = np.array(jd['longitude_indices']) #i = jd['longitude']
		if isinstance(i, str):
			rr = range(int(i[2:4]), int(i[-4:-1]))
			i = np.array([ii for ii in rr])
		j = np.array(jd['latitude_indices']) #j = jd['latitude']
		# shift indices of one cell downstream according to side of the river to avoid overwriting of
		# open boundary conditions (which are all zeros by construction for the E. coli concentrations)
		if jd['side'] == 'E':
			i -= 1
		elif jd['side'] == 'W':
			i += 1
		elif jd['side'] == 'S':
			j += 1
		elif jd['side'] == 'N':
			j -= 1
		k = np.argmin(np.abs(z.values - depth[j, i].values))
		try:
			cntr = 0
			while sewage_rivers[cntr]['name'] != jd['name']:
				cntr += 1
			if 'concentrations' in sewage_rivers[cntr].keys():
				if fixed_conc == None:
					c = uniform_concentration
				elif isinstance(fixed_conc, float):
					c = fixed_conc
				conc_list[riv_cntr][j,i] = c
				riv_cntr += 1
				#print('\t\t', discharge_json[cntr]['name'], jd['name'])
		except:
			print('Warning!')
		#relax_sal[k, j, i] = rel_S
		#conc_list[ns][j,i] = c
	return conc_list #relax_sal, xr.where(relax_sal == 0., 0., 1.), conc

#

def write_binary_files(
		relax_sal,
		S_mask,
		sst_mask,
		concs,
		out_dir : Path,
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
	relax_sal.values.astype('f4').tofile(filename)
	print(filename)

	for i, conc in enumerate(concs):
		filename = out_dir / f'conc{i+1:02}_bottom_fluxes.bin'
		print(filename)
		(conc.values).astype('f4').tofile(filename)

	filename = out_dir / "bottom_sources_S_mask.bin"
	S_mask.values.astype('f4').tofile(filename)
	print(filename)

	filename = out_dir / 'SST_mask.bin'
	sst_mask.values.astype('f4').tofile(filename)
	print(filename)

	# check
	filename  = out_dir / "check_fluxes.nc"
	print(filename)
	xr.merge([
		S_mask.rename('salinity_mask'),
		relax_sal.rename('relax_salinity'),
		sst_mask.rename('sst_mask')] +
		[conc.rename(f'CONC{i+1:02}') for i, conc in enumerate(concs)]
		).to_netcdf(filename)
#


###



if True: #def main():
	inputfile = args.domdir / ('MIT_static.nc')

	# sewage
	sewers_opensea = open_point_sources(args.sewage)

	relax_salt = initialise_sal(xr.open_dataset(inputfile).hFacC)
	tracer_conc = initialise_conc_fluxes(xr.open_dataset(inputfile).hFacC)
	coords = get_spatial_masks(xr.open_dataset(inputfile))

	SST_mask = fill_sst_mask(xr.open_dataset(inputfile).hFacC)
	relax_salt, mask_salt, opensea_sew_conc_list = get_opensea_swg_buoyant_plume(
		relax_sal=relax_salt,
		conc=tracer_conc,
		x=coords['xc'],
		y=coords['yc'],
        z=coords['zc'],
		depth=coords['depth'],
		json_data=sewers_opensea,
		fixed_conc=1.)

	print("len(opensea_sew_conc_list):", len(opensea_sew_conc_list))
	# rivers
	rivers_path = args.domdir / f'rivers_positions.json'
	all_rivers = open_river_sources(rivers_path)
	sewage_rivers = open_sewage_rivers(args.river)

	tracer_conc = initialise_conc_fluxes(xr.open_dataset(inputfile).hFacC)
	
	swg_river_conc_list=get_river_swg_plume(
		conc = tracer_conc,
		x=coords['xc'],
		y=coords['yc'],
		z=coords['zc'],
		depth=coords['depth'],
		all_rivers=all_rivers,
		sewage_rivers=sewage_rivers,
		fixed_conc=1.)
	print("len(swg_river_conc_list):", len(swg_river_conc_list))
	write_binary_files(relax_salt, mask_salt, SST_mask, opensea_sew_conc_list + swg_river_conc_list, out_dir=args.domdir)

	# build a dictionary with river concentrations for later use
	D={}
	i=51
	for s in sewers_opensea:
		i += 1
		name = 'TRAC' + f'{i:02}'
		D[name] = s['Nome_scarico']

	for r in all_rivers:
		i += 1
		name = 'TRAC' + f'{i:02}'
		D[name] = r['name']
    
	RBCS_tracer_names_file = args.domdir / 'RBCS_tracer_names.json'
	with open(RBCS_tracer_names_file, 'w') as jfile:
		print(RBCS_tracer_names_file)
		json.dump(D, jfile, indent=4)

#
#
#if __name__ == '__main__':
#	main()



