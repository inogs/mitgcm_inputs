import argparse
from bitsea.utilities.argparse_types import existing_dir_path


def argument():
	parser = argparse.ArgumentParser(description = """
	Generates RBCs files
	- bottom_sources_S_mask.bin
	- bottom_sources_S_relaxation.bin
	- conc01_bottom_fluxes.bin
	- conc02_bottom_fluxes.bin
	- ...
	- check_fluxes.nc
    """)
	parser.add_argument(
        '--inputdir', '-i',
        type=existing_dir_path,
        required=True,
        help="Output of scarichi_json_gen.py "
    )
	parser.add_argument(
        '--domain', '-d',
        type=str,
        required=True,
        help ="NAD, SAD, ..."
    )
	parser.add_argument(
        '--domdir',
        type=existing_dir_path,
        required=True,
    )
	return parser.parse_args()

args = argument()


import numpy as np
import xarray as xr
import sys
import json

#

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

def open_old_rivers(
		json_file,
):
	print(json_file)
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

def fill_sal_conc(
		relax_sal,
		conc,
		x,
		y,
		z,
		depth,
		json_data,
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
		#k = np.argmin(np.abs(z.values - jd['CMS_depth']))
		print(x[i].values, y[j].values)
		print(z[k].values, depth[j,i].values)
		#while z[k] > depth[j,i]:
		#	k =- 1
		print('•••••••••••••••••••••••••••v')
		print(jd['Codice_Scarico'])
		print(z.values)
		print(depth[j,i].values)
		print('•••••••••••••••••••••••••••v')
		rel_S = jd['CMS_avgS'] - water_freshener
		if fixed_conc == 'None':
			c = jd['Carico_Ingresso_AE'] * jd['Dilution_factor']
		elif isinstance(fixed_conc, float):
			c = fixed_conc
		relax_sal[k,j,i] = rel_S
		conc_list[ns][j,i] = c
	return relax_sal, xr.where(relax_sal == 0., 0., 1.).astype('f4'), conc_list

#

def fill_river_conc(
		conc,
		x,
		y,
		z,
		depth,
		json_data,
		discharge_json,
		uniform_concentration = 1000.,
		fixed_conc = None,
):
	"""
	Given the info for the point sources and the spatial static data,
	fills the relaxation salinity array, creates the salinity mask and
	fills the concentration array.
	"""

	#n_sources = len(json_data)
	n_sources = 0
	for dr in discharge_json:
		if 'concentrations' in dr.keys():
			n_sources += 1
	conc_list = []
	print('\n•••••••••••••••\n',n_sources,'\n')
	for ns in range(n_sources):
		conc_list.append(conc * 0)

	riv_cntr = 0
	for ns, jd in enumerate(json_data): #[direc]
		print('\n')
		print(ns)
		print(riv_cntr)
		print('—————————————————————————————————\n')
		i = np.array(jd['longitude_indices']) #i = jd['longitude']
		if isinstance(i, str):
			rr = range(int(i[2:4]), int(i[-4:-1]))
			i = np.array([ii for ii in rr])
		j = np.array(jd['latitude_indices']) #j = jd['latitude']
		if jd['side'] == 'E':
			i -= 1
		elif jd['side'] == 'W':
			i += 1
		elif jd['side'] == 'S':
			j += 1
		elif jd['side'] == 'N':
			j -= 1
		k = np.argmin(np.abs(z.values - depth[j, i].values))
		print(j, i)
		try:
			cntr = 0
			while discharge_json[cntr]['name'] != jd['name']:
				cntr += 1
			if 'concentrations' in discharge_json[cntr].keys():
				if fixed_conc == None:
					c = uniform_concentration
				elif isinstance(fixed_conc, float):
					c = fixed_conc
				conc_list[riv_cntr][j,i] = c
				riv_cntr += 1
				print('\t\t', discharge_json[cntr]['name'], jd['name'])
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
		out_dir,
		conc_in_time=368,
):
	"""
	Args:
	- relax_sal : ndarray [depth,lat,lon] values to relax to
	- S_mask    : ndarray [depth,lat,lon] of 0.0 1.0 (float)
	- sst_mask  : ndarray [depth,lat,lon] of 0.0 1.0 (float)
	- concs     : list of arrays [lat,lon] with sources concentrations


	Given the relaxation salinity, salinity mask and concentration
	arrays, writes them to MITgcm-appropriate binary files:
	- bottom_sources_S_mask.bin
	- bottom_sources_S_relaxation.bin
	- SST_mask.bin
	- conc01_bottom_fluxes.bin
	- conc02_bottom_fluxes.bin
	- ...
	- check_fluxes.nc

	Conc files are replicated 'conc_in_time' times
	"""

	relax_sal.values.astype('f4').tofile(out_dir / 'bottom_sources_S_relaxation.bin')
	for i, conc in enumerate(concs):
		if conc_in_time == 1:
			conc.values.astype('f4').tofile(out_dir / f'conc{i+1:02}_bottom_fluxes.bin')
		else:
			np.stack([conc.values] * conc_in_time).astype('f4').tofile(
				out_dir / f'conc{i+1:02}_bottom_fluxes.bin')

	S_mask.values.astype('f4').tofile(out_dir / 'bottom_sources_S_mask.bin')
	sst_mask.values.astype('f4').tofile(out_dir / 'SST_mask.bin')

    #check
	xr.merge([S_mask.rename('salinity_mask'), relax_sal.rename('relax_salinity'), sst_mask.rename('sst_mask')] +
		[conc.rename(f'CONC{i+1:02}') for i, conc in enumerate(concs)]).to_netcdf(out_dir / "check_fluxes.nc")


#


###


base_path = args.inputdir
domain = args.domain

def main():
	inputfile = args.domdir / ('MIT_static_' + domain + '.nc')

	# sewage
	sewage_path = args.inputdir / f'PointSource_{domain}.json'
	sewers = open_point_sources(sewage_path)

	relax_salt = initialise_sal(xr.open_dataset(inputfile).hFacC)
	tracer_conc = initialise_conc_fluxes(xr.open_dataset(inputfile).hFacC)
	coords = get_spatial_masks(xr.open_dataset(inputfile))

	SST_mask = fill_sst_mask(xr.open_dataset(inputfile).hFacC)
	relax_salt, mask_salt, concentrations = fill_sal_conc(relax_salt, tracer_conc, coords['xc'], coords['yc'], coords['zc'], coords['depth'], sewers, fixed_conc=1.)

	# rivers
	rivers_path = args.domdir / f'rivers_positions_{domain}.json'
	old_rivers_path = args.inputdir / f'{domain}.json'
	rivers = open_river_sources(rivers_path)
	old_rivers = open_old_rivers(old_rivers_path)

	tracer_conc = initialise_conc_fluxes(xr.open_dataset(inputfile).hFacC)
	concentrations = concentrations + fill_river_conc(tracer_conc, coords['xc'], coords['yc'], coords['zc'], coords['depth'], rivers, old_rivers, fixed_conc=1.)
	print(type(concentrations), len(concentrations))

	write_binary_files(relax_salt, mask_salt, SST_mask, concentrations, out_dir=args.domdir)




#
#
if __name__ == '__main__':
	main()



