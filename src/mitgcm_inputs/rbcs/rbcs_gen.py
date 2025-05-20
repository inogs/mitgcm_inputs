import argparse
from bitsea.utilities.argparse_types import existing_dir_path


def argument():
	parser = argparse.ArgumentParser(description = """
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

def open_old_rivers(
		json_file,
):
	print(json_file)
	with open(json_file, 'r') as jfile:
		jdata = json.load(jfile)
	return jdata['river_points']


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
		only_one = False
):
	"""
	Given the info for the point sources and the spatial static data,
	fills the relaxation salinity array, creates the salinity mask and
	fills the concentration array.
	"""

	n_sources = len(json_data)
	conc_list = [conc] * n_sources

	for ns, jd in enumerate(json_data):
		lon = jd['Long']
		lat = jd['Lat']
		i = np.argmin(np.abs(x.values - lon))
		j = np.argmin(np.abs(y.values - lat))
		k = np.argmin(np.abs(z.values + depth[j,i].values))
		rel_S = jd['CMS_avgS'] - water_freshener
		c = jd['Carico_Ingresso_AE'] * jd['Dilution_factor']
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

	for ns, jd in enumerate(json_data): #[direc]
		print('\n')
		print(ns)
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
			while discharge_json[cntr]['rivername'] != jd['name']:
				cntr += 1
			print(cntr)
			print(discharge_json[cntr]['rivername'])
			c = discharge_json[cntr]['MEAN_2011_2023'] * uniform_concentration
			print(c)
		except:
			c = 50.e3
		#relax_sal[k, j, i] = rel_S
		conc_list[ns][j,i] += c
		print(conc_list[ns].max().values)
	print(np.unique(conc_list[1] == conc_list[0]))
	return conc_list #relax_sal, xr.where(relax_sal == 0., 0., 1.), conc

#

def write_binary_files(
		relax_sal,
		S_mask,
		concs,
		out_dir,
		aggregated=True,
		relax_path='bottom_sources_S_relaxation.bin',
		conc_path='tracer_concetrations.bin',
		mask_path='bottom_sources_S_mask.bin',
		aggreg_path='sewage_discharges.nc',
		conc_in_time=365,
):
	"""
	Given the relaxation salinity, salinity mask and concentration
	arrays, writes them to MITgcm-appropriate binary files.
	"""

	relax_sal.values.astype('f4').tofile(out_dir / relax_path)
	for i, conc in enumerate(concs):
		if conc_in_time == 1:
			conc.values.astype('f4').tofile(out_dir / (conc_path[:-4] + f'{i:02}' + conc_path[-4:]))
		else:
			np.stack([conc.values] * conc_in_time).astype('f4').tofile(
				out_dir / (conc_path[:-4] + f'{i:02}' + conc_path[-4:]))

	S_mask.values.astype('f4').tofile(out_dir / mask_path)

	if aggregated:
		xr.merge([S_mask.rename('salinity_mask'), relax_sal.rename('relax_salinity')] +
				 [conc.rename(f'CONC{i:02}') for i, conc in enumerate(concs)]).to_netcdf(out_dir / aggreg_path)


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

	relax_salt, mask_salt, concentrations = fill_sal_conc(relax_salt, tracer_conc, coords['xc'], coords['yc'], coords['zc'], coords['depth'], sewers)#, only_one=True)
	#write_binary_files(relax_salt, mask_salt, tracer_conc, out_dir=args.domdir, conc_path='conc01_bottom_fluxes.bin',)

	# rivers
	rivers_path = args.domdir / f'rivers_positions_{domain}.json'
	old_rivers_path = args.inputdir / f'RiverSource_{domain}.json'
	rivers = open_river_sources(rivers_path)
	old_rivers = open_old_rivers(old_rivers_path)

	tracer_conc = initialise_conc_fluxes(xr.open_dataset(inputfile).hFacC)
	concentrations = concentrations + fill_river_conc(tracer_conc, coords['xc'], coords['yc'], coords['zc'], coords['depth'], rivers, old_rivers) #, only_one=True
	print(type(concentrations), len(concentrations))

	write_binary_files(relax_salt, mask_salt, concentrations, out_dir=args.domdir, conc_path='CONC.bin',
					   aggreg_path='check_fluxes.nc')




#
#
if __name__ == '__main__':
	main()
	print('DONE!')


