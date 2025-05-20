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
	if only_one:
		jd = json_data[0]
		lon = jd['Long']
		lat = jd['Lat']
		i = np.argmin(np.abs(x.values - lon))
		j = np.argmin(np.abs(y.values - lat))
		k = np.argmin(np.abs(z.values + depth[j, i].values))
		rel_S = jd['CMS_avgS'] - water_freshener
		c = jd['Carico_Ingresso_AE'] * jd['Dilution_factor']
		relax_sal[k, j, i] = rel_S
		conc[j, i] = c
	else:
		for jd in json_data:
			lon = jd['Long']
			lat = jd['Lat']
			i = np.argmin(np.abs(x.values - lon))
			j = np.argmin(np.abs(y.values - lat))
			k = np.argmin(np.abs(z.values + depth[j,i].values))
			rel_S = jd['CMS_avgS'] - water_freshener
			c = jd['Carico_Ingresso_AE'] * jd['Dilution_factor']
			relax_sal[k,j,i] = rel_S
			conc[j,i] = c
	return relax_sal, xr.where(relax_sal == 0., 0., 1.), conc

#

def write_binary_files(
		relax_sal,
		S_mask,
		conc,
		aggregated = True,
		out_dir = 'data/',
		relax_path = 'relax_salinity.bin',
		conc_path = 'tracer_concetrations.bin',
		mask_path = 'S_source_mask_MER_V3.dat',
		aggreg_path = 'sewage_discharges.nc',
		conc_in_time = 365,
):
	"""
	Given the the relaxation salinity, salinity mask and concentration
	arrays, writes them to MITgcm-appropriate binary files.
	"""
	
	relax_sal.values.astype('f4').tofile(out_dir / relax_path)
	if conc_in_time == 1:
		conc.values.astype('f4').tofile(out_dir / conc_path)
	else:
		np.stack([conc.values]*conc_in_time).astype('f4').tofile(out_dir / conc_path)
	S_mask.values.astype('f4').tofile(out_dir / mask_path)
	if aggregated:
		xr.merge([relax_sal.rename('relax_salinity'), conc.rename('tracer_conc')]).to_netcdf(out_dir / aggreg_path)

#

def fill_river_conc(
		conc,
		x,
		y,
		z,
		depth,
		json_data,
		water_freshener = 0.5,
		only_one = False,
):
	"""
	Given the info for the point sources and the spatial static data,
	fills the relaxation salinity array, creates the salinity mask and
	fills the concentration array.
	"""
	#TODO...
	if only_one:
		jd = json_data[0]
		lon = jd['Long']
		lat = jd['Lat']
		i = np.argmin(np.abs(x.values - lon))
		j = np.argmin(np.abs(y.values - lat))
		k = np.argmin(np.abs(z.values + depth[j, i].values))
		c = jd['Carico_Ingresso_AE'] * jd['Dilution_factor']
		relax_sal[k, j, i] = rel_S
		conc[j, i] = c
	else:
		for jd in json_data:
			lon = jd['Long']
			lat = jd['Lat']
			i = np.argmin(np.abs(x.values - lon))
			j = np.argmin(np.abs(y.values - lat))
			k = np.argmin(np.abs(z.values + depth[j,i].values))
			c = jd['Carico_Ingresso_AE'] * jd['Dilution_factor']
			relax_sal[k,j,i] = rel_S
			conc[j,i] = c
	return relax_sal, xr.where(relax_sal == 0., 0., 1.), conc

#


###


base_path = args.inputdir
domain = args.domain


inputfile=args.domdir / ('MIT_static_' + domain + '.nc')
relax_salt = initialise_sal(xr.open_dataset(inputfile).hFacC)
tracer_conc = initialise_conc_fluxes(xr.open_dataset(inputfile).hFacC)
coords = get_spatial_masks(xr.open_dataset(inputfile))

sewage_path = args.inputdir / f'PointSource_{domain}.json'
sewers = open_point_sources(sewage_path)

relax_salt, mask_salt, tracer_conc = fill_sal_conc(relax_salt, tracer_conc, coords['xc'], coords['yc'], coords['zc'], coords['depth'], sewers, only_one=True)
write_binary_files(relax_salt, mask_salt, tracer_conc, out_dir=args.domdir,
				   relax_path='bottom_sources_S_relaxation_'+domain+'.bin', conc_path='conc01_bottom_fluxes_'+domain+'.bin',
				   mask_path='bottom_sources_S_mask_'+domain+'.bin', aggreg_path='sewage_discharges_'+domain+'.bin')
#
#rivers_path = base_path + domain + f'/RiverSource_{domain}.json'
#rivers = open_river_sources(rivers_path)



#
#
# if __name__ == '__main__':
# 	main()
# 	print('DONE!')


