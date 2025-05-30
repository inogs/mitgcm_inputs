import argparse
from bitsea.utilities.argparse_types import existing_file_path
from bitsea.utilities.argparse_types import existing_dir_path



def argument():
    parser = argparse.ArgumentParser(description = """
    Reads seawage info on excel file and for all the domains
    build the json files with the coordinate of the pointsource for E.Coli tracers
    RiverSource_NAD.json
    RiverSource_SAD.json
    ...

    """)

    parser.add_argument(
        '--inputfile', '-i',
        type=existing_file_path,
        default="listaFIUMI_MER.xlsx",
        required=False,
        help="Fiumi, river discharge locations "
    )
    parser.add_argument(
        '--outdir', '-o',
        type=existing_dir_path,
        required=True,
    )


    return parser.parse_args()

args = argument()

import pandas as pd
import json

# Allegato 1 of sewage discharge locations 
excel_file = args.inputfile
OUTDIR = args.outdir

# select the n. of domain from 1 to 7 and name the json file

domains=['NAD','SAD','ION','SIC','TYR','LIG','SAR', 'GoT', 'GSN']

# Read the Excel file into a DataFrame
df = pd.read_excel(excel_file)

for id, namedomain in enumerate(domains):

    '''
    # Filter rows where the value in column 'Dominio' is id+1 (valid number from 1 to 7)
    if id < 7:
        filtered_df = df[df['Dominio'] == id+1]
    else:
        filtered_df = df[df['Sotto-Dominio'] == 8]  ### hardcoded horribly !!!
    '''
    filtered_df = df*1
    
    # Select only the interesting columns
    selected_columns = filtered_df[['rivername', 'lat_mouth', 'lon_mouth','MEAN_2011_2023', 'catchment','Region']]

    # Convert the selected rows to a list of dictionaries
    filtered_rows = selected_columns.to_dict(orient='records')

    # Create the final JSON structure
    output_data = {
        "file_name_origin": excel_file.name,
        "domain_number": id+1,
        "domain_name": domains[id],
        "n_points": len(filtered_rows),
        "river_points": filtered_rows
    }
    output_file = OUTDIR / ('RiverSource_' + domains[id] + '.json')
    # build the output file name
    # Write the final JSON structure to a JSON file
    with open(output_file, 'w') as json_file:
        json.dump(output_data, json_file, indent=4)

    print(f"JSON file created: {output_file}")

