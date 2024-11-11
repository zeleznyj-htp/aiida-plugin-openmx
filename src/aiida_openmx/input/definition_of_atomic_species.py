import pandas as pd
import os
from pathlib import Path
import numpy as np
from aiida_openmx import data
from pymatgen.core import Structure

def get_elements(structure):
    elements = set()
    for i in range(len(structure['sites'])):
        elements.add(structure['sites'][i]['species'][0]['element'])
    return list(elements)

def pseudopotentials_filenames(elements):
    weird = ['Fe','Co','Ni','Cu','Zn']
    pseudopotentials_names = [atom + '_PBE19H' if atom in weird else atom + '_PBE19' for atom in elements]
    return pseudopotentials_names

def pseudo_basis_names(elements, q):
    path = Path(os.path.dirname(os.path.realpath(__file__)))
    csv_path = path.parent / 'data' / 'pseudopotentials.csv'
    df = pd.read_csv(csv_path)

    # Determine which column to use based on the value of q
    column_map = {1: 'Quick', 2: 'Standard', 3: 'Precise'}

    if q not in column_map:
        raise ValueError("Invalid value for q. Only q=1, q=2, and q=3 are supported.")

    # Create a dictionary to map values in the first column to values in the second column
    value_map = dict(zip(df['VPS'], df[column_map[q]]))
    # Match values in input_list with the first column and retrieve corresponding values from the second column
    output_list = [value_map.get(value, None) for value in pseudopotentials_filenames(elements)]
    return output_list

def valence_electrons(elements_on_site):
    # Load the CSV file into a pandas DataFrame
    path = Path(os.path.dirname(os.path.realpath(__file__)))
    csv_path = path.parent / 'data' / 'pseudopotentials.csv'
    df = pd.read_csv(csv_path)

    # Return the list of valences from the selected column
    element_to_value = dict(zip(df[df.columns[0]], df[df.columns[1]]))

    # Build the output list by looking up each element in the dictionary
    valence = [element_to_value.get(element, None) for element in pseudopotentials_filenames(elements_on_site)]
    return valence

def atomic_species(structure, q):
    string = "<Definition.of.Atomic.Species\n"
    elements = get_elements(structure)
    for i in range(len(elements)):
        string += (str(elements[i]) + "\t" + str(pseudo_basis_names(elements, q)[i]) + "\t" + str(pseudopotentials_filenames(elements)[i]) + "\n")
    string += "Definition.of.Atomic.Species>"
    return string

def band_path(n, critical_p, k_path):
    string = "<Band.kpath\n"
    for sub_path in k_path:
        for i in range(len(sub_path)-1):
            string += "  "+str(n)+"  "+" ".join([str(s) for s in critical_p[sub_path[i]]])+"   "+" ".join([str(s) for s in critical_p[sub_path[i+1]]])+"   "+sub_path[i]+" "+sub_path[i+1]+"\n"
    string += "Band.kpath>"
    return string

def band_kpath_unit_cell(structure):
    string = "<Band.KPath.UnitCell\n"
    lengths = structure['lattice']['matrix']
    s1 = 2*np.pi/(lengths[0][0])
    s2 = 2*np.pi/(lengths[1][1])
    s3 = 2*np.pi/(lengths[2][2])
    string += f"{s1} 0 0\n0 {s2} 0\n0 0 {s3}\n"
    string += "Band.KPath.UnitCell>\n"
    return string
