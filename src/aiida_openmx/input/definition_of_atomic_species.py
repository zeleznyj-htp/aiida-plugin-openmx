import pandas as pd
import importlib.resources as pkg_resources
from aiida_openmx import data
from pymatgen.core import Structure

def struct_to_dict(structure):
    return structure.as_dict()

def get_elements(structure):
    structure_dict = struct_to_dict(structure)
    elements = set()
    for i in range(len(structure_dict['sites'])):
        elements.add(structure_dict['sites'][i]['species'][0]['element'])
    return list(elements)

def get_elements_ordered(structure):
    structure_dict = struct_to_dict(structure)
    elements = []
    for i in range(len(structure_dict['sites'])):
        elements.append(structure_dict['sites'][i]['species'][0]['element'])
    return elements

def pseudopotentials_filenames(structure):
    pseudopotentials_names = [atom + '_PBE19' for atom in get_elements(structure)]
    return pseudopotentials_names

def pseudopotentials_filenames_ordered(structure):
    pseudopotentials_names = [atom + '_PBE19' for atom in get_elements_ordered(structure)]
    return pseudopotentials_names

def pseudo_basis_names(structure, csv_file, q):
    # Load the CSV file into a pandas DataFrame
    with pkg_resources.files(data).joinpath(csv_file).open('r') as f:
        df = pd.read_csv(f)

    # Determine which column to use based on the value of q
    column_map = {1: df.columns[2], 2: df.columns[3], 3: df.columns[4]}

    if q not in column_map:
        raise ValueError("Invalid value for q. Only q=1, q=2, and q=3 are supported.")

    # Select the appropriate column based on q
    selected_column = column_map[q]

    # Filter rows based on names_list and get the corresponding filenames
    filtered_df = df[df[df.columns[0]].isin(pseudopotentials_filenames_ordered(structure))]

    # Return the list of filenames from the selected column
    return filtered_df[selected_column].tolist()


#filenames = get_filenames_from_csv(csv_file, input(atoms_example), q)
#print(filenames)
def valence_electrons(elements_on_site, csv_file):
    # Load the CSV file into a pandas DataFrame
    with pkg_resources.files(data).joinpath(csv_file).open('r') as f:
        df = pd.read_csv(f)
    elements_name = [atom + '_PBE19' for atom in elements_on_site]

    # Return the list of valences from the selected column
    element_to_value = dict(zip(df[df.columns[0]], df[df.columns[1]]))

    # Build the output list by looking up each element in the dictionary
    valence = [element_to_value.get(element, None) for element in elements_name]
    return valence

def atomic_species(structure, csv_file, q):
    string = "<Definition.of.Atomic.Species\n"
    for i in range(len(get_elements(structure))):
        string += (str(get_elements(structure)[i]) + "\t" + str(pseudo_basis_names(structure, csv_file, q)[i]) + "\t" + str(pseudopotentials_filenames(structure)[i]) + "\n")
    string += "Definition.of.Atomic.Species>"
    return string
