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

def pseudopotentials_filenames(elements):
    pseudopotentials_names = [atom + '_PBE19' for atom in elements]
    return pseudopotentials_names

def pseudo_basis_names(elements, csv_file, q):
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
    filtered_df = df[df[df.columns[0]].isin(pseudopotentials_filenames(elements))]
    a = filtered_df[selected_column].tolist()
    # Return the list of filenames from the selected column
    return a


#filenames = get_filenames_from_csv(csv_file, input(atoms_example), q)
#print(filenames)
def valence_electrons(elements_on_site, csv_file):
    # Load the CSV file into a pandas DataFrame
    with pkg_resources.files(data).joinpath(csv_file).open('r') as f:
        df = pd.read_csv(f)

    # Return the list of valences from the selected column
    element_to_value = dict(zip(df[df.columns[0]], df[df.columns[1]]))

    # Build the output list by looking up each element in the dictionary
    valence = [element_to_value.get(element, None) for element in pseudopotentials_filenames(elements_on_site)]
    return valence

def atomic_species(structure, csv_file, q):
    string = "<Definition.of.Atomic.Species\n"
    elements = get_elements(structure)
    for i in range(len(elements)):
        string += (str(elements[i]) + "\t" + str(pseudo_basis_names(elements, csv_file, q)[i]) + "\t" + str(pseudopotentials_filenames(elements)[i]) + "\n")
    string += "Definition.of.Atomic.Species>"
    return string
