import pandas as pd
from structure import cif_to_dict

def pseudopotentials_filenames(structure_file):
    pseudopotentials_names = [atom + '_PBE19' for atom in get_elements(structure_file)]
    return pseudopotentials_names

def pseudo_basis_names(structure_file, csv_file, q):
    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file)

    # Determine which column to use based on the value of q
    column_map = {1: df.columns[2], 2: df.columns[3], 3: df.columns[4]}

    if q not in column_map:
        raise ValueError("Invalid value for q. Only q=1, q=2, and q=3 are supported.")

    # Select the appropriate column based on q
    selected_column = column_map[q]

    # Filter rows based on names_list and get the corresponding filenames
    filtered_df = df[df[df.columns[0]].isin(pseudopotentials_filenames(structure_file))]

    # Return the list of filenames from the selected column
    return filtered_df[selected_column].tolist()


#filenames = get_filenames_from_csv(csv_file, input(atoms_example), q)
#print(filenames)

def get_elements(structure_file):
    structure = cif_to_dict(structure_file)
    elements = set()
    for i in range(len(structure['sites'])):
        elements.add(structure['sites'][i]['species'][0]['element'])
    return list(elements)

def atomic_species(structure_file, csv_file, q):
    string = "<Definition.of.Atomic.Species\n"
    for i in range(len(get_elements(structure_file))):
        string += (str(get_elements(structure_file)[i])+"\t"+str(pseudo_basis_names(structure_file, csv_file, q)[i])+"\t"+str(pseudopotentials_filenames(structure_file)[i])+"\n")
    string += "Definition.of.Atomic.Species>"
    return string

# Example usage:
csv_file = 'pseudopotentials.csv'  # Replace with your actual CSV file path
q = 2  # Example q value

#print(atomic_species("Crys-MnO.cif", csv_file, q))

