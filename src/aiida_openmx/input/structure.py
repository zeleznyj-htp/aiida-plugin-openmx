from pymatgen.io.cif import CifParser
from pymatgen.core import Structure
import os
from aiida_openmx.input.definition_of_atomic_species import valence_electrons

'''
def file_to_struct(filename):
    """
    Load a structure from a .cif or .json file.

    Parameters:
    - file_path (str): Path to the structure file (.cif or .json)

    Returns:
    - Structure: A pymatgen Structure object
    """
    # Get the file extension
    file_ext = os.path.splitext(filename)[1].lower()

    if file_ext == ".cif":
        structure = cif_to_struct(filename)
    elif file_ext == ".json":
        structure = Structure.from_file(filename)
    else:
        raise ValueError("Unsupported file format. Please use a .cif or .json file.")

    return structure

def cif_to_struct(filename):
    parser = CifParser(filename)
    structure = parser.parse_structures()[0]
    #print(type(structure))
    return structure
'''
def struct_to_dict(structure):
    #print(type(structure))
    return structure.as_dict()

def dict_to_struct(structure_dict):
    loaded_structure = Structure.from_dict(structure_dict)
    return loaded_structure

#print(cif_to_dict("Crys-MnO.cif")['charge'])


def atom_spec_coord(structure, csv_file):
    elements_on_site = []
    s_vectors = []
    for i in range(len(struct_to_dict(structure)['sites'])):
        elements_on_site.append(str(struct_to_dict(structure)['sites'][i]['species'][0]['element']))
        a_vectors = [str(a) for a in struct_to_dict(structure)['sites'][i]['abc']]
        s_vectors.append(" ".join(a_vectors))
    valence = valence_electrons(elements_on_site, csv_file)
    string = "<Atoms.SpeciesAndCoordinates\n"
    for i in range(len(elements_on_site)):
        string += (str(i+1) +"\t" + elements_on_site[i] + "\t" + s_vectors[i] +"\t"+ str(valence[i]) +"\n")
    string += "Atoms.SpeciesAndCoordinates>"
    return string

'''
def atom_spec_coord(structure):
    string = "<Atoms.SpeciesAndCoordinates\n"
    for i in range(len(struct_to_dict(structure)['sites'])):
        elements_on_site = struct_to_dict(structure)['sites'][i]['species'][0]['element']
        a_vectors = [str(a) for a in struct_to_dict(structure)['sites'][i]['abc']]
        s_vectors = " ".join(a_vectors)
        string += (str(i+1) +"\t" + str(struct_to_dict(structure)['sites'][i]['species'][0]['element']) + "\t" + s_vectors + "\n")
    string += "Atoms.SpeciesAndCoordinates>"
    return string
'''
def atom_unit_vectors(structure):
    string = "<Atoms.UnitVectors\n"
    for i in range(3):
        a_vectors = struct_to_dict(structure)['lattice']['matrix'][i]
        s_vectors = [str(a) for a in a_vectors]
        string += " ".join(s_vectors)+"\n"
    string += "Atoms.UnitVectors>"
    return string
