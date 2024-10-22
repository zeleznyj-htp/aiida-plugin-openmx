from pymatgen.io.cif import CifParser
from pymatgen.core import Structure

def cif_to_struct(filename):
    parser = CifParser(filename)
    structure = parser.parse_structures()[0]
    #print(type(structure))
    return structure

def struct_to_dict(structure):
    #print(type(structure))
    return structure.as_dict()

#print(cif_to_dict("Crys-MnO.cif")['charge'])

def atom_spec_coord(structure):
    string = "<Atoms.SpeciesAndCoordinates\n"
    for i in range(len(struct_to_dict(structure)['sites'])):
        a_vectors = [str(a) for a in struct_to_dict(structure)['sites'][i]['abc']]
        s_vectors = " ".join(a_vectors)
        string += (str(i+1) +"\t" + str(struct_to_dict(structure)['sites'][i]['species'][0]['element']) + "\t" + s_vectors + "\n")

    string += "Atoms.SpeciesAndCoordinates>"
    return string

def atom_unit_vectors(structure):
    string = "<Atoms.UnitVectors\n"
    for i in range(3):
        a_vectors = struct_to_dict(structure)['lattice']['matrix'][i]
        s_vectors = [str(a) for a in a_vectors]
        string += " ".join(s_vectors)+"\n"
    string += "Atoms.UnitVectors>"
    return string
