from pymatgen.io.cif import CifParser
def cif_to_dict(filename):
    parser = CifParser(filename)
    structure = parser.parse_structures()[0]
    return structure.as_dict()

#print(cif_to_dict("Crys-MnO.cif")['charge'])

def atom_spec_coord(filename):
    string = "<Atoms.SpeciesAndCoordinates\n"
    for i in range(len(cif_to_dict(filename)['sites'])):
        a_vectors = [str(a) for a in cif_to_dict(filename)['sites'][i]['abc']]
        s_vectors = " ".join(a_vectors)
        string += (str(i+1)+"\t"+str(cif_to_dict(filename)['sites'][i]['species'][0]['element'])+"\t"+s_vectors+"\n")

    string += "Atoms.SpeciesAndCoordinates>"
    return string

def atom_unit_vectors(filename):
    string = "<Atoms.UnitVectors\n"
    for i in range(3):
        a_vectors = cif_to_dict(filename)['lattice']['matrix'][i]
        s_vectors = [str(a) for a in a_vectors]
        string += " ".join(s_vectors)+"\n"
    string += "Atoms.UnitVectors>"
    return string
