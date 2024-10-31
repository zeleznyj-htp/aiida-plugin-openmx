from pymatgen.io.cif import CifParser
from pymatgen.core import Structure
import os
from aiida_openmx.input.definition_of_atomic_species import valence_electrons

def atom_spec_coord(structure):
    elements_on_site = []
    s_vectors = []
    for i in range(len(structure['sites'])):
        elements_on_site.append(str(structure['sites'][i]['species'][0]['element']))
        a_vectors = [str(a) for a in structure['sites'][i]['abc']]
        s_vectors.append(" ".join(a_vectors))
    valence = valence_electrons(elements_on_site)
    string = "<Atoms.SpeciesAndCoordinates\n"
    for i in range(len(elements_on_site)):
        string += (str(i+1) +"\t" + elements_on_site[i] + "\t" + s_vectors[i] +"\t"+ str(valence[i]) +"\n")
    string += "Atoms.SpeciesAndCoordinates>"
    return string

def atom_unit_vectors(structure):
    string = "<Atoms.UnitVectors\n"
    for i in range(3):
        a_vectors = structure['lattice']['matrix'][i]
        s_vectors = [str(a) for a in a_vectors]
        string += " ".join(s_vectors)+"\n"
    string += "Atoms.UnitVectors>"
    return string