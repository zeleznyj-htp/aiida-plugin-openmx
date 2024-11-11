from prompt_toolkit.input import Input
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure
import os
from aiida_openmx.input.definition_of_atomic_species import valence_electrons

def get_valence_split(structure,spin_split=None):
    elements_on_site = []
    for i in range(len(structure['sites'])):
        elements_on_site.append(str(structure['sites'][i]['species'][0]['element']))
    valence = valence_electrons(elements_on_site)

    if spin_split is None:
        spin_split = [0] * len(valence)

    valence_split = []
    for i, v in enumerate(valence):
        m = spin_split[i]
        n_up = (v + m) / 2
        n_down = (v - m) / 2

        if n_up < 0:
            n_up = 0
            n_down = v

        if n_down < 0:
            n_down = 0
            n_up = v

        valence_split.append(( n_up, n_down))

    return valence_split

def atom_spec_coord(structure, spin_split=None):
    elements_on_site = []
    s_vectors = []
    for i in range(len(structure['sites'])):
        elements_on_site.append(str(structure['sites'][i]['species'][0]['element']))
        a_vectors = [str(a) for a in structure['sites'][i]['abc']]
        s_vectors.append(" ".join(a_vectors))

    valence_split = get_valence_split(structure, spin_split)

    string = "<Atoms.SpeciesAndCoordinates\n"
    for i in range(len(elements_on_site)):
        string += (str(i+1) +"\t" + elements_on_site[i] + "\t" + s_vectors[i] +"\t"+
                   str(valence_split[i][0]) + "\t"+ str(valence_split[i][1])+"\n")
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