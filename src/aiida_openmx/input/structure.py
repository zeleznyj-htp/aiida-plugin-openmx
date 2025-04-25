from prompt_toolkit.input import Input
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure
import os
from aiida_openmx.input.definition_of_atomic_species import valence_electrons
from symmetr.magndata import get_magndata_structure
import numpy as np

def cartesian2spherical(x,y,z):
    r = np.linalg.norm([x,y,z])
    theta = np.arccos( z / r)
    phi = np.arctan2(y,x)
    return theta,phi

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

def atom_spec_coord(structure, spin_split=None, non_collinear=False, non_collinear_constraint=1, plusU_orbital=None):
    elements_on_site = []
    s_vectors = []
    for i in range(len(structure['sites'])):
        elements_on_site.append(str(structure['sites'][i]['species'][0]['element']))
        a_vectors = [str(a) for a in structure['sites'][i]['abc']]
        s_vectors.append(" ".join(a_vectors))

    bool_to_on = {True: 'on', False: 'off'}
    if plusU_orbital is not None:
        if isinstance(plusU_orbital,bool):
            plusU_orbital = [plusU_orbital] * len(elements_on_site)
        orbital_string = [bool_to_on[x] for x in plusU_orbital]
    else:
        orbital_string = [''] * len(elements_on_site)

    valence_split = get_valence_split(structure, spin_split)

    nc_string = [''] * 5
    if non_collinear:
        for i,site in enumerate((structure['sites'])):
            mom = site['properties']['magmom']
            theta, phi = cartesian2spherical(*mom)
            theta = theta * 180 / np.pi
            phi = phi * 180 / np.pi
            nc_string[i] += f"{theta:.4f} "
            nc_string[i] += f"{phi:.4f} "
            nc_string[i] += f"{theta:.4f} "
            nc_string[i] += f"{phi:.4f} "

        if isinstance(non_collinear_constraint,int):
            constraints = [non_collinear_constraint] * len(structure['sites'])
        else:
            constraints = non_collinear_constraint

        for i in range(len(structure['sites'])):
            nc_string[i] += str(constraints[i])

    string = "<Atoms.SpeciesAndCoordinates\n"
    for i in range(len(elements_on_site)):
        string += (str(i+1) +"\t" + elements_on_site[i] + "\t" + s_vectors[i] +"\t"+
                   str(valence_split[i][0]) + "\t"+ str(valence_split[i][1]) + "\t" + nc_string[i] + "\t" + orbital_string[i] + "\n")
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
