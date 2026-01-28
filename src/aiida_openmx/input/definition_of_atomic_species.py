import pandas as pd
import os
from pathlib import Path
import numpy as np
from aiida_openmx import data
from pymatgen.core import Structure
import re

def get_elements(structure):
    elements = set()
    for i in range(len(structure['sites'])):
        elements.add(structure['sites'][i]['species'][0]['element'])
    return sorted(list(elements))

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
    n_path = 0
    string = "<Band.kpath\n"
    for sub_path in k_path:
        for i in range(len(sub_path)-1):
            n_path += 1
            string += "  "+str(n)+"  "+" ".join([str(s) for s in critical_p[sub_path[i]]])+"   "+" ".join([str(s) for s in critical_p[sub_path[i+1]]])+"   "+sub_path[i]+" "+sub_path[i+1]+"\n"
    string += "Band.kpath>\n"
    string += "Band.Nkpath  {}\n".format(n_path)
    return string

def band_kpath_unit_cell(unit_cell):
    string = "<Band.KPath.UnitCell\n"
    for i in range(3):
        string += ' '.join([str(x) for x in unit_cell[i, :]])
        string += '\n'
    string += "Band.KPath.UnitCell>\n"
    return string

def get_elements_order(structure):
    return [
        site["species"][0]["element"]
        for site in structure["sites"]
    ]

def unique_preserve_order(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

def basis_to_hubbard_line(basis_label: str, element: str, U_specs=None) -> str:
    """
    basis_label e.g. 'Ni6.0H-s2p2d2f1'
    U_specs: {'d': 5.0}  -> 1d 5.0, 2d 0.0
    """
    if U_specs is None:
        U_specs = {}

    norm_U_specs = {orb.lower(): float(U) for orb, U in U_specs.items()}

    try:
        basis = basis_label.split("-", 1)[1]
    except IndexError:
        raise ValueError(f"Invalid basis label: {basis_label!r}")

    matches = re.findall(r"([spdf])(\d+)", basis)
    if not matches:
        raise ValueError(f"No orbitals found in basis label: {basis_label!r}")

    shells = []
    for orb, count_str in matches:
        count = int(count_str)
        U = norm_U_specs.get(orb, 0.0)
        for n in range(1, count + 1):
            # Apply U only to the first shell of that l-channel
            U_here = U if n == 1 else 0.0
            shells.append(f"{n}{orb} {U_here:.1f}")

    return f"{element} " + " ".join(shells)

def generate_hubbard_block(structure: dict, q: int, hubbard_map: dict) -> str:
    elements_full = get_elements_order(structure)
    elements = unique_preserve_order(elements_full)
    basis_list = pseudo_basis_names(elements, q)

    lines = []
    for elem, basis in zip(elements, basis_list):
        U_specs = hubbard_map.get(elem, {})
        line = basis_to_hubbard_line(basis, elem, U_specs)
        lines.append("  " + line)

    block = "<Hubbard.U.values   # eV\n"
    block += "\n".join(lines) + "\n"
    block += "Hubbard.U.values>"

    return block