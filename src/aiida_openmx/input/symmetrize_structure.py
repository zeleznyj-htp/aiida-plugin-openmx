import re
from pymatgen.core import Structure, Lattice


def parse_structure_update(text):
    # --- Extract new lattice vectors ---
    lattice_pattern = re.compile(
        r'a1\s*=\s*([\d\.\-Ee+]+)\s+([\d\.\-Ee+]+)\s+([\d\.\-Ee+]+).*?\n'
        r'\s*a2\s*=\s*([\d\.\-Ee+]+)\s+([\d\.\-Ee+]+)\s+([\d\.\-Ee+]+).*?\n'
        r'\s*a3\s*=\s*([\d\.\-Ee+]+)\s+([\d\.\-Ee+]+)\s+([\d\.\-Ee+]+)'
    )

    lattice_match = lattice_pattern.search(text)
    if not lattice_match:
        raise ValueError("Lattice vectors not found in text.")

    a1 = list(map(float, lattice_match.group(1, 2, 3)))
    a2 = list(map(float, lattice_match.group(4, 5, 6)))
    a3 = list(map(float, lattice_match.group(7, 8, 9)))
    lattice = [a1, a2, a3]

    # --- Extract fractional coordinates ---
    coord_pattern = re.compile(
        r'Fractional coordinates of the final structure.*?\n[*]+\n[*]+\n(.*?)\n\n',
        re.DOTALL
    )
    coord_block_match = coord_pattern.search(text)
    if not coord_block_match:
        raise ValueError("Fractional coordinates block not found.")

    coord_block = coord_block_match.group(1).strip()
    new_frac_coords = []
    for line in coord_block.splitlines():
        parts = line.strip().split()
        coords = list(map(float, parts[2:5]))  # skip index and element
        new_frac_coords.append(coords)

    return new_frac_coords, lattice


def update_structure_dict(structure_dict, new_frac_coords, new_lattice):
    # Update lattice info
    structure_dict['lattice']['matrix'] = new_lattice
    structure_dict['lattice']['a'] = new_lattice[0][0]
    structure_dict['lattice']['b'] = new_lattice[1][1]
    structure_dict['lattice']['c'] = new_lattice[2][2]
    structure_dict['lattice']['alpha'] = 90.0
    structure_dict['lattice']['beta'] = 90.0
    structure_dict['lattice']['gamma'] = 90.0
    # Optional: recompute volume if needed (or let pymatgen do it after loading)

    # Update fractional coordinates and recalculate cartesian
    from pymatgen.core import Lattice
    lattice = Lattice(new_lattice)

    for i, coord in enumerate(new_frac_coords):
        structure_dict['sites'][i]['abc'] = coord
        cart = lattice.get_cartesian_coords(coord)
        structure_dict['sites'][i]['xyz'] = list(cart)

    return structure_dict


def make_tetragonal(structure):
    structure_dict = structure.as_dict()
    # Step 1: Extract original lattice
    lattice_matrix = structure_dict['lattice']['matrix']
    a_vec, b_vec, c_vec = lattice_matrix
    a = a_vec[0]
    b = b_vec[1]
    c = c_vec[2]

    # Step 2: Create new lattice matrix with tetragonal symmetry
    m = (a + b) / 2
    new_lattice = [
        [m, 0, 0],
        [0, m, 0],
        [0, 0, c]
    ]

    # Step 3: Get original z-coordinates of first four atoms
    original_z = [site['abc'][2] for site in structure_dict['sites'][:4]]
    symmetrized_z = []
    symmetrized_z.append((original_z[0] + 1 - original_z[1]) / 2)
    symmetrized_z.append(1 - symmetrized_z[0])
    symmetrized_z.append((original_z[2] + 1 - original_z[3]) / 2)
    symmetrized_z.append(1 - symmetrized_z[2])

    # Step 4: Define new fractional coordinates
    new_frac_coords = [
        [0.0, 0.5, symmetrized_z[0]],
        [0.5, 1.0, symmetrized_z[1]],
        [1.0, 0.5, symmetrized_z[2]],
        [0.5, 0.0, symmetrized_z[3]],
        [0.5, 0.5, 0.0],
        [0.0, 0.0, 0.0],
    ]

    # Step 5: Update the dictionary using your method
    tetragonal_structure = update_structure_dict(structure_dict, new_frac_coords, new_lattice)
    return tetragonal_structure


def parse_and_change(text, structure):
    # Step 1: Your text input (from external source)
    #text = n.outputs.retrieved.get_object_content('aiida.out')
    # Step 2: Your original structure dict
    #structure_dict = structure.as_dict()
    # Step 3: Apply update
    new_frac_coords, new_lattice = parse_structure_update(text)
    updated_dict = update_structure_dict(structure, new_frac_coords, new_lattice)
    new_structure = Structure.from_dict(updated_dict)
    return new_structure


def structure_from_output(text, structure):
    new_structure = parse_and_change(text, structure)
    return make_tetragonal(new_structure)
