from aiida_openmx.input.structure import atom_spec_coord, atom_unit_vectors
from aiida_openmx.input.definition_of_atomic_species import atomic_species
from aiida_openmx.input.flat import replace_backslash
from aiida.orm import Dict
from pymatgen.core import Structure

def write_mixed_output(input_file, folder, data_backslash, structure):
    """
    Write key-value pairs and structure elements based on the specified sequence.
    :param output: Path to the output text file
    :param data: Dictionary with keys and their values
    :param structure_filename: List containing text to insert
    """
    data = replace_backslash(data_backslash)
    structure_string = {'Definition.of.Atomic.Species': atomic_species(structure, data['q']),
                        'Atoms.SpeciesAndCoordinates': atom_spec_coord(structure),
                        'Atoms.UnitVectors': atom_unit_vectors(structure)}
    with folder.open(input_file, 'w') as handle:
        for item in structure_string:
            # Write the corresponding structure item
            handle.write(f"{structure_string[item]}\n")
        for item in data:
                # Write the key-value pair from the dictionary
                value = data[item]
                handle.write(f"{item:<35} {value:<35}\n")