from aiida_openmx.input.structure import atom_spec_coord, atom_unit_vectors, file_to_struct
from aiida_openmx.input.definition_of_atomic_species import atomic_species
from aiida_openmx.input.flat import replace_backslash
from aiida.orm import Dict
from pymatgen.core import Structure


#structure_filename = 'Methane.cif'
#output_file = 'output.txt'
#csv_file = 'pseudopotentials.csv'  # Replace with your actual CSV file path

def write_mixed_output(input_file, folder, data_backslash, structure_filename):
    """
    Write key-value pairs and structure elements based on the specified sequence.

    :param output: Path to the output text file
    :param data: Dictionary with keys and their values
    :param structure_filename: List containing text to insert
    :param data_sequence: List specifying the order of dictionary keys and structure elements
    """
    csv_file = "pseudopotentials.csv"
    data = replace_backslash(data_backslash)
    structure = file_to_struct(structure_filename)
    structure_string = {'Definition.of.Atomic.Species': atomic_species(structure, csv_file, data['q']),
                        'Atoms.SpeciesAndCoordinates': atom_spec_coord(structure, csv_file),
                        'Atoms.UnitVectors': atom_unit_vectors(structure)}
    with folder.open(input_file, 'w') as handle:
        for item in structure_string:
            # Write the corresponding structure item
            handle.write(f"{structure_string[item]}\n")
        for item in data:
                # Write the key-value pair from the dictionary
                value = data[item]
                handle.write(f"{item:<35} {value:<35}\n")

# Example dictionary input with keys and values
'''
data = {
    'System.CurrrentDirectory': './',
    'System.Name': 'met',
    'level.of.stdout': 1,
    'level.of.fileout': 1,
    'Species.Number': 2,
    'Atoms.Number': 5,
    'Atoms.SpeciesAndCoordinates.Unit': 'Ang',
    'Atoms.UnitVectors.Unit': 'Ang',
    'scf.XcType': 'GGA-PBE',
    'scf.SpinPolarization': 'off',
    'scf.ElectronicTemperature': 300.0,
    'scf.energycutoff': 120.0,
    'scf.maxIter': 100,
    'scf.EigenvalueSolver': 'cluster',
    'scf.Kgrid': '1 1 1',
    'scf.Mixing.Type': 'rmm-diis',
    'scf.Init.Mixing.Weight': 0.200,
    'scf.Min.Mixing.Weight': 0.001,
    'scf.Max.Mixing.Weight': 0.200,
    'scf.Mixing.History': 7,
    'scf.Mixing.StartPulay': 4,
    'scf.criterion': 1.0e-10,
    'scf.lapack.dste': 'dstevx',
    'MD.Type': 'nomd',
    'MD.maxIter': 1,
    'MD.TimeStep': 1.0,
    'MD.Opt.criterion': 1.0e-4,
    'DATA.PATH': '/storage/praha1/home/parizev/openmx3.9/DFT_DATA19 #default=/storage/praha1/home/parizev/openmx3.9/DFT_DATA19',
    'q': 2,
}
'''

#write_mixed_output(output_file, data, structure_filename)
