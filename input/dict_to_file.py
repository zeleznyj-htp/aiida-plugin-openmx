import re
from structure import atom_spec_coord, atom_unit_vectors
from definition_of_atomic_species import atomic_species
filename = "Crys-MnO.cif"
csv_file = 'pseudopotentials.csv'  # Replace with your actual CSV file path
q = 2  # Example q value

def write_mixed_output(output_file, data, structure, data_sequence):
    """
    Write key-value pairs and structure elements based on the specified sequence.

    :param output_file: Path to the output text file
    :param data: Dictionary with keys and their values
    :param structure: List containing text to insert
    :param data_sequence: List specifying the order of dictionary keys and structure elements
    """
    with open(output_file, 'w') as file:
        for item in data_sequence:
            if item in data:
                # Write the key-value pair from the dictionary
                value = data[item]
                file.write(f"{item:<35} {value:<35}\n")
            elif isinstance(item, str) and item in structure:
                # Write the corresponding structure item
                file.write(f"{structure[item]}\n")
            else:
                # Write any other string directly
                file.write(f"{item}\n")


# Example usage:
output_file = 'output.txt'

#struct = structure("filename.cif")
structure = {'Definition.of.Atomic.Species':atomic_species(filename, csv_file, q),'Atoms.SpeciesAndCoordinates':atom_spec_coord(filename),
             'Atoms.UnitVectors':atom_unit_vectors(filename)}
# Example dictionary input with keys and values
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
}

data_sequence = [
    'System.CurrrentDirectory',
    'System.Name',
    'level.of.stdout',
    'level.of.fileout',
    'Species.Number',
    'Definition.of.Atomic.Species',
    'Atoms.Number',
    'Atoms.SpeciesAndCoordinates.Unit',
    'Atoms.SpeciesAndCoordinates',
    'Atoms.UnitVectors.Unit',
    'Atoms.UnitVectors',
    'scf.XcType',
    'scf.SpinPolarization',
    'scf.ElectronicTemperature',
    'scf.energycutoff',
    'scf.maxIter',
    'scf.EigenvalueSolver',
    'scf.Kgrid',
    'scf.Mixing.Type',
    'scf.Init.Mixing.Weight',
    'scf.Min.Mixing.Weight',
    'scf.Max.Mixing.Weight',
    'scf.Mixing.History',
    'scf.Mixing.StartPulay',
    'scf.criterion',
    'scf.lapack.dste',
    'MD.Type',
    'MD.maxIter',
    'MD.TimeStep',
    'MD.Opt.criterion',
    'DATA.PATH',
]

write_mixed_output(output_file, data, structure, data_sequence)
