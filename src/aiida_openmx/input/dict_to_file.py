from aiida_openmx.input.structure import atom_spec_coord, atom_unit_vectors
from aiida_openmx.input.definition_of_atomic_species import atomic_species, get_elements, band_path
from aiida_openmx.input.flat import replace_backslash
from aiida.orm import Dict
from pymatgen.core import Structure

def write_mixed_output(input_file, folder, data_backslash, structure, precision, spin_split=None, executable_path=None, n_bands = None, critical_points = None, k_path = None):
    """
    Write key-value pairs and structure elements based on the specified sequence.
    :param output: Path to the output text file
    :param data: Dictionary with keys and their values
    :param structure_filename: List containing text to insert
    """
    data = replace_backslash(data_backslash)
    structure_string = {'Definition.of.Atomic.Species': atomic_species(structure, precision),
                        'Atoms.SpeciesAndCoordinates': atom_spec_coord(structure, spin_split),
                        'Atoms.UnitVectors': atom_unit_vectors(structure)}
    with folder.open(input_file, 'w') as handle:
        for item in structure_string:
            # Write the corresponding structure item
            handle.write(f"{structure_string[item]}\n")
        handle.write('Atoms.SpeciesAndCoordinates.Unit   FRAC\n')
        handle.write('Atoms.UnitVectors.Unit  Ang\n')

        #these should not be changed by user since they are internal to how aiida does the calculation
        data['System.Name'] = 'aiida'
        data['System.CurrrentDir'] = './'

        #in case the 'Atoms.SpeciesAndCoordinates.Unit' is present in the dict we remove it
        data.pop('Atoms.SpeciesAndCoordinates.Unit',None)
        data.pop('Atoms.UnitVectors.Unit  Ang\n',None)

        #the DATA.PATH should be fine like this in most cases
        if executable_path is not None and 'DATA.PATH' not in data:
            data['DATA.PATH'] = str(executable_path.parents[1] / 'DFT_DATA19')

        if 'Species.Number' not in data:
            data['Species.Number'] = len(get_elements(structure))

        #these variables are not necessary to be set by user as they can be determined from the structure
        if 'Atoms.Number' not in data:
            data['Atoms.Number'] = len(structure['sites'])

        for item in sorted(data.keys()):
                # Write the key-value pair from the dictionary
                value = data[item]
                handle.write(f"{item:<35} {value:<35}\n")
        if ('Band.dispersion' in data) and (data['Band.dispersion'] in ['on','On','ON']):
            handle.write(band_path(n_bands, critical_points, k_path)+'\n')