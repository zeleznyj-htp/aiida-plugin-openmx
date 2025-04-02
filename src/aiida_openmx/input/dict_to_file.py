from aiida.common import InputValidationError
from aiida_openmx.input.structure import atom_spec_coord, atom_unit_vectors
from aiida_openmx.input.definition_of_atomic_species import atomic_species, get_elements, band_path, band_kpath_unit_cell
from aiida_openmx.input.flat import replace_backslash

import numpy as np

def dict_lowercase(d):
    d2 = {}
    for key in d:
        if isinstance(d[key],str):
            val = d[key].lower()
        else:
            val = d[key]
        d2[key.lower()] = val

    return d2

def convert_to_string(value):
    """
    Converts the input value to a string.
    - For numbers or strings, returns the string representation.
    - For 1D lists or numpy arrays, returns a space-separated string of elements.
    - For 2D lists or numpy arrays, returns a newline-separated string of space-separated rows.

    Args:
        value (int, float, str, list, numpy.ndarray): The input value to convert.

    Returns:
        str: The formatted string.
    """
    if isinstance(value, (int, float, str)):  # Handle numbers and strings
        return str(value)

    elif ( isinstance(value, list) and all(isinstance(item, (int, float, str)) for item in value) ) \
        or ( isinstance(value, np.ndarray) and value.ndim == 1 ):
        # Handle 1D lists and numpy arrays
        return ' '.join(str(item) for item in value)

    elif (isinstance(value, list) and all(isinstance(item, list) for item in value)) or \
         (isinstance(value, np.ndarray) and value.ndim == 2):
        # Handle 2D lists and numpy arrays
        return '\n'.join(' '.join(str(sub_item) for sub_item in row) for row in value)

    else:
        raise TypeError("Input must be an int, float, str, list, or numpy array.")


def write_mixed_output(input_file, folder, data_backslash, structure, precision, spin_split=None,
                       non_collinear = False, non_collinear_constraint = 1, executable_path=None,
                       n_bands = None, critical_points = None, k_path = None, unit_cell=None, plusU_orbital=False):
    """
    Write key-value pairs and structure elements based on the specified sequence.
    """
    data = replace_backslash(data_backslash)

    if 'scf.Hubbard.U' not in data or data['scf.Hubbard.U'] == 'off':
        plusU_orbital = None

    structure_string = {'Definition.of.Atomic.Species': atomic_species(structure, precision),
                        'Atoms.SpeciesAndCoordinates': atom_spec_coord(structure, spin_split, non_collinear, non_collinear_constraint, plusU_orbital),
                        'Atoms.UnitVectors': atom_unit_vectors(structure)}
    with folder.open(input_file, 'w') as handle:
        for item in structure_string:
            # Write the corresponding structure item
            handle.write(f"{structure_string[item]}\n")
        handle.write('Atoms.SpeciesAndCoordinates.Unit   FRAC\n')
        handle.write('Atoms.UnitVectors.Unit  Ang\n')

        #these should not be changed by user since they are internal to how aiida does the calculation
        data['System.Name'] = 'aiida'
        data['System.CurrentDir'] = './'

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
                if not item.startswith('<'):
                    handle.write(f"{item:<35} {convert_to_string(value):<35}\n")
                else:
                    handle.write(f"{item}\n")
                    handle.write(convert_to_string(value)+"\n")
                    handle.write(f"{item[1:]}>\n")

        if ('Band.dispersion' in data) and (data['Band.dispersion'] in ['on','On','ON']):
            if None in [n_bands, critical_points, k_path]:
                raise InputValidationError('For bands calcualation n_bands, critical_points and k_path must be defined.')
            if unit_cell is not None:
                handle.write(band_kpath_unit_cell(unit_cell))
            handle.write(band_path(n_bands, critical_points, k_path))