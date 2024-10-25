def flatten_dict(d, parent_key='', sep='.'):
    """
    Recursively flatten a nested dictionary.

    Args:
        d (dict): The nested dictionary to flatten.
        parent_key (str): The base key to use for concatenation (used internally in recursion).
        sep (str): Separator to use between concatenated keys (default is '.').

    Returns:
        dict: A flattened dictionary where nested keys are combined using the separator.
    """
    print(d)
    items = []
    for k, v in d.items():
        # Construct new key by joining parent key and current key
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            # Recursively flatten dictionaries
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            # Add the current key-value pair to the flattened dictionary
            items.append((new_key, v))
    return dict(items)

def unflatten_dict(d, sep='.'):
    """
    Convert a flat dictionary with dot-separated keys into a nested dictionary.

    Args:
        d (dict): The flattened dictionary to convert.
        sep (str): Separator used in the flat dictionary keys (default is '.').

    Returns:
        dict: A nested dictionary reconstructed from the flat dictionary.
    """
    result = {}
    for flat_key, value in d.items():
        keys = flat_key.split(sep)  # Split the flat key by the separator
        d_nested = result
        for key in keys[:-1]:  # Traverse and build the nested structure
            if key not in d_nested:
                d_nested[key] = {}
            d_nested = d_nested[key]
        d_nested[keys[-1]] = value  # Set the final value
    return result

'''
parameters_tmp = {"System.CurrrentDirectory": "./",
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
    'DATA.PATH': '/storage/praha1/home/parizev/openmx3.9/DFT_DATA19', #default=/storage/praha1/home/parizev/openmx3.9/DFT_DATA19'
    'q': 2
}


# Example usage:
input_parameters = {"a.b.c": "./.l", "c.d":"test"}
flattened = flatten_dict(unflatten_dict(parameters_tmp))
print("flattened dict = "+str(flattened))
print(parameters_tmp == flattened)
'''
