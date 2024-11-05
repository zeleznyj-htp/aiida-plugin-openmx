from aiida.engine import submit
from aiida_openmx.input.flat import replace_dots

SinglefileData = DataFactory('core.singlefile')
Dict = DataFactory('core.dict')
#List = DataFactory('core.list')
code = load_code('openmx@tarkil')
builder = code.get_builder()
#builder.structure_file = SinglefileData(file='/home/parizek/aiida-openmx/Methane.json')
#builder.csv_file = SinglefileData(file='/home/parizek/aiida-openmx/pseudopotentials.csv')
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
structure = {'@module': 'pymatgen.core.structure',
'@class': 'Structure', 'charge': 0.0,
'lattice': {'matrix': [[10.0, 0.0, 6.123233995736766e-16], [1.6081226496766364e-15, 10.0, 6.123233995736766e-16], [0.0, 0.0, 10.0]], 'pbc': (True, True, True), 'a': 10.0, 'b': 10.0, 'c': 10.0, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0, 'volume': 1000.0},
'properties': {},
'sites': [{'species': [{'element': 'H', 'occu': 1.0}], 'abc': [0.911, 0.93707, 0.0], 'properties': {}, 'label': 'A2', 'xyz': [9.110000000000001, 9.3707, 1.1316165050501245e-15]}, {'species': [{'element': 'H', 'occu': 1.0}], 'abc': [0.0, 0.06293, 0.911], 'properties': {}, 'label': 'A3', 'xyz': [1.0119915834415073e-16, 0.6293, 9.11]}, {'species': [{'element': 'H', 'occu': 1.0}], 'abc': [0.0, 0.06293, 0.089], 'properties': {}, 'label': 'A4', 'xyz': [1.0119915834415073e-16, 0.6293, 0.8899999999999999]}, {'species': [{'element': 'H', 'occu': 1.0}], 'abc': [0.089, 0.93707, 0.0], 'properties': {}, 'label': 'A5', 'xyz': [0.8900000000000015, 9.3707, 6.282866706005624e-16]}, {'species': [{'element': 'C', 'occu': 1.0}], 'abc': [0.0, 0.0, 0.0], 'properties': {}, 'label': 'A1', 'xyz': [0.0, 0.0, 0.0]}]
}
builder.parameters = Dict(dict=replace_dots(parameters_tmp))
builder.structure = Dict(dict=structure)
builder.metadata.options.withmpi = False
builder.metadata.options.resources = {
    'num_machines': 1,
    'num_mpiprocs_per_machine': 1,
}
builder.metadata.description = "test openmx"
builder.metadata.options.output_filename = "met.std"

print(submit(builder))
