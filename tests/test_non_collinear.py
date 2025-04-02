from .common import Mn3NiN, load_config
from aiida.orm import load_code, load_node
from aiida.engine import run_get_pk
from aiida_openmx.input.structure import atom_spec_coord

from aiida import load_profile
load_profile()

def test_atom_spec_coord():
    string = atom_spec_coord(Mn3NiN, [3,3,3,0,0], True, [1,1,1,0,0])
    print(string)
    lines = string.split('\n')
    assert(abs(float(lines[1].split()[7])) - 114.0948 < 1e-3)
    assert(abs(float(lines[2].split()[8])) - 135 < 1e-3)
    assert(abs(float(lines[3].split()[8])) - 116.5651 < 1e-3)
    assert(len(lines) == 7)

def test_run_noncollinear():
    parameters = {
        # "System.CurrrentDirectory": "./",
        'level.of.stdout': 1,
        'level.of.fileout': 1,
        'scf.XcType': 'GGA-PBE',
        'scf.SpinPolarization': 'NC',
        'scf.ElectronicTemperature': 1500.0,
        'scf.energycutoff': 500,
        'scf.maxIter': 2,
        'scf.EigenvalueSolver': 'band',
        'scf.Kgrid': '1 1 1',
        'scf.Constraint.NC.Spin ': 'on',
        'scf.Constraint.NC.Spin.v': 5
    }

    code_openmx, code_jx, options = load_config()
    code = load_code(code_openmx)
    builder = code.get_builder()

    builder.parameters = parameters
    builder.structure = Mn3NiN
    builder.precision = 1
    builder.non_collinear_constraint = [1, 1, 1, 0, 0]
    builder.spin_splits = [3, 3, 3, 0, 0]
    builder.metadata.options = options
    builder.metadata.description = "Mn3NiN test non-collinear"

    _, pk = run_get_pk(builder)

    n = load_node(pk)
    assert n.is_finished_ok

    properties = n.outputs.properties.get_dict()
    assert abs(properties['mulliken']['1']['sum'][0] - 14.1) < 1e-3
    assert abs(properties['mulliken']['3']['diff'][1] - 1.37) < 1e-3
    assert abs(properties['mulliken']['4']['euler'][1][0] - 147.78) < 1e-3




