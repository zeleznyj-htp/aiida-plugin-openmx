from .common import bccFe, parameters, load_config
from aiida.orm import load_code, load_node
from aiida.engine import run_get_pk

from aiida import load_profile
load_profile()

def test_bands():
    parameters['Band.dispersion'] = 'on'

    code_openmx, code_jx, options = load_config()
    code = load_code(code_openmx)
    builder = code.get_builder()

    builder.parameters = parameters
    builder.structure = bccFe.as_dict()
    builder.precision = 1
    builder.spin_splits = [2]
    builder.metadata.options = options
    builder.metadata.description = "test openmx"
    builder.bands.critical_points = {
        'G': [0, 0, 0],
        'H': [-0.5, 0.5, 0.5],
        'N': [0, 0.5, 0],
        'P': [0.25, 0.25, 0.25]
    }
    builder.bands.k_path = [['G', 'H', 'N', 'G', 'P', 'H', 'N', 'P']]

    _, pk = run_get_pk(builder)

    n = load_node(pk)
    assert n.is_finished_ok
    assert 'aiida.Band' in n.outputs.retrieved.list_object_names()



