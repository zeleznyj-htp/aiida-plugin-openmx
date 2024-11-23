from numpy import array

import pytest

from aiida.orm import load_code, load_node
from aiida.engine import run_get_pk
from aiida import load_profile

from aiida_openmx.calculations import find_nns, define_ij_pairs

from .common import bccFe, parameters, load_config

load_profile()



@pytest.fixture(scope="module")
def openmx_calc():

    code_openmx, code_jx, options = load_config()
    code = load_code(code_openmx)
    builder = code.get_builder()

    builder.parameters = parameters
    builder.structure = bccFe.as_dict()
    builder.precision = 1
    builder.spin_splits = [2]
    builder.metadata.options = options
    builder.metadata.description = "test openmx"

    _, pk = run_get_pk(builder)

    return pk

def test_jx(openmx_calc):

    code_openmx, code_jx, options = load_config()
    code = load_code(code_jx)

    nns = find_nns(bccFe, nn_prec=1e-3)
    ij_pairs = define_ij_pairs([[0, 0, [1,2,3]]], nns)

    parameters = {
        'Flag.PeriodicSum': 'off',
        'Num.Poles': '60',
        'Num.Kgrid': '1 1 1'
    }

    builder = code.get_builder()
    builder.parameters = parameters
    builder.ij_pairs = ij_pairs
    builder.scfout = load_node(openmx_calc).outputs.retrieved
    builder.metadata.options = options

    _, pk = run_get_pk(builder)
    n = load_node(pk)

    assert len(n.outputs.Jijs['Js']) == 3
