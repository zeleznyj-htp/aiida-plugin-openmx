from numpy import array

import pytest

from aiida.orm import load_code, load_node
from aiida.engine import run_get_pk
from aiida import load_profile

from aiida_openmx.calculations import find_nns, define_ij_pairs

from pymatgen.core import Lattice,Structure
import numpy as np

load_profile()

lattice = Lattice.from_parameters(2.459,2.459,2.459,109.471,109.471,109.471)
bccFe = Structure(lattice,['Fe'],[[0.0,0.0,0.0]],site_properties={'magmom':[np.array([1.0,0.0,0.0])]})

parameters = {
    #"System.CurrrentDirectory": "./",
    'level.of.stdout': 1,
    'level.of.fileout': 1,
    'scf.XcType': 'GGA-PBE',
    'scf.SpinPolarization': 'on',
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
    'HS.fileout': 'on',
}

def load_config():

    try:
        from .test_config import code_openmx, code_jx, options
    except:
        raise Exception("Cannot load test config.")

    return code_openmx, code_jx, options

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
