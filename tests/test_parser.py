from . import TEST_DIR
from aiida_openmx.parsers import parse_std, parse_jx_std

def test_parse_Febcc():
    with open(TEST_DIR / 'example_outputs' / 'Febcc2.std') as f:
        lines = f.readlines()
    properties, version = parse_std(lines)

    assert abs( properties['Utot'] - -179.142122856536) < 1e-6
    assert properties['finished_normally']
    assert properties['converged']

    assert len(properties['convergence']) == 28
    assert None not in properties['convergence']

    assert len(properties['convergence_moment']) == 28
    assert None not in properties['convergence_moment']

    assert version == '3.9.9'

def test_parse_empty():
    lines = ['']
    properties, version = parse_std(lines)
    assert version is None

def test_parse_jx():

    with open(TEST_DIR / 'example_outputs' / 'jx.std') as f:
        lines = f.readlines()
    finished_ok, Jijs = parse_jx_std(lines)

    assert finished_ok
    assert len(Jijs['Js']) == 236
    assert len(Jijs['pairs']) == 236
    assert None not in Jijs['Js']
    assert None not in Jijs['pairs']

def test_parse_jx_empty():
    lines = ['']
    finished_ok, Jijs = parse_jx_std(lines)
    assert not finished_ok

def test_parse_non_collinear():
    with open(TEST_DIR / 'example_outputs' / 'Mn3NiN_non_collinear.std') as f:
        lines = f.readlines()
    properties, version = parse_std(lines)

    assert len(properties['mulliken']) == 5
    assert len(properties['mulliken'][1]['diff']) == 200
    assert abs(properties['mulliken'][1]['diff'][-1] - 3.25) < 1e-5
    assert abs(properties['mulliken'][1]['sum'][-1] - 14.62) < 1e-5
    assert abs(properties['mulliken'][1]['euler'][-1][0] - 114.09) < 1e-5
    assert abs(properties['mulliken'][3]['euler'][-1][1] - 116.57) < 1e-5

    assert None not in properties['convergence_moment']
    assert (properties['convergence_moment'][43] - 0.000298734) / 0.000298734  < 1e-5


