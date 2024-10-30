from . import TEST_DIR
from aiida_openmx.parsers import parse_std

def test_parse_Febcc():
    with open(TEST_DIR / 'example_outputs' / 'Febcc2.std') as f:
        lines = f.readlines()
    properties = parse_std(lines)

    assert abs( properties['Utot'] - -179.142122856536) < 1e-6
    assert properties['finished_normally']
    assert properties['converged']

    assert len(properties['convergence']) == 28
    assert None not in properties['convergence']

    assert len(properties['convergence_moment']) == 28
    assert None not in properties['convergence_moment']
