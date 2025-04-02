from pymatgen.core import Lattice,Structure
import numpy as np

lattice = Lattice.from_parameters(2.459,2.459,2.459,109.471,109.471,109.471)
bccFe = Structure(lattice,['Fe'],[[0.0,0.0,0.0]],site_properties={'magmom':[np.array([1.0,0.0,0.0])]})
Mn3NiN = {'@module': 'pymatgen.core.structure',
 '@class': 'Structure',
 'charge': 0,
 'lattice': {'matrix': [[3.838489, 0.0, 0.0],
   [0.0, 3.838489, 0.0],
   [0.0, 0.0, 3.838489]],
  'pbc': (True, True, True),
  'a': 3.838489,
  'b': 3.838489,
  'c': 3.838489,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 90.0,
  'volume': 56.55628849330412},
 'properties': {},
 'sites': [{'species': [{'element': 'Mn', 'occu': 1}],
   'abc': [0.0, 0.5, 0.5],
   'properties': {'magmom': [2, -1, -1], 'magmom_basis': 'cart'},
   'label': 'Mn',
   'xyz': [0.0, 1.9192445, 1.9192445]},
  {'species': [{'element': 'Mn', 'occu': 1}],
   'abc': [0.5, 0.5, 0.0],
   'properties': {'magmom': [-1, -1, 2], 'magmom_basis': 'cart'},
   'label': 'Mn',
   'xyz': [1.9192445, 1.9192445, 0.0]},
  {'species': [{'element': 'Mn', 'occu': 1}],
   'abc': [0.5, 0.0, 0.5],
   'properties': {'magmom': [-1, 2, -1], 'magmom_basis': 'cart'},
   'label': 'Mn',
   'xyz': [1.9192445, 0.0, 1.9192445]},
  {'species': [{'element': 'Ni', 'occu': 1}],
   'abc': [0.0, 0.0, 0.0],
   'properties': {'magmom': [0, 0, 1], 'magmom_basis': 'cart'},
   'label': 'Ni',
   'xyz': [0.0, 0.0, 0.0]},
  {'species': [{'element': 'N', 'occu': 1}],
   'abc': [0.5, 0.5, 0.5],
   'properties': {'magmom': [0, 0, 1], 'magmom_basis': 'cart'},
   'label': 'N',
   'xyz': [1.9192445, 1.9192445, 1.9192445]}]}

parameters_init = {
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