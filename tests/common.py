from pymatgen.core import Lattice,Structure
import numpy as np

lattice = Lattice.from_parameters(2.459,2.459,2.459,109.471,109.471,109.471)
bccFe = Structure(lattice,['Fe'],[[0.0,0.0,0.0]],site_properties={'magmom':[np.array([1.0,0.0,0.0])]})

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