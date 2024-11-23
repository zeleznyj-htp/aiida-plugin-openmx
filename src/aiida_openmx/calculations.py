"""
Calculations provided by aiida_openmx.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from aiida.common import datastructures
from aiida.common.folders import Folder
from aiida.engine import CalcJob
from aiida.orm import SinglefileData, Dict, Int, to_aiida_type, List, FolderData, ArrayData
from aiida.plugins import DataFactory

from aiida_openmx.input.dict_to_file import write_mixed_output
from aiida_openmx.input.flat import replace_dots
from aiida_openmx.input.structure import get_valence_split
from aiida_openmx.input.jx_input import write_jx_input

from aiida_openmx.sd_model.sd_model import Model, input_from_pymatgen
from pymatgen.core import Structure

import numbers

class OpenMX(CalcJob):
    """
    AiiDA calculation plugin for the OpenMX DFT code.

    :param structure: structure in the pymatgen Structure format, represented as dict.
    :param parameters: Dict containing the parameters of the calculation for the OpenMX input
                        in the format of keyword: value.
    :param precision: Integer that controls the size of the basis. 1 is smallest, 3 is the largest.
    :param spin_splits: A list of initial spin-polarizations for each atom.
    :param bands.critical_points: Dict defining The points in reciprocal space that define the paths in the k-space.
        Example:
            {
                'G': [0,0,0],
                'H': [-0.5,0.5,0.5],
                'N': [0,0.5,0],
                'P': [0.25,0.25,0.25]
            }
    :param bands.k_path: A List of list containg the paths along which the bands are calculated using the defintion of
        the critical points.
        Example: [['G','H','N'],['G','P','H']]
    :param bands.n_band: Number of k-points along each line in the path.
    :param bands.unit_cell: 3x3 Array. Definition of the reciprocal unit cell. IF not present then the recirpocal cell
        corresponding to the lattice unit vectors is used.

    Parameters can contain any openmx keywords with the following exceptions:

    - System.CurrrentDir, System.Name, Atoms.SpeciesAndCoordinates.Unit:
      These should not be specified. These will be overwritten if they are in the input
    - Atoms.UnitVectors.Unit is Ang by default, but this can be specified in the input
    - Species.Number, Atoms.Number: Need not to be specified though they can be.
    - DATA.PATH: this is by default assumed to be ../../DATA_DFT19 with respect to the openmx executable

    """

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        super().define(spec)

        # set default values for AiiDA options
        spec.input('metadata.options.resources', valid_type=dict,
                   default={'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        spec.input('metadata.options.parser_name', valid_type=str, default='openmx')
        # new ports
        spec.input("metadata.options.output_filename", valid_type=str, default="openmx.std")
        spec.input("structure", valid_type=Dict, help="Structure of the material",
                   serializer=to_aiida_type)
        spec.input("parameters", valid_type=Dict, help="Parameters of the calculation",
                   serializer=dict_dot_serializer)
        spec.input("precision", valid_type=Int, default=lambda: Int(2),
                   help="Controls the size of the basis, 1 is the smallest and 3 the largest.",
                   serializer=to_aiida_type)
        spec.input("spin_splits", valid_type=List , default=None,
                    help="The initial spin-polarization for each atom.",
                    serializer=to_aiida_type,
                    validator=spin_split_validator,
                    required=False)
        spec.input_namespace('bands')
        spec.input("bands.critical_points", valid_type=Dict, default=None,
                   help="The coordinates of th critical points",
                   serializer=to_aiida_type,
                   required=False)
        spec.input("bands.k_path", valid_type=List, default=None,
                   help="The path defined by the critical points",
                   serializer=to_aiida_type,
                   required=False)
        spec.input("bands.n_band", valid_type=Int, default=lambda: Int(15),
                   help="Number of plotted points between two critical points in the band plot",
                   serializer=to_aiida_type)
        spec.input("bands.unit_cell", valid_type=ArrayData, default=None,
                   help="The reciprocal unit cell with respect to which the critical points are defined."
                        "If not present, openmx will use by default the unit cell corresponding to the lattice unit cell.",
                   serializer=to_aiida_type,
                   required=False)
        spec.output("output_file", valid_type=SinglefileData, help="output_file")
        spec.output("properties", valid_type=Dict, help="Output properties of the calculation")
        spec.output("calculation_info", valid_type=Dict, help="Shows versions of the software used to run the calculation.")

        spec.exit_code(
            300,
            "ERROR_MISSING_OUTPUT_FILES",
            message="Calculation did not produce all expected output files.",
        )

    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files
            needed by the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        input_filename = "input_file"
        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = [input_filename]
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename

        if self.inputs.spin_splits is None:
            spin_splits = None
        else:
            spin_splits = self.inputs.spin_splits.get_list()

        if self.inputs.bands.critical_points is None:
            critical_points = None
        else:
            critical_points = self.inputs.bands.critical_points.get_dict()

        if self.inputs.bands.k_path is None:
            k_path = None
        else:
            k_path = self.inputs.bands.k_path.get_list()

        if self.inputs.bands.unit_cell is None:
            unit_cell = None
        else:
            unit_cell = self.inputs.bands.unit_cell.get_array()

        write_mixed_output(input_filename,
                           folder,
                           self.inputs.parameters.get_dict(),
                           self.inputs.structure.get_dict(),
                           self.inputs.precision.value,
                           spin_splits,
                           self.inputs.code.filepath_executable,
                           self.inputs.bands.n_band.value,
                           critical_points,
                           k_path,
                           unit_cell)
        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.retrieve_list = ['input_file','*.std','*.xyz','*.out','*.md','*.Band','*.Dos.*','*.scfout',
                                  self.metadata.options.output_filename]

        return calcinfo

class JxCalculation(CalcJob):
    """
    Calculates the exchange constants from an OpenMX calculations using the jx code.

    -Works for collinear calcualtions only!
    - Requires a spin-polarized calculation
    - Requires setting 'HS.fileout' to 'on' in the OpenMX calculation.

    :param scfout: The .scfout file from the OpenMX  calculation. It can be either a SinglefileData or FolderData,
        in which case the folder must contain a single scfout file.
    :param parameters: The parameters for the jx calculation.
    :param ij_pairs: The ij pairs to calculate in the format [(i,j,R_x,R_y,R_z),...], where i,j are the atom labels
        and Rx,Ry,Rz is the unit cell translation vector of the j atom. Here i,j start from 1!!!
    """

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        super().define(spec)

        # set default values for AiiDA options
        spec.input('metadata.options.resources', valid_type=dict,
                   default={'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        spec.input('metadata.options.parser_name', valid_type=str, default='jx')
        # new ports
        spec.input("metadata.options.output_filename", valid_type=str, default="openmx.std")
        spec.input("scfout", valid_type=(SinglefileData,FolderData), help="scfout from OpenMX calculation",
                   validator=scfout_validator)
        spec.input("parameters", valid_type=Dict, help="Parameters of the calculation",
                   serializer=dict_dot_serializer)
        spec.input("ij_pairs", valid_type=List, help="The ij pairs to calculate",
                   serializer=to_aiida_type)
        spec.output("output_file", valid_type=SinglefileData, help="output_file")
        spec.output("Jijs", valid_type=Dict, help="Output properties of the calculation")
        spec.output("calculation_info", valid_type=Dict, help="Shows versions of the software used to run the calculation.")

        spec.exit_code(
            300,
            "ERROR_MISSING_OUTPUT_FILES",
            message="Calculation did not produce all expected output files.",
        )

    def prepare_for_submission(self, folder):
        """
        Create input files.

        :param folder: an `aiida.common.folders.Folder` where the plugin should temporarily place all files
            needed by the calculation.
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        input_filename = "jx.config"
        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = ['aiida.scfout', input_filename]
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.retrieve_list = [input_filename,self.metadata.options.output_filename]

        if isinstance(self.inputs.scfout, SinglefileData):
            calcinfo.local_copy_list = [(self.inputs.scfout.uuid, self.inputs.scfout.filename, 'aiida.scfout')]
        else:
            file_names = self.inputs.scfout.list_object_names()
            scfout_filename = [x for x in file_names if 'scfout' in x][0]
            calcinfo.local_copy_list = [(self.inputs.scfout.uuid, scfout_filename, 'aiida.scfout')]

        write_jx_input(input_filename,
                       folder,
                       self.inputs.parameters.get_dict(),
                       self.inputs.ij_pairs,
                       )


        return calcinfo

def find_nns(structure,ncells=2,nns=10,nn_prec=None):
    """
    Returns the nearest neighbors for a given structure.
    """
    if isinstance(structure,dict):
        structure = Structure.from_dict(structure)
    config = input_from_pymatgen(structure,nns=nns)
    m = Model(config=config)
    m.nn_n_cells = ncells
    if nn_prec is not None:
        m.nn_prec = nn_prec

    nns = m.find_nn_new()
    return nns

def define_ij_pairs(pairs,nns):
    """
    Creates a definition of ij pairs that can be used as an input in JxCalculation.

    :param pairs: Defines which neighbors to consider. The format is [[i,j,order],...], where i,j are the atom
        labels and order is the order of the hopping. Here i,j start from zero!!!
    :param nns: nearest neighbors as returns by find_nns.
    """

    out = []

    for pair in pairs:

        i = pair[0]
        j = pair[1]
        ij_done = []
        for nn in nns[i]:
            if nn[0] == j and nn[3] in pair[2]:
                if nn[3] not in ij_done:
                    out.append((i+1,j+1,*nn[1]))
                    ij_done.append(nn[3])

    return out

def dict_dot_serializer(d):
    d_slashes = replace_dots(d)
    return Dict(d_slashes)

def spin_split_validator(spin_split,port):
    if spin_split is not None:
        for s in spin_split:
            if not isinstance(s,numbers.Real):
                return 'spin_splits must be list of real numbers'

def scfout_validator(scfout,port):
    if isinstance(scfout,FolderData):
        file_names = scfout.list_object_names()
        scfouts = [x for x in file_names if 'scfout' in x]
        if len(scfouts) != 1:
            return 'the scfout folder must contain exactly one scfout file'
