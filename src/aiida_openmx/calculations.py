"""
Calculations provided by aiida_openmx.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from aiida.common import datastructures, AttributeDict
from aiida.common.folders import Folder
from aiida.engine import CalcJob, WorkChain, ToContext, if_, ProcessSpec
from aiida.orm import SinglefileData, Dict, Int, Bool, to_aiida_type, List, FolderData, ArrayData, Code
from aiida.plugins import DataFactory, CalculationFactory

from aiida_openmx.input.dict_to_file import write_mixed_output, dict_lowercase
from aiida_openmx.input.flat import replace_dots
from aiida_openmx.input.structure import get_valence_split
from aiida_openmx.input.jx_input import write_jx_input
from aiida_openmx.input.symmetrize_structure import structure_from_output

from aiida_openmx.sd_model.sd_model import Model, input_from_pymatgen
from pymatgen.core import Structure

import numbers
import re


class OpenMX(CalcJob):
    """
    AiiDA calculation plugin for the OpenMX DFT code.

    :param structure: structure in the pymatgen Structure format, represented as dict.
    :param parameters: Dict containing the parameters of the calculation for the OpenMX input
                        in the format of keyword: value.
    :param precision: Integer that controls the size of the basis. 1 is smallest, 3 is the largest.
    :param spin_splits: A list of initial spin-polarizations for each atom.
    :param non_collinear_constraint: Int or [Int], controls whethere the non-collinear constraint is applied
    :param bands.critical_points: Dict defining The points in reciprocal space that define the paths in the k-space.
        Example:
            {
                'G': [0,0,0],
                'H': [-0.5,0.5,0.5],
                'N': [0,0.5,0],
                'P': [0.25,0.25,0.25]
            }
    :param plusU_orbital: Bool or [Bool], specifies whether the yes switch should be used for DFT+U which is used to get
        a state with orbital polarization. If it's a list then this specifies it for each atom.
    :param retrieve_rst: Bool that controls whether the _rst files are retrieved.
    :param rst_files: Folder containing the _rst files used for restarting the calculation.
    :param bands.k_path: A List of list containg the paths along which the bands are calculated using the defintion of
        the critical points.
        Example: [['G','H','N'],['G','P','H']]
    :param bands.n_band: Number of k-points along each line in the path.
    :param bands.unit_cell: 3x3 Array. Definition of the reciprocal unit cell. IF not present then the recirpocal cell
        corresponding to the lattice unit vectors is used.

    When 'scf.spinpolarization' is set to 'NC' then the magnetic moments are read from the structure site property 'magmom'.

    Parameters can contain any openmx keywords with the following exceptions:

    - System.CurrrentDir, System.Name, Atoms.SpeciesAndCoordinates.Unit:
      These should not be specified. These will be overwritten if they are in the input
    - Atoms.UnitVectors.Unit is Ang by default, but this can be specified in the input
    - Species.Number, Atoms.Number: Need not to be specified though they can be.
    - DATA.PATH: this is by default assumed to be ../../DATA_DFT19 with respect to the openmx executable

    The values in parameters can be numbers or strings, lists of numbers or strings, nested strings or numpy arrays.

    If parameter starts with <, then the value will be put into new line and closed with >parameter. So for example:



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
        spec.input("non_collinear_constraint", valid_type=(Int,List), default=lambda: Int(0),
                   serializer=to_aiida_type,
                   help="Can be 1 or 0 or a list of 1s or 0s. Controls which atoms have constrained directions in"
                   "non-collinear calculations.")
        spec.input("plusU_orbital", valid_type=(Bool,List) , default=lambda: Bool(False),
                   help="Controls whether to use the on switch to find orbital polarization with DFT+U."
                        "Can be either boolean or a list of booleans, which then controls this switch for every atom.",
                   serializer=to_aiida_type,
                   )
        spec.input("retrieve_rst", valid_type=Bool, default=lambda: Bool(False),
                   help="Controls whether the _rst files are downloaded, which allow for continuing the calculation.",
                   serializer=to_aiida_type)
        spec.input("rst_files", valid_type=FolderData, required=False, default=None,
                   help="Folder containing the rst files needed to restart OpenMX calculation.")
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

        if isinstance(self.inputs.plusU_orbital,Bool):
            plusU_orbital = self.inputs.plusU_orbital.value
        else:
            plusU_orbital = self.inputs.plusU_orbital.get_list()

        parameters = self.inputs.parameters.get_dict()

        parameters_lower = dict_lowercase(parameters)

        if 'scf\\spinpolarization' in parameters_lower and parameters_lower['scf\\spinpolarization'] == 'nc':
            non_collinear = True
        else:
            non_collinear = False

        if isinstance(self.inputs.non_collinear_constraint,Int):
            non_collinear_constraint = self.inputs.non_collinear_constraint.value
        else:
            non_collinear_constraint = self.inputs.non_collinear_constraint.get_list()

        write_mixed_output(input_filename,
                           folder,
                           parameters,
                           self.inputs.structure.get_dict(),
                           self.inputs.precision.value,
                           spin_splits,
                           non_collinear,
                           non_collinear_constraint,
                           self.inputs.code.filepath_executable,
                           self.inputs.bands.n_band.value,
                           critical_points,
                           k_path,
                           unit_cell,
                           plusU_orbital)
        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = []

        if self.inputs.rst_files is not None:

            #It seems it's not possible to use pattern matching in local_copy_list
            file_names = self.inputs.rst_files.list_object_names()
            rst_directory = None
            for fname in file_names:
                if re.match('.+_rst',fname) is not None:
                    rst_directory = fname
                    break

            if rst_directory is None:
                raise Exception('No rst_files present in the folder.')

            calcinfo.local_copy_list.append((self.inputs.rst_files.uuid, rst_directory, '.'))

            if 'scf.restart' not in parameters:
                parameters['scf.restart'] = 'on'

        calcinfo.retrieve_list = ['input_file','*.std','*.xyz','*.out','*.md','*.Band','*.Dos.*','*.scfout',
                                  self.metadata.options.output_filename]
        if self.inputs.retrieve_rst:
            calcinfo.retrieve_list.append(('*_rst/*','.',2))

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
        if isinstance(pair[2],int):
            order = [pair[2]]
        else:
            order = pair[2]
        for nn in nns[i]:
            if nn[0] == j and nn[3] in order:
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

def safe_assign(builder, name, value, serializer=None):
    if serializer and isinstance(value, dict):
        setattr(builder, name, serializer(value))
    else:
        setattr(builder, name, value)

class OpenMXWorkchain(WorkChain):
    """
    A Workchain that runs an OpenMX calculation, then uses its output to run another
    OpenMX calculation with optionally updated parameters.
    """

    @classmethod
    def define(cls, spec: ProcessSpec):
        super().define(spec)

        # Inputs for the first calculation
        # Metadata
        spec.input('code', valid_type=Code)
        spec.input('metadata1', valid_type=Dict, serializer=dict_dot_serializer)
        #spec.input("metadata.options.output_filename", valid_type=str, default="openmx.std")
        spec.input("structure", valid_type=Dict, help="Structure of the material",
                   serializer=to_aiida_type)
        spec.input("parameters", valid_type=Dict, help="Parameters of the calculation",
                   serializer=dict_dot_serializer)
        spec.input("precision", valid_type=Int, default=lambda: Int(2),
                   help="Controls the size of the basis, 1 is the smallest and 3 the largest.",
                   serializer=to_aiida_type)
        spec.input("spin_splits", valid_type=List, default=None,
                   help="The initial spin-polarization for each atom.",
                   serializer=to_aiida_type,
                   validator=spin_split_validator,
                   required=False)
        spec.input("non_collinear_constraint", valid_type=(Int, List), default=lambda: Int(0),
                   serializer=to_aiida_type,
                   help="Can be 1 or 0 or a list of 1s or 0s. Controls which atoms have constrained directions in"
                        "non-collinear calculations.")
        spec.input("plusU_orbital", valid_type=(Bool, List), default=lambda: Bool(False),
                   help="Controls whether to use the on switch to find orbital polarization with DFT+U."
                        "Can be either boolean or a list of booleans, which then controls this switch for every atom.",
                   serializer=to_aiida_type,
                   )
        spec.input("retrieve_rst", valid_type=Bool, default=lambda: Bool(False),
                   help="Controls whether the _rst files are downloaded, which allow for continuing the calculation.",
                   serializer=to_aiida_type)
        spec.input("rst_files", valid_type=FolderData, required=False, default=None,
                   help="Folder containing the rst files needed to restart OpenMX calculation.")

        # Bands namespace
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

        # Control logic
        spec.input("modify_for_second", valid_type=Bool, default=lambda: Bool(False), serializer=to_aiida_type)
        spec.input("override_inputs", valid_type=Dict, serializer=dict_dot_serializer, help="Input overrides for the second calculation")
        spec.output("output_file1", valid_type=SinglefileData, help="output_file")
        spec.output("properties1", valid_type=Dict, help="Output properties of the calculation")
        spec.output("calculation_info1", valid_type=Dict, help="Shows versions of the software used to run the calculation.")
        spec.output("output_file2", valid_type=SinglefileData, help="output_file")
        spec.output("properties2", valid_type=Dict, help="Output properties of the calculation")
        spec.output("calculation_info2", valid_type=Dict, help="Shows versions of the software used to run the calculation.")

        # Outline
        spec.outline(
            cls.run_first_calc,
            cls.inspect_first,
            if_(cls.should_run_second)(
                cls.run_second_calc,
                cls.inspect_second
            ),
            cls.finalize
        )

        # Exit codes
        spec.exit_code(100, 'ERROR_FIRST_CALC_FAILED', message='The first OpenMX calculation failed.')
        spec.exit_code(101, 'ERROR_SECOND_CALC_FAILED', message='The second OpenMX calculation failed.')

    def build_openmx_inputs(self, override_inputs=None):
        inputs = AttributeDict({
            'structure': self.inputs.structure,
            'parameters': self.inputs.parameters,
            'precision': self.inputs.precision,
            'metadata': dict(self.inputs.metadata1),
            'code': self.inputs.code,
            'retrieve_rst': self.inputs.retrieve_rst,
        })

        optional_inputs = ['spin_splits', 'non_collinear_constraint', 'plusU_orbital', 'rst_files']
        for key in optional_inputs:
            if key in self.inputs:
                inputs[key] = self.inputs[key]

        inputs['bands'] = AttributeDict()
        for key in ['critical_points', 'k_path', 'n_band', 'unit_cell']:
            if key in self.inputs.bands:
                inputs['bands'][key] = self.inputs.bands[key]

        # Apply clean overrides
        if override_inputs:
            data = override_inputs.get_dict()

            # Override parameters
            if 'parameters' in data:
                param_dict = inputs['parameters'].get_dict()
                param_dict.update(data['parameters'])
                inputs['parameters'] = dict_dot_serializer(param_dict)

            # Override bands
            if 'bands' in data:
                for k, v in data['bands'].items():
                    inputs['bands'][k] = v

            # Override any top-level simple keys
            for key in ['precision', 'retrieve_rst', 'spin_splits', 'non_collinear_constraint', 'plusU_orbital']:
                if key in data:
                    inputs[key] = data[key]

        return inputs

    def run_first_calc(self):
        self.report('Running first OpenMX calculation...')
        inputs = self.build_openmx_inputs()
        future = self.submit(OpenMX, **inputs)
        return ToContext(first_calc=future)

    def inspect_first(self):
        if not self.ctx.first_calc.is_finished_ok:
            self.report('First OpenMX calculation failed.')
            return self.exit_codes.ERROR_FIRST_CALC_FAILED
        self.report('First OpenMX calculation finished successfully.')

    def should_run_second(self):
        return self.inputs.modify_for_second

    def run_second_calc(self):
        self.report('Preparing second OpenMX calculation with modified parameters...')

        # Read the aiida.out content
        with self.ctx.first_calc.outputs.retrieved.open('aiida.out', "r") as handle:
            output_text = handle.read()
        # Extract the original structure (as dict)
        original_structure_dict = self.inputs.structure.get_dict()

        # Generate new structure using your existing function
        structure = structure_from_output(output_text, original_structure_dict)  # returns a StructureData node

        override_inputs = self.inputs.override_inputs if 'override_inputs' in self.inputs else {}
        inputs = self.build_openmx_inputs(override_inputs=override_inputs)
        inputs['structure'] = structure  # the AiiDA StructureData node
        future = self.submit(OpenMX, **inputs)
        return ToContext(second_calc=future)

    def inspect_second(self):
        if not self.ctx.second_calc.is_finished_ok:
            self.report('Second OpenMX calculation failed.')
            return self.exit_codes.ERROR_SECOND_CALC_FAILED
        self.report('Second OpenMX calculation finished successfully.')

    def finalize(self):
        self.out('properties1', self.ctx.first_calc.outputs.properties)
        self.out('output_file1', self.ctx.first_calc.outputs.output_file)
        self.out('calculation_info1', self.ctx.first_calc.outputs.calculation_info)
        if 'second_calc' in self.ctx and 'properties' in self.ctx.second_calc.outputs:
            self.out('properties2', self.ctx.second_calc.outputs.properties)
            self.out('output_file2', self.ctx.second_calc.outputs.output_file)
            self.out('calculation_info2', self.ctx.second_calc.outputs.calculation_info)
