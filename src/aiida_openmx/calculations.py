"""
Calculations provided by aiida_openmx.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida.orm import SinglefileData, Dict, Int, to_aiida_type, List
from aiida.plugins import DataFactory

from aiida_openmx.input.dict_to_file import write_mixed_output
from aiida_openmx.input.flat import replace_dots
from aiida_openmx.input.structure import get_valence_split

import numbers

class OpenMX(CalcJob):
    """
    AiiDA calculation plugin for the OpenMX DFT code.

    :param structure: structure in the pymatgen Structure format, represented as dict.
    :param parameters: Dict containing the parameters of the calculation for the OpenMX input
                        in the format of keyword: value.
    :param precision: Integer that controls the size of the basis. 1 is smallest, 3 is the largest.
    :param spin_splits: A list of initial spin-polarizations for each atom.

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

        write_mixed_output(input_filename,
                           folder,
                           self.inputs.parameters.get_dict(),
                           self.inputs.structure.get_dict(),
                           self.inputs.precision.value,
                           spin_splits,
                           self.inputs.code.filepath_executable)
        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.retrieve_list = ['input_file','*.std','*.xyz','*.out','*.md','*.Band','*.Dos.*']

        return calcinfo

def dict_dot_serializer(d):
    d_slashes = replace_dots(d)
    return Dict(d_slashes)

def spin_split_validator(spin_split,port):
    if spin_split is not None:
        for s in spin_split:
            if not isinstance(s,numbers.Real):
                return 'spin_splits must be list of real numbers'