"""
Calculations provided by aiida_openmx.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida.orm import SinglefileData, Dict, List
from aiida.plugins import DataFactory

from aiida_openmx.input.dict_to_file import write_mixed_output

class OpenMXInputFile(CalcJob):
    """
    AiiDA calculation plugin wrapping the diff executable.

    Simple AiiDA plugin wrapper for 'diffing' two files.
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
        spec.input("metadata.options.output_filename", valid_type=str, default="met.std")
        spec.input("structure", valid_type=Dict, help="Structure of the material")
        spec.input("parameters", valid_type=Dict, help="Parameters of the calculation")
        spec.output(
            "output_file",
            valid_type=SinglefileData,
            help="output_file",
        )

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
        write_mixed_output(input_filename, folder, self.inputs.parameters.get_dict(), self.inputs.structure.get_dict())
        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]

        #structure_file = self.inputs.structure_file
        #structure_filename = structure_file.filename

        # Store the file in the calculation folder
        #with structure_file.open(mode='rb') as fhandle1:
        #    folder.create_file_from_filelike(fhandle1, structure_filename)

        # Construct the full path to the file in the calculation folder
        #structure_file_path = folder.get_abs_path(structure_filename)
        calcinfo.retrieve_list = [self.metadata.options.output_filename]

        return calcinfo
