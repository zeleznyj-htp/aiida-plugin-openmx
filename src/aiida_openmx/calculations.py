"""
Calculations provided by aiida_openmx.

Register calculations via the "aiida.calculations" entry point in setup.json.
"""
from aiida.common import datastructures
from aiida.engine import CalcJob
from aiida.orm import SinglefileData, Dict, List
from aiida.plugins import DataFactory

from aiida_openmx.input.dict_to_file import write_mixed_output

#DiffParameters = DataFactory("openmx")


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
        spec.inputs["metadata"]["options"]["resources"].default = {
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        }
        spec.inputs["metadata"]["options"]["parser_name"].default = "openmx"

        # new ports
        spec.input("metadata.options.output_filename", valid_type=str, default="met.out")
        #spec.input("input_file", valid_type=SinglefileData, help="Input file")
        spec.input("structure_filename", valid_type=SinglefileData, help="Structure file")
        spec.input("csv_file", valid_type=SinglefileData, help="File with PAO basis functions")
        spec.input("parameters", valid_type=Dict, help="Parameters of the calculation")
        spec.input("data_sequence", valid_type=List, help="Sequence of the parameters written to the input file.")
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
        codeinfo.cmdline_params = [input_filename, self.metadata.options.output_filename]
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdout_name = self.metadata.options.output_filename

        # Prepare a `CalcInfo` to be returned to the engine
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        write_mixed_output(input_filename, folder, "parameters", "structure_filename", "data_sequence", "csv_file")
        calcinfo.retrieve_list = [self.metadata.options.output_filename]

        return calcinfo
