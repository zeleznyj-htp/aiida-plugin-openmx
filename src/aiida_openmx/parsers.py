"""
Parsers provided by aiida_openmx.

Register parsers via the "aiida.parsers" entry point in setup.json.
"""

from aiida.common import exceptions
from aiida.engine import ExitCode
from aiida.orm import SinglefileData
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory

import re

OpenMXInputFile = CalculationFactory("openmx")


class OpenMXParser(Parser):
    """
    Parser class for parsing output of calculation.
    """

    def __init__(self, node):
        """
        Initialize Parser instance

        Checks that the ProcessNode being passed was produced by a DiffCalculation.

        :param node: ProcessNode of calculation
        :param type node: :class:`aiida.orm.nodes.process.process.ProcessNode`
        """
        super().__init__(node)
        if not issubclass(node.process_class, OpenMXInputFile):
            raise exceptions.ParsingError("Can only parse OpenMXCalculation")

    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.

        :returns: an exit code, if parsing fails (or nothing if parsing succeeds)
        """
        output_filename = self.node.get_option("output_filename")

        # Check that folder content is as expected
        files_retrieved = self.retrieved.list_object_names()
        files_expected = [output_filename]
        # Note: set(A) <= set(B) checks whether A is a subset of B
        if not set(files_expected) <= set(files_retrieved):
            self.logger.error(f"Found files '{files_retrieved}', expected to find '{files_expected}'")
            return self.exit_codes.ERROR_MISSING_OUTPUT_FILES

        # add output file
        self.logger.info(f"Parsing '{output_filename}'")
        with self.retrieved.open(output_filename, "rb") as handle:
            output_node = SinglefileData(file=handle)
            output_lines = handle.readlines()
        properties = parse_std(output_lines)
        self.out("openmx", output_node)

        return ExitCode(0)

def parse_std(lines):
    """
    This parses the std of openmx. With MD this contains just the last scf run!
    """

    properties = {}

    properties['finished_normally'] = False
    properties['converged'] = False
    properties['Utot'] = None

    for line in lines:
        if 'The calculation was normally finished.' in line:
            properties['finished_normally'] = True
        if 'The SCF was achieved' in line:
            properties['converged'] = True
        if 'Utot' in line:
            match = re.search(r"Utot\s*=\s*(-?\d+\.\d+)", line)
            if match is not None:
                properties['Utot'] = float(match.group(1))

    convergence = []
    for line in lines:
        if 'dUele' in line:
            try:
                convergence.append(float(line.split('=')[2]))
            except:
                convergence.append(None)
    properties['convergence'] = convergence

    convergence_moment = []
    for line in lines:
        if 'Total Spin Moment' in line:
            try:
                convergence_moment.append(float(line.split('=')[1]))
            except:
                convergence_moment.append(None)
    properties['convergence_moment'] = convergence_moment

    properties['spin_moment'] = convergence_moment[-1]

    return properties
