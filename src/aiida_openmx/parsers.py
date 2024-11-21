"""
Parsers provided by aiida_openmx.

Register parsers via the "aiida.parsers" entry point in setup.json.
"""

from aiida.common import exceptions
from aiida.engine import ExitCode
from aiida.orm import SinglefileData, Dict
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory

from aiida_openmx import __version__

import logging
import re

from aiida_openmx.calculations import JxCalculation

OpenMXInputFile = CalculationFactory("openmx")
logger = logging.getLogger(__name__)
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
        with self.retrieved.open(output_filename, "rb") as handle1:
            output_node = SinglefileData(file=handle1)
        self.out("output_file", output_node)
        
        with self.retrieved.open(output_filename, "r") as handle2:
            output_lines = handle2.readlines()
        properties, openmx_version = parse_std(output_lines)

        output_dict = Dict(dict=properties)
        self.out("properties", output_dict)

        calculation_info = {
            'openmx_version' : openmx_version,
            'aiida_openmx_plugin_version' : __version__
                            }

        self.out('calculation_info', Dict(calculation_info))
        return ExitCode(0)


class JxParser(Parser):
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
        if not issubclass(node.process_class, JxCalculation):
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
        self.out("output_file", output_node)

        with self.retrieved.open(output_filename, "r") as handle:
            output_lines = handle.readlines()
        finished_ok, Jijs = parse_jx_std(output_lines)

        self.out("Jijs",Dict(Jijs))

        calculation_info = {
            'aiida_openmx_plugin_version': __version__
        }

        self.out('calculation_info', Dict(calculation_info))
        return ExitCode(0)

def parse_std(lines):
    """
    This parses the std of openmx. With MD this contains just the last scf run!
    """

    properties = {}

    properties['finished_normally'] = False
    properties['converged'] = False
    properties['Utot'] = None

    version = None

    for line in lines:
        if 'The calculation was normally finished.' in line:
            properties['finished_normally'] = True
        if 'The SCF was achieved' in line:
            properties['converged'] = True
        if 'Utot' in line:
            match = re.search(r"Utot\s*=\s*(-?\d+\.\d+)", line)
            if match is not None:
                properties['Utot'] = float(match.group(1))
        if 'Ver.' in line:
            version = line.split()[-1]

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
            logger.report("Total SPin Moment found")
            try:
                convergence_moment.append(float(line.split('=')[1]))
            except:
                convergence_moment.append(None)
    properties['convergence_moment'] = convergence_moment

    if len(convergence_moment) > 0:
        properties['spin_moment'] = convergence_moment[-1]
    else:
        properties['spin_moment'] = None


    return properties, version

def parse_jx_std(lines):
    out = {}
    out['pairs'] = []
    out['Js'] = []

    finished_ok = False

    parse = False
    for line in lines:
        if 'Elapsed time' in line:
            finished_ok = True
        if not parse and '----' in line:
            parse = True
            continue
        if parse:
            if len(line) < 2:
                parse = False
                continue
            else:
                try:
                    line_s = line.split()
                    pair = line_s[0:5]
                    pair = tuple([int(i) for i in pair])
                    J = float(line_s[5])
                    out['pairs'].append(pair)
                    out['Js'].append(J)
                except:
                    out['pairs'].append(None)
                    out['Js'].append(None)

    return finished_ok, out

