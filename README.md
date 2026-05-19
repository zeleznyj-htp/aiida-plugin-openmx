# aiida-openmx

AiiDA plugin for the [OpenMX](http://www.openmx-square.org/) DFT code.

## Features

This plugin provides the following AiiDA components for interacting with OpenMX:

### Calculations

* **`OpenMX`**: A calculation plugin for running standard OpenMX DFT calculations. It supports various features including:
  * Structure input via AiiDA `StructureData` or pymatgen `Structure`.
  * OpenMX parameters provided as a nested dictionary (`Dict`).
  * Collinear and non-collinear spin-polarization (`spin_splits`, `non_collinear_constraint`).
  * Band structure calculations along high-symmetry paths (`bands.critical_points`, `bands.k_path`, `bands.n_band`, `bands.unit_cell`).
  * DFT+U with orbital polarization (`plusU_orbital`, `hubbard_orbital_map`).
  * Restarting calculations from previous states (`retrieve_rst`, `rst_files`).

* **`JxCalculation`**: A calculation plugin for calculating exchange constants (J) from an OpenMX calculation using the `jx` code. 
  * Note: works for collinear calculations only.
  * Requires a spin-polarized calculation and setting `'HS.fileout': 'on'` in the OpenMX calculation.

### Workchains

* **`OpenMXWorkchain`**: A workchain that smoothly runs an OpenMX calculation, then uses its output to run another OpenMX calculation with optionally updated parameters. Useful for performing structure optimizations followed by a static SCF calculation or a band structure calculation. Can also reuse restart files (`rst_files`) between steps.

## Installation

```shell
pip install aiida-openmx
verdi quicksetup  # better to set up a new profile
verdi plugin list aiida.calculations  # should now show openmx and jx
```

## Usage

The `aiida-openmx` plugin allows you to configure calculations using native AiiDA data types while exposing the full power of the OpenMX DFT code. The AiiDA `Dict` parameter maps directly to the OpenMX keyword specifications found in the [OpenMX user manual](http://www.openmx-square.org/openmx_man3.9/index.html).

### 1. Basic SCF Calculation

Here is an example of setting up a basic SCF calculation:

```python
from aiida import load_profile
from aiida.plugins import CalculationFactory
from aiida.orm import Dict, StructureData, Int, List, load_code

load_profile()

OpenMX = CalculationFactory('openmx')
builder = OpenMX.get_builder()

# Set the code
builder.code = load_code('openmx@localhost')

# Set the structure (can also be a standard python dict of a pymatgen structure)
builder.structure = StructureData(pymatgen_structure=...)

# OpenMX parameters are passed as a python dictionary (AiiDA Dict). 
# Keys are exactly as defined in the OpenMX manual (dots in keys are supported).
builder.parameters = Dict(dict={
    'scf.XcType': 'GGA-PBE',
    'scf.energycutoff': 150.0,
    'scf.maxiter': 100,
    'scf.criterion': '1.0e-6',
    'scf.spinpolarization': 'on',
})

# Set precision (integer 1-3). Controls the size of the PAO basis set.
# 1 is the smallest, 3 is the largest.
builder.precision = Int(2)

# For a spin-polarized calculation, you can set the initial spin-splitting 
# for each atom in the structure:
builder.spin_splits = List(list=[2.0, -2.0]) # Assuming a 2-atom structure

# Submit the calculation
from aiida.engine import submit
submit(builder)
```

### 2. Band Structure Calculation

To compute band dispersions, you can define critical points and the $k$-path directly through the `bands` namespace, rather than specifying them manually in the `parameters` dictionary:

```python
# Enable band dispersion in the OpenMX parameters
builder.parameters = Dict(dict={
    # ... other parameters ...
    'Band.dispersion': 'on',
})

# Define the high-symmetry points (critical points)
builder.bands.critical_points = Dict(dict={
    'G': [0.0, 0.0, 0.0],
    'H': [-0.5, 0.5, 0.5],
    'N': [0.0, 0.5, 0.0],
    'P': [0.25, 0.25, 0.25]
})

# Define the paths between the critical points
builder.bands.k_path = List(list=[['G', 'H', 'N', 'G', 'P', 'H', 'N', 'P']])

# Set the number of k-points along each path segment
builder.bands.n_band = Int(15)
```

### 3. Non-Collinear DFT Calculations

OpenMX supports non-collinear DFT calculations, and the plugin provides convenient inputs for setting the constraint flags per atom. 

The magnetic moments have to be provided as site_parameter called "magmom" of the pymatgen Structure. They are specified in the cartesian coordinates. This is used to generate the theta and phi angles for the OpenMX input file. The magnitude of the magnetic moments is specified through the spin_splits.

```python
builder.parameters = Dict(dict={
    # ... other parameters ...
    'scf.SpinPolarization': 'NC',        # Turn on Non-Collinear calculation
    'scf.Constraint.NC.Spin': 'on',      # Apply constraints
    'scf.Constraint.NC.Spin.v': 5.0      # Penalty parameter for the constraint
})

# Specify the initial spin moments for each atom
builder.spin_splits = List(list=[3, 3, 3, 0, 0])

# Specify which atoms should be constrained (1 = constrained, 0 = unconstrained)
builder.non_collinear_constraint = List(list=[1, 1, 1, 0, 0])
```

## Outputs

The AiiDA parsers for this plugin automatically extract useful data from the standard output and make them available as AiiDA output nodes.

### OpenMX Calculation Outputs
- **`output_file`** (`SinglefileData`): The raw standard output file from the OpenMX run.
- **`calculation_info`** (`Dict`): Information about the calculation, including the `openmx_version` and the `aiida_openmx_plugin_version`.
- **`properties`** (`Dict`): Parsed physical properties and convergence metrics, including:
  - `finished_normally` (bool): Whether the calculation reached normal completion.
  - `converged` (bool): Whether the SCF cycle converged.
  - `Utot` (float): The final total energy.
  - `convergence` (list): The history of the SCF convergence criterion (`dUele`).
  - `energy_ele` (list): The electronic energy at each SCF step.
  - `spin_moment` (float): The final total spin moment of the system.
  - `mulliken` (dict): Mulliken population analysis, mapped by atom index (1-based). Each atom contains lists of the `sum` (total charge), `diff` (spin difference), and `euler` angles (theta, phi, for non-collinear calculations) at each SCF step.

### Jx Calculation Outputs
- **`output_file`** (`SinglefileData`): The raw standard output file from the `jx` run.
- **`calculation_info`** (`Dict`): Contains the `aiida_openmx_plugin_version`.
- **`Jijs`** (`Dict`): The computed exchange coupling constants. It contains:
  - `pairs`: A list of tuples defining the interacting pairs `(i, j, Rx, Ry, Rz)`.
  - `Js`: A list of the calculated exchange coupling constants $J$ corresponding to each pair.


## License

MIT

## Contact

vojta.parizek@gmail.com
