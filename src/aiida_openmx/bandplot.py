import re
import numpy as np
import matplotlib.pyplot as plt
from aiida_shell import launch_shell_job
from aiida.plugins import CalculationFactory
from aiida.orm import SinglefileData

def run_bandgnu(bandgnu_path, retrieved):
    results, node = launch_shell_job(
        bandgnu_path,
        arguments='aiida.Band',
        nodes={
            'directory': retrieved,
        },
        outputs=['aiida.BANDDAT*', 'aiida.GNUBAND']
    )
    return results, node

# Function to parse the Gnuplot script
def parse_gnuplot_script(file_lines):
    config = {}
    for line in file_lines:
        # Extract x and y range
        if line.startswith("set xra"):
            config['xrange'] = list(map(float, re.findall(r"[-+]?\d*\.?\d+|\d+", line)))
        elif line.startswith("set yra"):
            config['yrange'] = list(map(float, re.findall(r"[-+]?\d*\.?\d+|\d+", line)))
        # Extract y-label
        elif line.startswith("set ylabel"):
            config['ylabel'] = re.search(r'"(.*?)"', line).group(1)
        # Extract xtics
        elif line.startswith("set xtics"):
            matches = re.findall(r'"(.*?)"\s+([-+]?\d*\.?\d+|\d+)', line)
            config['xticks'] = [(label, float(pos)) for label, pos in matches]
        # Extract data files
        elif line.startswith("plot"):
            config['data_files'] = re.findall(r'"(.*?)"', line)
    return config

def read_band_data(file_lines):
    data = []
    block = []
    for line in file_lines:
        if line.strip():  # Non-empty line
            x, y = map(float, line.split())
            block.append((x, y))
        else:  # Empty line indicates a new block
            if block:
                data.append(np.array(block))
                block = []
    if block:  # Add the last block if it exists
        data.append(np.array(block))
    return data

def plot_from_gnuplot_config(results, config, ylim=None):
    fig, ax = plt.subplots()

    # Plot data from each file
    colors = ['red', 'black']
    for i, file in enumerate(config['data_files']):
        file = file.replace('.', '_')
        file_lines = results[file].get_content().split('\n')
        band_data = read_band_data(file_lines)
        for band in band_data:
            ax.plot(band[:, 0], band[:, 1], linewidth=1.5, color=colors[i])

    # Configure plot
    ax.set_xlim(config['xrange'])
    if ylim is None:
        ax.set_ylim(config['yrange'])
    else:
        ax.set_ylim(ylim)
    ax.set_ylabel(config['ylabel'])

    # Configure x-ticks
    xtick_positions = [pos for _, pos in config['xticks']]
    xtick_labels = [label for label, _ in config['xticks']]
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels)

    # Add a zero axis
    ax.axhline(0, color="black", linewidth=0.8)

    # Hide the legend
    ax.legend().set_visible(False)

    # Display the plot
    plt.tight_layout()
    plt.show()
    return fig, ax

def plot_bands(n, bandgnu_path=None, ylim=None):
    """
    This will plot the bands for a given calcultion node n. If bandgnu has not been run, it will run it, otherwise it
    will utilize previous bandgnu result.
    """

    nr = n.outputs.retrieved
    outlinks = nr.base.links.get_outgoing(node_class=CalculationFactory('core.shell')).all_nodes()

    bandgnu_node = None
    for link in outlinks:
        if 'bandgnu13' in link.base.attributes.get('process_label'):
            if 'aiida_GNUBAND' in link.outputs:
                bandgnu_node = link
                results = link.outputs

    if bandgnu_node is None:
        if bandgnu_path is None:
            raise Exception('bandgnu_path must be specified')
        results, bandgnu_node = run_bandgnu(bandgnu_path, n.outputs.retrieved)

    config = parse_gnuplot_script(results['aiida_GNUBAND'].get_content().split('\n'))
    fig, ax = plot_from_gnuplot_config(results, config, ylim=ylim)
    return fig, ax