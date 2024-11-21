from aiida_openmx.input.flat import replace_backslash

def write_jx_input(filename,folder,parameters,ij_pairs):

    parameters = replace_backslash(parameters)
    if 'num.ij.pairs' not in parameters:
       parameters['num.ij.pairs'] = len (ij_pairs)

    with folder.open(filename, 'w') as handle:
        for key in parameters:
            value = parameters[key]
            handle.write(f"{key:<35} {value:<35}\n")

        handle.write('<ijpairs.cellid\n')
        for ij in ij_pairs:
            ij_str = ' '.join([str(i) for i in ij])
            handle.write(ij_str + '\n')
        handle.write('ijpairs.cellid>\n')