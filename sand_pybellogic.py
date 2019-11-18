
import pybel
import subprocess


def convert(inp_type: str, inp_path: str, out_type: str, out_path: str) -> True:
    """"""
    inp = next(pybel.readfile(inp_type, inp_path))
    out = pybel.Outputfile(out_type, out_path, overwrite=False)
    out.write(inp)
    return True


def geneticsearch(inp_file: str, out_file: str, n: int = 0) -> subprocess.run:
    """"""
    """
    OpenBabel doc: One of the ops
    conformer    Conformer Searching (not displayed in GUI)
    Typical usage: obabel infile.xxx -O outfile.yy --conformer --nconf
    options:             description
    --log            output a log of the energies (default = no log)
    --nconf #        number of conformers to generate
    forcefield based methods for finding stable conformers:
    --systematic     systematically generate all conformers
    --fast           fast systematic search from central bonds
    --random         randomly generate conformers
    --weighted       weighted rotor search for lowest energy conformer
    --ff #           select a forcefield (default = MMFF94)
    --rings          sample ring torsions
    genetic algorithm based methods (default):
    --children #     number of children to generate for each parent (default = 5)
    --mutability #   mutation frequency (default = 5)
    --converge #     number of identical generations before convergence is reached
    --score #        scoring function [rmsd|energy|minrmsd|minenergy] (default = rmsd)
    """
    command = ['obabel', inp_file, '-O', out_file, '--conformer', f'--nconf {n}', '--writeconformers']

    return subprocess.run(command)


'''
def confabserch(inp_file: str, out_file: str) -> subprocess.run:
    """"""
    """
    OpenBabel doc: One of the ops
    confab    Confab, the diverse conformer generator
    Typical usage: obabel infile.xxx -O outfile.yyy --confab --conf 1000000
    options:
    --conf #     Max number of conformers to test (default is 1 million)
    --rcutoff #  RMSD cutoff (default 0.5 Angstrom)
    --ecutoff #  Energy cutoff (default 50.0 kcal/mol)
    --original   Include the input conformation as the first conformer
    --verbose    Verbose output
    """
    command = ['obabel', inp_file, '-O', out_file, '--confab', '--writeconformers']

    return subprocess.run(command)


def confabreport(out_file: str, inp_file: str, n: float) -> subprocess.run:
    """"""
    command = ['obabel', out_file, '-oconfabreport', '-xf', inp_file, '-xr', str(n)]

    return subprocess.run(command)


def read(file_type: str, file_path: str):
    """"""
    return pybel.readfile(file_type, file_path)
'''
