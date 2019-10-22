import ntpath
import os
import sys
from shutil import copyfile
import pandas as pd
import pysnooper

import click
from rdkit import Chem
from rdkit.Chem import AllChem


def generate_conformations_RDkit(mol1, n, m):
    ids = AllChem.EmbedMultipleConfs(mol1,
                                     numConfs=n,
                                     numThreads=m)
    return mol1, list(ids)


def xyz_generator(atoms_names: tuple, coord: list) -> pd.DataFrame:
    df = pd.DataFrame(columns=list('xyz'), data=coord)
    df.insert(0, 'Atom', atoms_names)
    return df


def write_xyz(xyz: pd.DataFrame, file_name: str, n_atoms: int):
    """
    Write XYZ file

    Arguments
    ---------
    xyz: pd.DataFrame contains cartesian coordinates and atoms labels
    file_name: path to new file
    n_atoms: molecules len
    """
    with open(file_name, 'w+') as f:
        f.write(str(n_atoms) + '\n')
        # f.writelines([
        #     '\n',
        #     str(xyz['Atom']), ' ',
        #     str(xyz['x']), ' ',
        #     str(xyz['y']), ' ',
        #     str(xyz['z'])
        # ])


def xyz_to_mop(xyz: pd.DataFrame, file_name: str, charge=0):
    with open(file_name, 'w+') as f:
        f.write(' AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + '\n')
        _ = pd.Series([1 for i in range(len(xyz.loc[:, ('z')]))])
        # f.writelines([
        #     '\n',
        #     str(xyz['Atom']), ' ',
        #     str(xyz['x']), ' ',
        #     str(_), ' ',
        #     str(xyz['y']), ' ',
        #     str(_), ' ',
        #     str(xyz['z']), ' ', 
        #     str(_)
        # ])
        f.write('\n')
        f.write('  OLDGEO AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + '  FORCE THERMO\n\n\n')
        f.write(' OLDGEO AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + ' COSWRT\n')


def conformers_rmsd_matrix():
    pass


"""@click.command()
@click.argument('smi_input_file', nargs=-1)
@click.argument('sdf_output_file')
@click.option('-n', default=1, help='number of conformers')
def main(smi_input_file, sdf_output_file, n):
    n = int(n)
    #TODO Исправить
    print(sdf_output_file, smi_input_file)
    f = open(sdf_output_file, "w+")
    f.close()

    writer = Chem.SDWriter(sdf_output_file)

    # suppl = Chem.SmilesMolSupplier(smi_input_file, titleLine=False)

    with futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit a set of asynchronous jobs
        jobs = []
        for mol in suppl:
            if mol:
                job = executor.submit(generate_conformations, mol, n)
                jobs.append(job)

        widgets = [
            "Generating conformations; ",
            progressbar.Percentage(), " ",
            progressbar.ETA(), " ",
            progressbar.Bar()
        ]
        pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(jobs))
        for job in pbar(futures.as_completed(jobs)):
            mol, ids = job.result()
            for id in ids:
                writer.write(mol, confId=id)"""
@click.command()
@click.option('--input', help='mol input file', prompt='Input file')
# @click.option('--smiles',
#                 help='Input smiles string',
#                 default=None,
#                 prompt='smiles string')
@click.option('--output',
                help='out files format',
                default='sdf',
                prompt='out files format')
@click.option('-n',
                help='Conformers number',
                default=50,
                prompt='Conformers number')
@click.option('-m', help='Max cores', default=1)
@pysnooper.snoop()
def main(input, output, n, m):
    file_name = 'test'
    if input or smiles:
        if input:
            file_name = ntpath.basename(input).split('.')[0]
            mol = Chem.MolFromMolFile(input)
        # smiles path
        else:
            file_name = smiles
            mol = Chem.AtomFromSmiles()
        # Add Hydrogens
        mol1 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol1)
        i = 0
        while AllChem.MMFFOptimizeMolecule(mol1) != 0:
            i += 1
            if i >= 25:
                raise ('MMFF94 optimization had a bad end')

        mol1, ids = generate_conformations_RDkit(mol1=mol1, n=n, m=m)

        pwd = os.getcwd()
        if os.path.exists(os.path.join(pwd, file_name)):
            pass
        else:
            os.mkdir(os.path.join(pwd, file_name))

        if output == 'sdf':
            sdf_output_file = os.path.join(pwd, file_name, file_name) + '.sdf'
            open(sdf_output_file, 'a').close()
            writer = Chem.SDWriter(sdf_output_file)
            for i in ids:
                writer.write(mol1, confId=i)
            writer.close()

        if output == 'xyz':
            atoms = [atom.GetSymbol() for atom in mol1.GetAtoms()]
            n_atoms = len(atoms)
            path = os.path.join(pwd, file_name, file_name)
            for i in range(len(ids)):
                file_name = path + f'{i}.xyz'
                coord = xyz_generator(
                    atoms_names=atoms,
                    coord=mol1.GetConformer(i).GetPositions())
                write_xyz(n_atoms=n_atoms, xyz=coord, file_name=file_name)

        if output == 'mop':
            atoms = [atom.GetSymbol() for atom in mol1.GetAtoms()]
            path = os.path.join(pwd, file_name, file_name)
            for i in range(len(ids)):
                file_name = path + f'{i}.mop'
                coord = xyz_generator(
                    atoms_names=atoms,
                    coord=mol1.GetConformer(i).GetPositions())
                xyz_to_mop(xyz=coord, file_name=file_name)
    else:
        raise ('Please choose a input file or smiles')


if __name__ == '__main__':
    main()
    """
mol = ....
mol1 = Chem.AddHs(mol)
AllChem.MMFFOptimizeMolecule(mol1)
cids = AllChem.EmbedMultipleConfs(mol1, numConfs=500, numThreads=3)
ids = list(cids)

sdf_output_file = '/home/anton/Documents/Backstreet-Science-Project/test/data/asd.sdf'
writer = Chem.SDWriter(sdf_output_file)
for i in ids:
    writer.write(mol1, confId=i)
writer.close()

# GET XYZ for conformer
for i, j in list(zip(atoms ,mol1.GetConformer(1).GetPositions())):
    print(i, j)
    """