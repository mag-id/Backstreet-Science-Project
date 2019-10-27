import ntpath
import os
import sys
from shutil import copyfile
import pandas as pd
import pysnooper

import click
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


def generate_conformations_RDkit(mol1, n, m):
    # ids - list with conformer indexes
    ids = AllChem.EmbedMultipleConfs(mol1, numConfs=n, numThreads=m)
    return mol1, list(ids)


def xyz_generator(atoms_names: tuple, coord: list) -> pd.DataFrame:
    df = pd.DataFrame(columns=list('xyz'), data=coord)
    df.insert(0, 'Atom', atoms_names)
    return df


#TODO fixit
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
        f.write(str(n_atoms) + '\n\n')
        for index, row in xyz.iterrows():
            line = '{} {} {} {}\n'.format(row['Atom'], row['x'], row['y'],
                                          row['z'])
            f.write(line)


#TODO fixit and charge
def xyz_to_mop(xyz: pd.DataFrame, file_name: str, charge=0):
    """
    Write mop file

    Arguments
    ---------
    xyz: pd.DataFrame
        contains cartesian coordinates and atoms labels
    file_name: str
        path to new file
    charge: int (0)
        molecule charge
    """
    with open(file_name, 'w+') as f:
        f.write(' AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + '\n')
        f.write('generated mop input\n\n')
        for index, row in xyz.iterrows():
            line = '{}\t{}\t1\t{}\t1\t{}\t1\n'.format(row['Atom'], row['x'],
                                                      row['y'], row['z'])
            f.write(line)
        f.write('\n')
        f.write('  OLDGEO AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + '  FORCE THERMO\n\n\n')
        f.write(' OLDGEO AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + ' COSWRT\n\n')


#TODO write it
def conformers_rmsd_matrix():
    pass


@click.command()
@click.option('--input',
              help='mol input file',
              prompt='Input file',
              type=click.Path(exists=True))
# @click.option('--smiles',
#                 help='Input smiles string',
#                 default=None,
#                 prompt='smiles string')
@click.option('--output',
              help='out files format',
              type=click.Choice(['sdf', 'xyz', 'mop'], case_sensitive=False),
              prompt='out files format')
@click.option('-n',
              help='Conformers number',
              default=50,
              prompt='Conformers number')
@click.option('-m', help='Max cores', prompt='Max cores', default=1)
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

        # Molecule optimiztaion
        i = 0
        while AllChem.MMFFOptimizeMolecule(mol1) != 0:
            i += 1
            if i >= 25:
                raise ('MMFF94 optimization had a bad end')

        pwd = os.getcwd()
        if os.path.exists(os.path.join(pwd, file_name)):
            pass
        else:
            os.mkdir(os.path.join(pwd, file_name))

        path = os.path.join(pwd, file_name, file_name)

        mol2 = Chem.RemoveHs(mol1)
        AllChem.Compute2DCoords(mol2)
        Draw.MolToFile(mol2, path + '.png', size=(300, 300))

        # Conformers generation
        mol1, ids = generate_conformations_RDkit(mol1=mol1, n=n, m=m)

        if output == 'sdf':
            sdf_output_file = path + '.sdf'
            open(sdf_output_file, 'a').close()
            writer = Chem.SDWriter(sdf_output_file)
            for i in ids:
                writer.write(mol1, confId=i)
            writer.close()

        if output == 'xyz':
            atoms = [atom.GetSymbol() for atom in mol1.GetAtoms()]
            n_atoms = len(atoms)
            for i in range(len(ids)):
                file_name = path + f'{i}.xyz'
                coord = xyz_generator(
                    atoms_names=atoms,
                    coord=mol1.GetConformer(i).GetPositions())
                write_xyz(n_atoms=n_atoms, xyz=coord, file_name=file_name)

        if output == 'mop':
            atoms = [atom.GetSymbol() for atom in mol1.GetAtoms()]
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
