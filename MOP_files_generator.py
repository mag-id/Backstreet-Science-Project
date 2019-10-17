import glob
import os
import pandas as pd
from tqdm import tqdm
import datetime
import click
import tempfile


'''
Функция принимает на вход путь до файла mol2, 
создает генератор, возвращающий название следующей молекулы
'''
def conformer_detector(way):
    with open(way, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if '@<TRIPOS>MOLECULE' in lines[i]:
                conformer = (lines[i + 1])
                i += 1
                yield (conformer)

'''
Функция принимает на вход путь до файла mol2, 
создает генератор, возвращающий массив кортежей
с названием и коорднатами атомов каждой молекулы
'''
def structure_detector(way):
    with open(way, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if '@<TRIPOS>ATOM' in lines[i]:
                i += 1
                lines_stock = []
                while '@<TRIPOS>BOND' not in lines[i]:
                    cur_coord = (lines[i].split())
                    if len(cur_coord) > 0:
                        cur_data = (cur_coord[1][0], float(cur_coord[2]),
                                    float(cur_coord[3]), float(cur_coord[4]))
                        lines_stock.append(cur_data)
                    i += 1
                yield (lines_stock)


def get_files(folder: str):
    """Получает список *.mol2 файлов в данной папке"""
    return [files for files in glob.glob(folder + '/*.mol2')]


def new_folder_create(path: str):
    """
    Получает на вход путь и создает на его основое новую директорию
    в которую, в дальнейшем, заливаем сгенерированные mop файлы
    """
    now = datetime.datetime.now()
    new_path = path + '/MOP_inp_gen/' + str(now.strftime("%d-%m-%Y_%H-%M"))
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    return new_path


def mop_writer(new_path: str, data: pd.DataFrame, i: int,
               atom_coordinates: list, charge: str):
    """
    Получает на вход параметры, и записывает mop файл

    Arguments
    ---------
    new_path : str
        Путь к новому файлу
    data: pd.DataFrame
        Данные из файла mol2 в виде DataFrame
    i: int
        Порядковый номер конформера
    atom_coordinates: list
        xyz list
    charge: int
        Заряд молекулы
    """
    charge = charge
    trashhold = 0.5
    with open(new_path + '/test' + str(i) + '.mop', 'w+') as d:
        d.write(' AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + '\n')
        d.write('Conformer ' + str(i) + '\n\n')
        for index, row in data.iterrows():
            n_row = [str(j) for j in row]
            for crd in range(len(atom_coordinates)):
                distance = ((row['x'] - atom_coordinates[crd][0])**2 +
                            (row['y'] - atom_coordinates[crd][1])**2 +
                            (row['z'] - atom_coordinates[crd][2])**2)
                if distance >= trashhold**2 or distance == 0:
                    if crd == len(atom_coordinates) - 1:
                        d.write("{}{}{}".format(' ' * 3, '\t'.join(n_row),
                                                '\n'))
                        pass
                    else:
                        continue
                else:
                    break
        d.write('\n')
        d.write('  OLDGEO AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + '  FORCE THERMO\n\n\n')
        d.write(' OLDGEO AUX LARGE COSCCH NSPA=92 EPS=78.4 PM6-DH2X CHARGE=' +
                str(charge) + ' COSWRT\n')
    data.to_csv(str(new_path) + '/test' + str(i) + '.csv')


def mol2_prereader(file: str):
    """Исправляет косяк и записывает tmp файл, возвращает путь к временному, исправленному, файлу"""
    with tempfile.NamedTemporaryFile(delete=False, mode='w+') as tmp:
        with open(file, 'r') as f:
            for line in f.readlines():
                if '1<O>11' in line:
                    tmp.write(line.replace('1<O>11', '1 UNK'))
                else:
                    tmp.write(line)
    return (tmp.name, new_folder_create(file[:-5]))


def mol2_reader(temp_file: str, new_path: str):
    """
    Arguments
    ---------
    
    """
    i = 0
    # mol2_list = list(split_multimol2(str(temp_file)))
    mol2_list = list(
        zip(conformer_detector(temp_file), structure_detector(temp_file)))
    charge = input('Set charge of the molecule:\n')
    print('\n' + 'Parsing file in progress...\n')
    tt = tqdm(total=len(mol2_list))
    for name, mol2 in mol2_list:
        try:
            # cur_molecule = PandasMol2().read_mol2_from_list(mol2_lines=mol2[1],
            #                                                 mol2_code=mol2[0])
            atom_coordinates = []
            data = pd.DataFrame(columns=('atom_name', 'x', 'y', 'z'),
                                data=mol2)
            data.astype(float, errors='ignore', copy=False)
            # data = cur_molecule.df[['atom_name', 'x', 'y', 'z']]
            for index, row in data.iterrows():
                atom_coordinates.append((row['x'], row['y'], row['z']))
            data.loc[:,
                     ('atom_name')] = data['atom_name'].apply(lambda x: x[0])
            data.insert(2, "1_1",
                        [1 for i in range(len(data.loc[:, ('atom_name')]))],
                        True)
            data.insert(4, "1_2",
                        [1 for i in range(len(data.loc[:, ('atom_name')]))],
                        True)
            data.insert(6, "1_3",
                        [1 for i in range(len(data.loc[:, ('atom_name')]))],
                        True)
            i += 1
            mop_writer(new_path=new_path,
                       data=data,
                       i=i,
                       atom_coordinates=atom_coordinates,
                       charge=charge)
        except Exception as e:
            print(e)
            print(mol2)
            pass
        finally:
            tt.update(n=1)


def pdb_prereader(file: str):
    pass


def pdb_reader(file: str):
    pass


@click.command()
@click.argument('filename', type=click.Path(exists=True))
def main(filename=None):
    if filename:
        temp_file, new_path = mol2_prereader(filename)
        mol2_reader(temp_file=temp_file, new_path=new_path)
        os.unlink(temp_file)


if __name__ == '__main__':
    main()

#%%
