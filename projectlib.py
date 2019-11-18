
# TODO: write doc for initialization funcs

import pathlib
import pybel
import json


def deconstruct(path: str) -> (str, str, str):
    """"""
    pop = pathlib.Path(path)

    folder = pop.parent
    stem = pop.stem
    suffix = pop.suffix

    return str(folder), stem, suffix


def construct(folder: str = pathlib.Path.cwd(), stem: str = '', suffix: str = '') -> str:
    """"""
    folder = pathlib.Path(folder)
    name = pathlib.Path(stem + suffix)

    return str(pathlib.Path.joinpath(folder, name))


def create(path: str) -> True:
    """"""
    path = pathlib.Path(path)
    pathlib.Path.mkdir(path, parents=True, exist_ok=False)
    return True


def convert(inp_type: str, inp_path: str, out_type: str, out_path: str) -> True:
    """"""
    inp = next(pybel.readfile(inp_type, inp_path))
    out = pybel.Outputfile(out_type, out_path, overwrite=False)
    out.write(inp)
    return True


def confslice(inp_type: str, inp_path: str, out_type: str, out_path: str) -> True:
    """"""
    for number, structure in enumerate(pybel.readfile(inp_type, inp_path)):

        path, stem, suffix = deconstruct(out_path)
        new_path = construct(path, f'{stem}_{number}', suffix)

        conformer = pybel.Outputfile(out_type, new_path, overwrite=False)
        conformer.write(structure)

    return True


def json_save(path: str, dictionary: dict) -> True:
    """"""
    with open(path, 'a') as f:
        json.dump(dictionary, f)
    return True


def json_load(path: str) -> dict:
    """"""
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary


# TODO: to refactor out-to-mop and clean-and-rewrite funcs

import pandas
import copy
import os
import re


def find(directory, slash='/', pattern=r'.+\.out'):
    """
    Function to yield path for my_test.out file(s) by default.

    :param directory: path to the start directory, string
    :param slash: the path delimiter, '/' by default, string
    :param pattern: regular expression for parsing file names, was determined by default, raw string
    :return: string containing path to file like 'directory/subdirectory/my_test.out'
    """
    for directory, subdirectories, files in os.walk(directory):
        for file in files:
            if re.findall(pattern, str(file)):
                yield str(directory + slash + file)


def extract(path,
            top_flag='first_line',
            bottom_flag='last_line',
            pattern=r'([\d]*)\s+([a-zA-Z]*)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)',
            columns=('num', 'atom_name', 'x', 'y', 'z')):
    """
    Function is extract cartesian coordinates.

    Returns pandas DataFrame with columns where:
    'num' is atom number
    'atom_name' is chemical symbol
    'x', 'y', 'z' are Cartesian coordinates respectively

    :param path: path to file, string
    :param top_flag: string(s) above of the last cartesian block,
    used 'first_line' by default takes first line in file, string
    :param bottom_flag: string(s) below of the last cartesian block,
    used 'last_line' by default takes last line in file, string
    :param pattern: regular expression for parsing lines in cartesian block, was determined by default, raw string
    :param columns: the columns to passed for pandas DataFrame formation, tuple with strings
    :return: pandas DataFrame
    """
    top = int()
    bottom = int()
    rows = list()

    # flags counter
    with open(path, 'r') as p:
        for num, line in enumerate(p, 1):

            # top_flag block
            if top_flag is 'first_line':
                top = 1
            elif top_flag in line:
                top = num

            # bottom_flag block
            if bottom_flag is 'last_line':
                bottom = num
            elif bottom_flag in line:
                bottom = num

    # strings extraction
    with open(path, 'r') as p:
        for num, line in enumerate(p, 1):
            if (top <= num) and (bottom >= num):
                rows += re.findall(pattern, line)

    # pandas DataFrame
    return pandas.DataFrame(rows, columns=columns)


def save(cartesian, charge, name, comment='',
         mop_dir=str(os.getcwd()), csv_dir=str(os.getcwd()),
         mop_file=True, csv_file=True):
    """
    Function write the pandas DataFrame into the MOPAC input file .mop format
    and into .csv file for the subsequent MOPAC calculations

    Write 'name'.mop and 'name'.csv files into 'mop_dir' and 'csv_dir' directories

    Returns nothing

    :param cartesian: pandas DataFrame from extract() function
    :param charge: specify charge of the molecule, string
    :param name: name of the file from named() function, string
    :param comment: add comment into .mop input file, empty by default, string
    :param mop_dir: path to the directory to save .mop files, current work directory by default, string
    :param csv_dir: path to the directory to save .csv files, current work directory by default, string
    :param mop_file: save .mop file, True by default, bool
    :param csv_file: save .csv file, True by default, bool
    :return: nothing
    """
    # Preparation
    del cartesian['num']
    coordinates = 'x', 'y', 'z'
    columns = {'1_1': 2, '1_2': 4, '1_3': 6}
    for coord in coordinates:
        cartesian[coord] = cartesian[coord].astype(float)

    # .mop block
    if mop_file is True:
        cartesian = cartesian.round(decimals=5)
        for column in columns:
            cartesian.insert(columns[column], column, 0, True)

        # .mop constructor
        with open(mop_dir + 'sp_' + name + '.mop', 'a') as mop:
            mop.writelines('AUX LARGE CHARGE=' + charge + ' SINGLET NOOPT PM6-DH2X' + '\n')
            mop.writelines(comment + '\n'*2)
            mop.writelines(cartesian.to_string(header=False, index=False))
            mop.writelines('\n'*2)

    # .csv block
    if csv_file is True:
        cartesian = cartesian.round(decimals=4)
        for column in columns:
            if mop_file is True:
                del cartesian[column]
            cartesian.insert(columns[column], column, 1, True)
        cartesian.to_csv(csv_dir + name + '.csv')


def softcheck(template=list(), compared=list(), logging_mode=False):
    """
    Function is check elements of the template's list and compared list,
    returns list of the tuples with 'equality' (str) and 'index' (int) where:

    logging_mode=True print all cases in format 'index, element, equality, element,' where equality means:
    '<<' that the compared element is redundant and template element is absent
    '>>' the template element is redundant, while compared element is absent
    '!=' the elements are different
    '==' the elements are the same

    :param template: pandas DataFrame
    :param compared: pandas DataFrame
    :param logging_mode: print log
    :return: list of the tuples with string and integer
    """
    template = copy.deepcopy(template)
    compared = copy.deepcopy(compared)

    if len(template) >= len(compared):
        indexes = range(len(template))
    else:
        indexes = range(len(compared))

    errors = list()

    for index in indexes:
        try:
            if template[0] != compared[0]:
                equality = '!='
                template_label = template.pop(0)
                compared_label = compared.pop(0)

            else:
                equality = '=='
                template_label = template.pop(0)
                compared_label = compared.pop(0)

        except IndexError:

            if len(template) > len(compared):
                equality = '>>'
                template_label = template.pop(0)
                compared_label = '*'

            else:
                equality = '<<'
                template_label = '*'
                compared_label = compared.pop(0)

        if equality != '==':
            errors.append((equality, index))

        if logging_mode is True:
            print('{} {} {} {}'.format(index, template_label, equality, compared_label))

    if errors is None:
        return True
    else:
        return errors


def hardcheck(template=list(), compared=list(), logging_mode=False):
    """
    Function is check elements of the template's list and compared list,
    returns index of the missing element of the template's list, int.

    Raises IndexError when compared list has less element, than template.
    Raises ValueError when element in lists are different.

    logging_mode=True print all cases in format 'index, element, equality, element,' where equality means:
    '<<' that the compared element is redundant and template element is absent
    '>>' the template element is redundant, while compared element is absent
    '!=' the elements are different
    '==' the elements are the same

    :param template: pandas DataFrame
    :param compared: pandas DataFrame
    :param logging_mode: print log, disable ValueError and IndexError
    :return: index, int
    """

    template_label = str
    compared_label = str
    equality = str

    if len(template) >= len(compared):
        indexes = range(len(template))
    else:
        indexes = range(len(compared))

    for index in indexes:
        try:
            template_label = template[index]
            compared_label = compared[index]

            if template_label != compared_label:
                equality = '!='

            else:
                equality = '=='

        except IndexError:

            if len(template) > len(compared):
                equality = '>>'
                template_label = template[index]
                compared_label = '*'

            if len(template) < len(compared):
                equality = '<<'
                template_label = '*'
                compared_label = compared[index]

        message = '{} {} {} {}'.format(index, template_label, equality, compared_label)

        if logging_mode is True:
            print(message)
        else:
            if equality == '<<':
                yield index
            elif equality == '>>':
                raise IndexError('less labels at ' + message + ', try logging_mode=True')
            elif equality == '!=':
                raise ValueError('different labels at ' + message + ', try logging_mode=True')


def atomdrop(template=pandas.DataFrame(), compared=pandas.DataFrame(), logging_mode=False):
    """
        Function is check elements of the template's list and compared list,
        returns list with indexes (int) of the conflicting elements,
        also can returns list with indexes (int) and False (bool) within,
        if compared list has less element, than template.

        logging_mode=True print all cases in format 'element, equality, element, template's index',
        where equality means:
            '<<' that the compared element is redundant and template element is absent
            '>>' the template element is redundant, while compared element is absent
            '!=' the elements are different
            '==' the elements are the same

        :param template: pandas DataFrame
        :param compared: pandas DataFrame
        :param logging_mode: print log, disable ValueError
        :return: list with integers and bool within
        """
    # Initialize lists
    template_atoms = list(template['atom_name'])
    compared_atoms = list(compared['atom_name'])
    errorlist = list()

    # Define count
    if len(template_atoms) >= len(compared_atoms):
        indexes = len(template_atoms)
    else:
        indexes = len(compared_atoms)

    # pop from stacks
    for index in range(indexes):
        try:
            if template_atoms[0] == compared_atoms[0]:
                arg_one = template_atoms.pop(0)
                arg_two = '=='
                arg_three = compared_atoms.pop(0)
            else:
                arg_one = template_atoms[0]
                arg_two = '!='
                arg_three = compared_atoms.pop(0)
                errorlist.append(index)

        # if template > compared list: yield False without IndexError
        except IndexError:
            if len(template_atoms) > len(compared_atoms):
                arg_one = template_atoms.pop(0)
                arg_two = '>>'
                arg_three = '*'
                errorlist.append(False)
            else:
                arg_one = '*'
                arg_two = '<<'
                arg_three = compared_atoms.pop(0)
                errorlist.append(index)

        # log message
        if logging_mode is True:
            print('{}\t{} {} {}'.format(index, arg_one, arg_two, arg_three))

    return errorlist


def rewrite(cartesian, path, pattern):
    """
    Function rewrites .mop file using given cartesian coordinates

    :param cartesian: cartesian: pandas DataFrame from extract() function
    :param path: path to file, string
    :param pattern: regular expression for parsing lines in cartesian block, raw string
    :return: nothing
    """
    # Preparation
    coordinates = 'x', 'y', 'z'
    columns = {'1_1': 2, '1_2': 4, '1_3': 6}
    for coord in coordinates:
        cartesian[coord] = cartesian[coord].astype(float)

    cartesian = cartesian.round(decimals=4)
    for column in columns:
        cartesian.insert(columns[column], column, 1, True)

    # rewrite .mop
    flag = True
    with open(path, 'r') as old, open(path + '_tmp', 'a') as new:
        for line in old:
            if re.match(pattern, line) is None:
                new.write(line)
            elif re.match(pattern, line) and flag is True:
                new.writelines(cartesian.to_string(header=False, index=False))
                new.write('\n')
                flag = False

    # Rewrite
    os.remove(path)
    os.rename(path + '_tmp', path)


def named(path, pattern):
    """
    Function takes path to file and return its name.

    :param path: path to file, string
    :param pattern: regular expression for parsing file names, raw string
    :return: name of the file, string
    """
    return re.findall(pattern, path)[0]


def charged(path, pattern=r'.*CHARGE ON SYSTEM = +\+?(\-?\d+).*'):
    """
    Function takes path to .out file and return molecule charge

    :param path: path to file, string
    :param pattern: regular expression for parsing molecule charge value, was determined by default, raw string
    :return: last element of the list as string
    """
    with open(path, 'r') as p:
        charges = list(re.findall(pattern, (p.read())))

    return str(charges[-1])


def commented(path):
    """
    Function takes path to .out file and return comment

    :param path: path to file, string
    :return: comment, string
    """
    return 'Conformer ' + re.findall(r'.*test(\d+)\.out', path)[0]
