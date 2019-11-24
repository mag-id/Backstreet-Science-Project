
from os import walk, remove, rename
from copy import deepcopy

import pathlib
import pandas
import pybel
import json
import re


def deconstruct(path: str) -> (str, str, str):
    """
    The function deconstructs the path to the file into folder, stem, and suffix and returns them as strings.

    :param path: a path to file, string
    :return: str(folder), str(stem), str(suffix)
    """
    pop = pathlib.Path(path)

    folder = pop.parent
    stem = pop.stem
    suffix = pop.suffix

    return str(folder), stem, suffix


def construct(folder: str = str(pathlib.Path.cwd()), stem: str = '', suffix: str = '') -> str:
    """
    The function constructs the path to the file from the path to the folder as string, stem, and suffix.
    By default folder is defined as str(pathlib.Path.cwd()), stem and suffix as empty strings ('').
    Returns path to file as a string.

    :param folder: a path to folder, str
    :param stem: name of a file, str
    :param suffix: suffix of a file, str
    :return: str
    """
    folder = pathlib.Path(folder)
    name = pathlib.Path(stem + suffix)

    return str(pathlib.Path.joinpath(folder, name))


def create(*paths: str) -> True:
    """
    The function creates the path.

    :param paths: path(s), str
    :return: True, bool
    """
    for path in paths:
        path = pathlib.Path(path)
        pathlib.Path.mkdir(path, parents=True, exist_ok=False)

    return True


def convert(inp_type: str, inp_path: str, out_type: str, out_path: str) -> True:
    """
    The function takes the file according to the OpenBabel input type from the input path
    and converts it into the file according to the OpenBabel output type into output path.

    :param inp_type: type of the input file according to the OpenBabel, str
    :param inp_path: path to the input file, str
    :param out_type: type of the output file according to the OpenBabel, str
    :param out_path: path to the writing of the output file
    :return: True
    """
    inp = next(pybel.readfile(inp_type, inp_path))
    out = pybel.Outputfile(out_type, out_path, overwrite=False)
    out.write(inp)

    return True


def confslice(inp_type: str, inp_path: str, out_type: str, out_path: str) -> True:
    """
    The function slices the file of the OpenBabel type with multiple structures
    to the new files of the OpenBabel type with single structures.
    New files are named according to stem_number.suffix rule.

    :param inp_type: type of the input file according to the OpenBabel, str
    :param inp_path: a path to the input file, str
    :param out_type: type of the output files according to the OpenBabel, str
    :param out_path: a path to the writing of the output files
    :return: True
    """
    for number, structure in enumerate(pybel.readfile(inp_type, inp_path)):

        path, stem, suffix = deconstruct(out_path)
        new_path = construct(path, f'{stem}_{number}', suffix)

        conformer = pybel.Outputfile(out_type, new_path, overwrite=False)
        conformer.write(structure)

    return True


def json_save(path: str, dictionary: dict) -> True:
    """
    The function saves python dictionary at JSON format to the file according to the path, str.

    :param path: the path to the file, str
    :param dictionary: dict
    :return: True
    """
    with open(path, 'a') as f:
        json.dump(dictionary, f)

    return True


def json_load(path: str) -> dict:
    """
    The function returns Python dictionary from JSON.

    :param path: the path to the JSON file, str
    :return: dict
    """
    with open(path) as f:
        dictionary = json.load(f)

    return dictionary


def find(directory: str, slash: str = '/', pattern: str = r'.+\.out') -> str:
    """
    The function to yield path for my_test.out file(s) by default.

    :param directory: path to the start directory, string
    :param slash: the path delimiter, '/' by default, string
    :param pattern: regular expression for parsing file names, was determined by default, raw string
    :return: string containing path to file like 'directory/subdirectory/my_test.out'
    """

    for directory, subdirectories, files in walk(directory):
        for file in files:
            if re.findall(pattern, str(file)):
                yield str(directory + slash + file)


def xyzextract(
        path: str, start: str = 'start', end: str = 'end',
        pattern: str = r'([\d]*)\s+([a-zA-Z]*)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)',
        columns: tuple = ('num', 'atom_name', 'x', 'y', 'z')
            ) -> pandas.DataFrame or str:
    """
    The function is extract cartesian coordinates.

    Returns pandas DataFrame with columns where:
    'num' is atom number
    'atom_name' is chemical symbol
    'x', 'y', 'z' are Cartesian coordinates respectively
    Returns string with path if start or stop parameters do not found

    :param path: path to file, string
    :param start: string(s) above of the last cartesian block,
    used 'start' by default takes first line in file, string
    :param end: string(s) below of the last cartesian block,
    used 'end' by default takes last line in file, string
    :param pattern: regular expression for parsing lines in cartesian block, was determined by default, raw string
    :param columns: the columns to passed for pandas DataFrame formation, tuple with strings
    :return: pandas DataFrame or string
    """
    # flags counter
    top, bottom = int, int
    with open(path, 'r') as p:
        for num, line in enumerate(p, 1):

            # top_flag
            if start == 'start':
                top = 1
            elif start in line:
                top = num
            else:
                return f'start={start} not found for {path}'

            # bottom_flag
            if end == 'end':
                bottom = num
            elif end in line:
                bottom = num
            else:
                return f'start={end} not found for {path}'

    # strings extraction
    with open(path, 'r') as p:
        rows = list()
        for num, line in enumerate(p, 1):
            if (top <= num) and (bottom >= num):
                rows.append(re.findall(pattern, line))

    # pandas DataFrame
    return pandas.DataFrame(rows, columns=columns)


def xyz2mop(
        dataframe: pandas.DataFrame, charge: str, name: str, comment: str = '',
        mop_dir: str = str(pathlib.Path.cwd()), csv_dir: str = str(pathlib.Path.cwd()),
        mop_file=True, csv_file=True
        ) -> True:
    """
    The function write the pandas DataFrame into the MOPAC input file .mop format
    and into .csv file for the subsequent MOPAC calculations

    Write 'name'.mop and 'name'.csv files into 'mop_dir' and 'csv_dir' directories

    :param dataframe: pandas DataFrame from extract() function
    :param charge: specify charge of the molecule, string
    :param name: name of the file from named() function, string
    :param comment: add comment into .mop input file, empty by default, string
    :param mop_dir: path to the directory to save .mop files, current work directory by default, string
    :param csv_dir: path to the directory to save .csv files, current work directory by default, string
    :param mop_file: save .mop file, True by default, bool
    :param csv_file: save .csv file, True by default, bool
    :return: True
    """
    # Preparation
    del dataframe['num']
    coordinates = 'x', 'y', 'z'
    columns = {'1_1': 2, '1_2': 4, '1_3': 6}
    for coord in coordinates:
        dataframe[coord] = dataframe[coord].astype(float)

    # .mop block
    if mop_file is True:
        dataframe = dataframe.round(decimals=5)
        for column in columns:
            dataframe.insert(columns[column], column, 0, True)

        # .mop constructor
        with open(mop_dir + 'sp_' + name + '.mop', 'a') as mop:
            mop.writelines('AUX LARGE CHARGE=' + charge + ' SINGLET NOOPT PM6-DH2X' + '\n')
            mop.writelines(comment + '\n'*2)
            mop.writelines(dataframe.to_string(header=False, index=False))
            mop.writelines('\n'*2)

    # .csv block
    if csv_file is True:
        cartesian = dataframe.round(decimals=4)
        for column in columns:
            if mop_file is True:
                del cartesian[column]
            cartesian.insert(columns[column], column, 1, True)
        cartesian.to_csv(csv_dir + name + '.csv')

    return True


def softcheck(template: list, compared: list, logging_mode=False) -> True or list:
    """
    The function is check elements of the template's list and compared list,
    returns True or list of the tuples with 'equality' (str) and 'index' (int) where:

    logging_mode=True print all cases in format 'index, element, equality, element,' where equality means:
    '<<' that the compared element is redundant and template element is absent
    '>>' the template element is redundant, while compared element is absent
    '!=' the elements are different
    '==' the elements are the same

    :param template: pandas DataFrame
    :param compared: pandas DataFrame
    :param logging_mode: print log
    :return: True or list of the tuples with string and integer
    """
    template = deepcopy(template)
    compared = deepcopy(compared)

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
            print(f'{index} {template_label} {equality} {compared_label}')

    if errors is None:
        return True
    else:
        return errors


def hardcheck(template: list, compared: list, logging_mode=False) -> str or IndexError or ValueError:
    """
    The function is check elements of the template's list and compared list,
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

        message = f'{index} {template_label} {equality} {compared_label}'

        if logging_mode is True:
            print(message)
        else:
            if equality == '<<':
                yield index
            elif equality == '>>':
                raise IndexError(f'less labels at {message} , try logging_mode=True')
            elif equality == '!=':
                raise ValueError(f'different labels at {message} , try logging_mode=True')


def atomdrop(template: pandas.DataFrame(), compared: pandas.DataFrame(), logging_mode=False) -> list:
    """
    The function is check elements of the template's list and compared list,
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
            print(f'{index}\t{arg_one} {arg_two} {arg_three}')

    return errorlist


def xyzrewrite(cartesian: pandas.DataFrame, path: str, pattern: str) -> True:
    """
    The function rewrites .mop file using given cartesian coordinates

    :param cartesian: cartesian: pandas DataFrame from extract() function
    :param path: path to file, string
    :param pattern: regular expression for parsing lines in cartesian block, raw string
    :return: True
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
    remove(path)
    rename(path + '_tmp', path)

    return True


def charged(path: str, pattern: str = r'.*CHARGE ON SYSTEM = +\+?(\-?\d+).*') -> str:
    """
    The function takes path to .out file and return molecule charge

    :param path: path to file, string
    :param pattern: regular expression for parsing molecule charge value, was determined by default, raw string
    :return: last element of the list as string
    """
    with open(path, 'r') as p:
        charges = list(re.findall(pattern, (p.read())))

    return str(charges[-1])
