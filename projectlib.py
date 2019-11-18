
import os
import re
import pandas


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


def hardcheck(template=list(), compared=list(), logging_mode=False):
    """
    Function is check elements of the template's list and compared list,
    returns index of the missing element of the template's list, int.

    Raises ValueError when:
    1) compared list has less element, than template
    2) element in lists are different

    logging_mode=True print all cases in format 'index, element, equality, element,' where equality means:
    '<<' that the compared element is redundant and template element is absent
    '>>' the template element is redundant, while compared element is absent
    '!=' the elements are different
    '==' the elements are the same

    :param template: pandas DataFrame
    :param compared: pandas DataFrame
    :param logging_mode: print log, disable ValueError
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


def named(path):
    """
    Function takes path to .out file and return its name

    :param path:  path to file, string
    :return: name of the file, string
    """
    return re.findall(r'.*(test\d+)\.out', path)[0]


def commented(path):
    """
    Function takes path to .out file and return comment

    :param path: path to file, string
    :return: comment, string
    """
    return 'Conformer ' + re.findall(r'.*test(\d+)\.out', path)[0]
