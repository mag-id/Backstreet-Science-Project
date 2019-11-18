import os
import re
import pandas


def extract(path,
            top_flag='CARTESIAN COORDINATES',
            bottom_flag='Empirical Formula:',
            pattern=r'([\d]*)\s+([a-zA-Z]*)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)'):
    """
    Function is extract cartesian coordinates from last block of the MOPAC output file by default.

    Returns pandas DataFrame with columns where:
    'num' is atom number
    'atom_name' is chemical symbol
    'x', 'y', 'z' are Cartesian coordinates respectively

    :param path: path to file, string
    :param top_flag: string(s) above of the last cartesian block, used 'CARTESIAN COORDINATES' by default, string
    :param bottom_flag: string(s) below of the last cartesian block, used 'Empirical Formula:' by default, string
    :param pattern: regular expression for parsing lines in cartesian block, was determined by default, raw string
    :return: pandas DataFrame
    """
    top = tuple()
    bottom = tuple()
    rows = list()

    with open(path, 'r') as p:
        for num, line in enumerate(p, 1):
            if top_flag in line:
                top = num
            if bottom_flag in line:
                bottom = num

    with open(path, 'r') as p:
        for num, line in enumerate(p, 1):
            if (top <= num) and (bottom >= num):
                rows += re.findall(pattern, line)

    return pandas.DataFrame(rows, columns=['num', 'atom_name', 'x', 'y', 'z'])


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
