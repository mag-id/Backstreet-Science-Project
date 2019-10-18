
import pandas
import sys
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


if __name__ == '__main__':

    # directories
    template_path = str(sys.argv[1])
    workdir = str(sys.argv[2])

    # arguments for extract()
    mop_pattern = r'\s+([a-zA-Z]*)\s+(\-?\d+\.\d+)\s+\d\s+(\-?\d+\.\d+)\s+\d\s+(\-?\d+\.\d+)\s+\d\s+'
    mop_columns = 'atom_name', 'x', 'y', 'z'

    template = extract(template_path, pattern=mop_pattern, columns=mop_columns)

    for path in find(workdir, pattern=r'.+\.mop'):

        file = extract(path, pattern=mop_pattern, columns=mop_columns)

        # Initialize labels lists
        file_labels = list(file['atom_name'])
        template_labels = list(template['atom_name'])

        # Check and drop
        for string in hardcheck(template_labels, file_labels):
            file = file.drop(string)

        # Rewrite
        rewrite(file, path, pattern=mop_pattern)
