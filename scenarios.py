
import sys

'''out-to-mop scenario'''

directory = str(sys.argv[1])

os.mkdir(directory + '/SP')
os.mkdir(directory + '/SP/_MOPs')
os.mkdir(directory + '/SP/_CSVs')

mop_dir = directory + '/SP/_MOPs/'
csv_dir = directory + '/SP/_CSVs/'

print(mop_dir, 'created')
print(csv_dir, 'created')

print('start')

for path in find(directory, pattern=r'test\d+\.out'):
    cartesian = extract(path)
    charge = charged(path)
    name = named(path)
    comment = commented(path)
    save(cartesian, charge, name, comment, mop_dir=mop_dir, csv_dir=csv_dir)
    print(path, 'is processed...')

print('done')


'''clean-and-rewrite scenario'''

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