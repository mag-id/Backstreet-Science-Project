
import klick
import projectlib as my

'''
import os
import sys
import projectlib as pl

#
# out-to-mop scenario
#

directory = str(sys.argv[1])

os.mkdir(directory + '/SP')
os.mkdir(directory + '/SP/_MOPs')
os.mkdir(directory + '/SP/_CSVs')

mop_dir = directory + '/SP/_MOPs/'
csv_dir = directory + '/SP/_CSVs/'

print(mop_dir, 'created')
print(csv_dir, 'created')

print('start')

for path in pl.find(directory, pattern=r'test\d+\.out'):
    cartesian = pl.extract(path)
    charge = pl.charged(path)
    name = pl.named(path)
    comment = pl.commented(path)
    pl.save(cartesian, charge, name, comment, mop_dir=mop_dir, csv_dir=csv_dir)
    print(path, 'is processed...')

print('done')


#
# clean-and-write scenario
#

# directories
template_path = str(sys.argv[1])
workdir = str(sys.argv[2])

# arguments for extract()
mop_pattern = r'\s+([a-zA-Z]*)\s+(\-?\d+\.\d+)\s+\d\s+(\-?\d+\.\d+)\s+\d\s+(\-?\d+\.\d+)\s+\d\s+'
mop_columns = 'atom_name', 'x', 'y', 'z'

template = pl.extract(template_path, pattern=mop_pattern, columns=mop_columns)

for path in pl.find(workdir, pattern=r'.+\.mop'):

    file = pl.extract(path, pattern=mop_pattern, columns=mop_columns)

    # Initialize labels lists
    file_labels = list(file['atom_name'])
    template_labels = list(template['atom_name'])

    # Check and drop
    for string in pl.hardcheck(template_labels, file_labels):
        file = file.drop(string)

    # Rewrite
    pl.rewrite(file, path, pattern=mop_pattern)


#
# initialization scenario
#

# deconstruct path to molecule
path, stem, suffix = pl.deconstruct(str(sys.argv[1]))

# initialize project directory
projectdir = pl.construct(stem='test_project')

# construct and create molecule folder
moldir = pl.construct(projectdir, stem)
pl.create(moldir)

# convert molecule into molecule folder
pl.convert('mopout', pl.construct(path, stem, suffix),
           'sdf', pl.construct(moldir, stem, '.sdf'))

# construct and create confsearch/genetic_default folder
confdir = pl.construct(moldir, 'confsearch/genetic_default')
pl.create(confdir)

# ...

# slice multiple-structure file to singe-structure file
pl.confslice('sdf', pl.construct(confdir, stem+'_conformers', '.sdf'),
          'sdf', pl.construct(confdir, stem+'_conformer', '.sdf'))
'''