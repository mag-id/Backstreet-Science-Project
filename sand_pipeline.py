
import filesys as fs
import pybellogic as pb
from confslice import confslice


# deconstruct path to molecule
path, stem, suffix = fs.deconstruct('/Users/mag/Code/pycharm/work/test_test1/test1.out')

# initialize project directory
projectdir = fs.construct(stem='test_project')

# construct and create molecule folder
moldir = fs.construct(projectdir, stem)
fs.create(moldir)

# convert molecule into molecule folder
pb.convert('mopout', fs.construct(path, stem, suffix),
           'sdf', fs.construct(moldir, stem, '.sdf'))

# construct and create confsearch/genetic_default folder
confdir = fs.construct(moldir, 'confsearch/genetic_default')
fs.create(confdir)

# default genetic conformer search
pb.geneticsearch(fs.construct(moldir, stem, '.sdf'),
                 fs.construct(confdir, stem+'_conformers', '.sdf'))

# slice multiple-structure file to singe-structure file
confslice('sdf', fs.construct(confdir, stem+'_conformers', '.sdf'),
          'sdf', fs.construct(confdir, stem+'_conformer', '.sdf'))
