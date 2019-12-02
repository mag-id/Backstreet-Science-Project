
# TODO: Write tests for Structure class instances with internaltests for sdf_parser
# TODO: Write tests for sdf_parser

from gaocs import structure_class as my
from gaocs import sdf_parser as pars
import unittest
import numpy


class StructureTest(unittest.TestCase):

    def test__add__(self):
        mol_one = mol + iso_one
        self.assertEqual(mol_one.show(),
                         [['mol', '0 H [-1. -1.  0.]', '1 O [0. 0. 0.]', '2 H [1. 1. 0.]'], ['None']],
                         'test__add__ for mol_one.show()')

        self.assertEqual(mol.show(),
                         [['mol', '0 H [-1. -1.  0.]', '1 O [0. 0. 0.]', '2 H [ 1. -1.  0.]'], ['None']],
                         'test__add__ for mol.show()')

    def test__sub__(self):
        mol_one = mol + iso_one
        mol_new_one = mol_one - iso_one
        self.assertEqual(mol_new_one.show(),
                         [['mol', '0 H [-1. -1.  0.]', '1 O [0. 0. 0.]', '2 H [ 1. -1.  0.]'], ['None']],
                         'test__add__ for mol_one.show()')

        self.assertEqual(mol_new_one is mol, False, 'mol_new_one is mol')


mol = my.Structure(
    numpy.array(
        [[-1.0, -1.0, 0.0],
         [0.0, 0.0, 0.0],
         [1.0, -1.0, 0.0]]),
    {0: 'H', 1: 'O', 2: 'H'},
    'mol')

iso_one = my.Structure(
    numpy.array(
        [[0.0, 0.0, 0.0],
         [0.0, 2.0, 0.0]]),
    {0: 'H', 2: 'H'},
    'iso_one')

iso_two = my.Structure(
    numpy.array(
        [[0.0, 0.0, 0.0],
         [0.0, 1.0, 0.0]]),
    {0: 'H', 2: 'H'},
    'iso_two')


if __name__ == '__main__':
    unittest.main()
