
import numpy
import copy


# TODO: write docs!
class Structure:
    """Structure class"""
    def __init__(self, cartesian: numpy.array, atoms: dict, name: str) -> None:
        """"""
        self.cartesian: numpy.array = cartesian
        self.atoms: dict = atoms
        self.name: str = name

        if len(self.cartesian) != len(self.atoms):
            raise IndexError(f'len(cartesian) != len(atoms) at {self.name}')

    def __add__(self, addend) -> object:
        """"""
        for row_index, label_index in enumerate(addend.atoms):
            if addend.atoms[label_index] == self.atoms[label_index]:
                self.cartesian[label_index] += addend.cartesian[row_index]
            else:
                raise IndexError(f'{self.name} and {addend.name} atoms are not complementary')
        return Structure(self.cartesian, self.atoms, self.name)

    def __sub__(self, addend) -> object:
        """"""
        for row_index, label_index in enumerate(addend.atoms):
            if addend.atoms[label_index] == self.atoms[label_index]:
                self.cartesian[label_index] -= addend.cartesian[row_index]
            else:
                raise IndexError(f'{self.name} and {addend.name} atoms are not complementary')
        return Structure(self.cartesian, self.atoms, self.name)

    def show(self) -> list:
        """"""
        lines = list([self.name])
        for row_index, label_index in enumerate(self.atoms):
            lines.append(f'{label_index} {self.atoms[label_index]} {self.cartesian[row_index]}')
        return lines


# For future tests
"""
mol = Structure(
    numpy.array(
        [[-1.0, -1.0, 0.0],
         [0.0, 0.0, 0.0],
         [1.0, -1.0, 0.0]]),
    {0: 'H', 1: 'O', 2: 'H'},
    'mol')

iso_one = Structure(
    numpy.array(
        [[0.0, 0.0, 0.0],
         [0.0, 2.0, 0.0]]),
    {0: 'H', 2: 'H'},
    'iso_one')

iso_two = Structure(
    numpy.array(
        [[0.0, 0.0, 0.0],
         [0.0, 1.0, 0.0]]),
    {0: 'H', 2: 'H'},
    'iso_two')

a = 1
b = 1
c = a + b
print(c, a)

structures = [mol, iso_one, iso_two]
for structure in structures:
    print(structure.show())

mol2_1 = mol + iso_one
mol2_2 = mol + iso_two

print(mol2_1.show())
print(mol2_2.show())
"""
