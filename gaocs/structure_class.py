
import numpy
import copy


# TODO: write docs!
class Structure:
    """Structure class"""
    def __init__(self, cartesian: numpy.array, atoms: dict, name: str, bonds: numpy.array = None) -> None:
        """
        :param cartesian: the 2D numpy.array of the float numbers with size 3 columns and n rows
        :param atoms: the dictionary with indexes as keys and atom labels as values
        :param name: the string with the name of a structure
        :param bonds: the 2D numpy.array of the integer numbers with size 3 columns and m rows, None by default
        :raise IndexError if length of the 'cartesian' parameter is not equal to length of the 'atoms' parameter
        """
        self.cartesian: numpy.array = cartesian
        self.atoms: dict = atoms
        self.name: str = name
        self.bonds: numpy.array = bonds

        if len(self.cartesian) != len(self.atoms):
            raise IndexError(f'len(cartesian) != len(atoms) at {self.name}')

    def __add__(self, addend) -> object:
        """
        :param addend: Structure class instance
        :raise IndexError if addend.atoms is not subset of the self.atoms
        :return: new instance of the Structure class
        """
        cartesian = copy.deepcopy(self.cartesian)
        atoms = copy.deepcopy(self.atoms)
        name = copy.deepcopy(self.name)
        bonds = copy.deepcopy(self.bonds)

        for row_index, label_index in enumerate(addend.atoms):
            if addend.atoms[label_index] == atoms[label_index]:
                cartesian[label_index] += addend.cartesian[row_index]
            else:
                raise IndexError(f'{addend.name}.atoms is not subset of the {name}.atoms')
        return Structure(cartesian, atoms, name, bonds)

    def __sub__(self, addend) -> object:
        """
        :param addend: addend: Structure class instance
        :raise IndexError if addend.atoms is not subset of the self.atoms
        :return: new instance of the Structure class
        """
        cartesian = copy.deepcopy(self.cartesian)
        atoms = copy.deepcopy(self.atoms)
        name = copy.deepcopy(self.name)
        bonds = copy.deepcopy(self.bonds)

        for row_index, label_index in enumerate(addend.atoms):
            if addend.atoms[label_index] == atoms[label_index]:
                cartesian[label_index] -= addend.cartesian[row_index]
            else:
                raise IndexError(f'{addend.name}.atoms is not subset of the {name}.atoms')
        return Structure(cartesian, atoms, name, bonds)

    def show(self) -> list:
        """
        Show the new which stored into instance of the Structure class
        :return: list
        """
        lines = [[self.name], []]

        for row_index, label_index in enumerate(self.atoms):
            lines[0].append(f'{label_index} {self.atoms[label_index]} {self.cartesian[row_index]}')

        if self.bonds is not None:
            lines[1].append(self.bonds)
        else:
            lines[1].append('None')

        return lines
