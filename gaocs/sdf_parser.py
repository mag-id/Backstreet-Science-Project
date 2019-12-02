
"""sdf parser"""

# TODO: Write docs!

from gaocs import structure_class as my
import openbabel
import pybel
import numpy


def pick(molecule: str) -> my.Structure:

    obj = initialise(molecule)

    cartesian = parsed_cartesian_from(obj)
    atoms = parsed_atoms_from(obj)
    name = parsed_name_from(obj)
    bonds = parsed_bonds_from(obj)

    return my.Structure(cartesian, atoms, name, bonds)


def initialise(sdf: str) -> openbabel.OBMol:
    """
    :param sdf: path to the file.sdf, str
    :return: OpenBabel.OBMol object
    """

    obconversion = openbabel.OBConversion()
    obconversion.SetInFormat('sdf')

    obj = openbabel.OBMol()
    obconversion.ReadFile(obj, sdf)

    return obj


def parsed_cartesian_from(obj: openbabel.OBMol) -> numpy.array:
    xyz = list()
    molecule = pybel.Molecule(obj)

    for atom in molecule:
        xyz.append(atom.coords)

    return numpy.array(xyz, dtype='float16')


def parsed_atoms_from(obj: openbabel.OBMol) -> dict:
    dictionary = dict()
    molecule = pybel.Molecule(obj)

    for atom in molecule:
        dictionary.update({atom.idx - 1: atom.type})

    return dictionary


def parsed_name_from(obj: openbabel.OBMol) -> str:
    molecule = pybel.Molecule(obj)
    return str(molecule.title)


def parsed_bonds_from(obj: openbabel.OBMol) -> numpy.array:
    bonds = list()
    molecule = openbabel.OBMolBondIter(obj)

    for bond in molecule:
        bonds.append((bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1, bond.GetBondOrder()))

    return numpy.array(bonds, dtype='int16')


if __name__ == '__main__':
    for line in (pick('small_test.sdf').show()):
        print(line)
