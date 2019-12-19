
import pandas
import networkx

import openbabel
import pybel

from alt_gaocs import variables as var


class Structure:
    def __init__(self, name, coordinates, graph) -> None:

        self.name: str = name
        self.coordinates: pandas.DataFrame = coordinates
        self.graph: networkx.Graph = graph

    def write(self, path_to_write: str) -> None:

        openbabel_obj = openbabel.OBMol()

        for index, row in self.coordinates.iterrows():
            atom = openbabel_obj.NewAtom()
            atom.SetAtomicNum(int(row['atomic_num']))
            atom.SetVector(float(row['x']), float(row['y']), float(row['z']))

        pybel_obj = pybel.Molecule(openbabel_obj)
        pybel_obj.write(var._write, path_to_write)

    @classmethod
    def read(cls, path_to_sdf: str) -> object:

        openbabel_obj = cls.__initialise(path_to_sdf)

        name = cls.__get_name(openbabel_obj)
        coordinates = cls.__get_coordinates(openbabel_obj)
        graph = cls.__get_graph(openbabel_obj)

        return cls(name, coordinates, graph)

    @staticmethod
    def __initialise(path_to_sdf: str) -> openbabel.OBMol:

        obconversion = openbabel.OBConversion()
        obconversion.SetInFormat(var._initialize)

        openbabel_obj = openbabel.OBMol()
        obconversion.ReadFile(openbabel_obj, path_to_sdf)

        return openbabel_obj

    @staticmethod
    def __get_name(openbabel_obj: openbabel.OBMol) -> str:
        pybel_obj = pybel.Molecule(openbabel_obj)
        return str(pybel_obj.title)

    @staticmethod
    def __get_coordinates(openbabel_obj: openbabel.OBMol) -> pandas.DataFrame:

        atom_coords = list()
        atom_types = list()

        pybel_obj = pybel.Molecule(openbabel_obj)

        for atom in pybel_obj:
            atom_coords.append(atom.coords)
            atom_types.append(atom.atomicnum)

        dataframe = pandas.DataFrame(atom_coords, columns=['x', 'y', 'z'], dtype=var._cooridates)
        dataframe.insert(0, 'atomic_num', atom_types)
        dataframe['atomic_num'] = dataframe['atomic_num'].astype(var._atomic_num)

        return dataframe

    @staticmethod
    def __get_graph(openbabel_obj: openbabel.OBMol) -> networkx.Graph:

        openbabel_iter_obj = openbabel.OBMolBondIter(openbabel_obj)
        molecular_graph = networkx.Graph()

        for bond in openbabel_iter_obj:
            molecular_graph.add_edge(bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1, weight=bond.GetBondOrder())

        return molecular_graph
