
"""molecule to graph"""

from gaocs import sdf_parser as pars
import networkx
import numpy

import matplotlib.pyplot as plt
import sys


def graph_from(dictionary: dict, array: numpy.array) -> networkx.Graph:
    mol_graph = networkx.Graph()

    for item in dictionary:
        mol_graph.add_node(item, label=dictionary[item])

    for item in array:
        mol_graph.add_edge(item[0], item[1], weight=item[2])

    return mol_graph


def graph_plot(graph: networkx.Graph) -> plt.show():
    networkx.draw(graph, with_labels=True)
    return plt.show()


if __name__ == '__main__':

    numpy.set_printoptions(threshold=sys.maxsize)

    ob_obj = pars.initialise('small_test.sdf')

    bonds = pars.parsed_bonds_from(ob_obj)
    atoms = pars.parsed_atoms_from(ob_obj)

    print('atoms: ', '\n', atoms, '\n')
    print('bonds: ', '\n', bonds, '\n')

    mol = graph_from(atoms, bonds)

    print('nodes: ', '\n', mol.nodes.data(), type(mol.nodes.data()), '\n')
    print('edges: ', '\n', mol.edges.data(), type(mol.edges.data()), '\n')

    print('cycles: ', '\n', networkx.cycle_basis(mol)), type(networkx.cycle_basis(mol))

    graph_plot(mol)
