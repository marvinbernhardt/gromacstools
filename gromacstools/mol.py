import numpy as np


def bounding_box(mol):
    """returns the bounding box of a atom list as two verticies"""

    bb_vertex0 = np.array(
        [
            min([atom["r"][0] for atom in mol]),
            min([atom["r"][1] for atom in mol]),
            min([atom["r"][2] for atom in mol]),
        ]
    )
    bb_vertex1 = np.array(
        [
            max([atom["r"][0] for atom in mol]),
            max([atom["r"][1] for atom in mol]),
            max([atom["r"][2] for atom in mol]),
        ]
    )

    return bb_vertex0, bb_vertex1


def calc_com(mol):
    """returns the center of mass of a atom list"""
    com = np.zeros(3)
    mass = 0
    for atom in mol:
        com += atom["r"] * atom["mass"]
        mass += atom["mass"]
    com /= mass
    return com
