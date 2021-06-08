import re
from itertools import combinations
from itertools import permutations

"""
Note:
Theese tools do reading, generating and writing tables from and to itp files.
However: NO PARAMETERS ARE SUPPORTED. I wrote this mainly to deal with opls
itp files.
"""


def get_section_lines(filepath, section):
    """
    Parse itp file for lines of a section

    Parameters
    ----------
    filepath : str
        Filepath for file to be parsed
    section : str
        section to parse for

    Returns
    -------
    lines : list
        list of lines

    """
    _reg_section = re.compile(r'\s*\[\s*([a-zA-Z0-9-_]+)\s*\]\s*')
    _reg_comment = re.compile(r'\s*;')
    _reg_statement = re.compile(r'\s*#')
    _reg_empty = re.compile(r'\s*' + '\n')

    section_lines = []
    with open(filepath, 'r') as f:
        lines = f.readlines()

        cur_section = None
        for line in lines:

            # if comment or empty
            if (_reg_comment.match(line) or _reg_empty.match(line)
                    or _reg_statement.match(line)):
                continue

            # if section header
            match_section = _reg_section.match(line)
            if match_section:
                cur_section = match_section.group(1)
                continue

            # if in wanted section
            if cur_section == section:
                section_lines.append(line)
    return section_lines


def get_bonds(filepath):
    """
    Parse itp file for bonds

    Parameters
    ----------
    filepath : str
        Filepath for file to be parsed

    Returns
    -------
    bonds : set

    """
    _reg_bond = re.compile(r'\s*([0-9]+)\s*([0-9]+)\s*([0-9]+)')

    bonds_lines = get_section_lines(filepath, 'bonds')

    bonds = set()
    for line in bonds_lines:
        match_bond = _reg_bond.match(line)
        if match_bond is None:
            raise Exception('unable to parse line into bond:', line)
        ai, aj, func = map(int, match_bond.group(1, 2, 3))
        bond = frozenset({ai, aj})
        bonds.add(bond)
    return bonds


def get_pairs(filepath):
    """
    Parse itp file for pairs

    Parameters
    ----------
    filepath : str
        Filepath for file to be parsed

    Returns
    -------
    pairs : set

    """
    _reg_pair = re.compile(r'\s*([0-9]+)\s*([0-9]+)\s*([0-9]+)')

    pairs_lines = get_section_lines(filepath, 'pairs')

    pairs = set()
    for line in pairs_lines:
        match_pair = _reg_pair.match(line)
        if match_pair is None:
            raise Exception('unable to parse line into pair:', line)
        ai, aj, func = map(int, match_pair.group(1, 2, 3))
        pair = frozenset({ai, aj})
        pairs.add(pair)
    return pairs


def get_angles(filepath):
    """
    Parse itp file for angles

    Parameters
    ----------
    filepath : str
        Filepath for file to be parsed

    Returns
    -------
    angles : set

    """
    _reg_angle = re.compile(r'\s*([0-9]+)\s*([0-9]+)\s*([0-9]+)\s*([0-9]+)')

    angles_lines = get_section_lines(filepath, 'angles')

    angles = set()
    for line in angles_lines:
        match_angle = _reg_angle.match(line)
        if match_angle is None:
            raise Exception('unable to parse line into angle:', line)
        # angle is an tuple of (middle_atom, set of {outer_atom1, outer_atom2})
        ai, aj, ak, func = map(int, match_angle.group(1, 2, 3, 4))
        angle = (aj, frozenset({ai, ak}))
        angles.add(angle)
    return angles


def get_dihedrals(filepath, ignore_funcs=(2, 4)):
    """
    Parse itp file for dihedrals

    Parameters
    ----------
    filepath : str
        Filepath for file to be parsed
    ignore_funcs : tuple
        Dihedral funcs to be ignored. Default is (2, 4) which are improper
        dihedrals in Gromacs.

    Returns
    -------
    dihedrals : set

    """
    _reg_dihedral = re.compile(r'\s*([0-9]+)\s*([0-9]+)'
                               + r'\s*([0-9]+)\s*([0-9]+)\s*([0-9]+)')

    dihedrals_lines = get_section_lines(filepath, 'dihedrals')

    dihedrals = set()
    for line in dihedrals_lines:
        match_dihedral = _reg_dihedral.match(line)
        if match_dihedral is None:
            raise Exception('unable to parse line into dihedral:', line)
        ai, aj, ak, al, func = map(int, match_dihedral.group(1, 2, 3, 4, 5))
        # exclude type 2 and 4 which are improper dihedrals
        if func in ignore_funcs:
            continue
        # dihedral is a set of {tuple of (inner_atom, outer_atom), tuple of
        # (inner_atom, outer_atom)}
        dihedral = frozenset({(aj, ai), (ak, al)})
        dihedrals.add(dihedral)
    return dihedrals


def dihedrals_from_bonds(bonds):
    """
    Generate dihedrals from bonds

    Parameters
    ----------
    bonds : set
        bonds describing a molecule

    Returns
    -------
    dihedrals : set

    """
    dihedrals = set()
    # iterate over all combinations of three bonds
    for bond_combo in permutations(bonds, 3):
        bond0 = bond_combo[0]
        bond1 = bond_combo[1]
        bond2 = bond_combo[2]
        common_atoms = (
            bond0.intersection(bond1)
            | bond1.intersection(bond2)
            | bond2.intersection(bond0)
        )
        if len(common_atoms) == 3:  # cycle of 3 atoms
            continue

        # if they have two common atoms
        if len(common_atoms) == 2:
            outer_bonds = [
                bond for bond in [bond0, bond1, bond2] if bond != common_atoms
            ]

            inner_atom0 = next(iter(common_atoms & outer_bonds[0]))
            inner_atom1 = next(iter(common_atoms & outer_bonds[1]))
            outer_atom0 = next(iter(outer_bonds[0] - common_atoms))
            outer_atom1 = next(iter(outer_bonds[1] - common_atoms))
            # dihedral is a set of {tuple of (inner_atom, outer_atom), tuple of
            # (inner_atom, outer_atom)}
            dihedral = frozenset(
                {(inner_atom0, outer_atom0), (inner_atom1, outer_atom1)}
            )
            dihedrals.add(dihedral)

    return dihedrals


def pairs_from_bonds(bonds, distance=3):
    """
    Generate pairs from bonds

    Parameters
    ----------
    bonds : set
        bonds describing a molecule
    distance : int
        distance (in bonds) the pairs should have

    Returns
    -------
    pairs : set

    """
    pairs = set()
    if distance == 3:
        # iterate over all combinations of two bonds
        for dihedral in dihedrals_from_bonds(bonds):
            dihedral_iter = iter(dihedral)
            pair = frozenset({next(dihedral_iter)[1], next(dihedral_iter)[1]})
            pairs.add(pair)
    else:
        raise Exception(f"Pair list for distance of {distance} not implemented yet")

    return pairs


def angles_from_bonds(bonds):
    """
    Generate angles from bonds

    Parameters
    ----------
    bonds : set
        bonds describing a molecule

    Returns
    -------
    angles : set

    """
    angles = set()
    # iterate over all combinations of two bonds
    for bond_combo in combinations(bonds, 2):
        bond0 = bond_combo[0]
        bond1 = bond_combo[1]
        common_atom = bond0.intersection(bond1)
        assert len(common_atom) < 2  # double bond?

        # if they have common atom
        if len(common_atom) == 1:
            angle = (next(iter(common_atom)), bond0 - common_atom | bond1 - common_atom)
            angles.add(angle)

    return angles


def angles_lines_from_angles(angles, func=1):
    """
    Generate angles lines for a topology file

    Parameters
    ----------
    angles : set

    Returns
    -------
    angles_lines : str
    """

    angles_lines = []

    for angle in angles:
        inner_atom = angle[0]
        outer_atoms = list(angle[1])
        angle_line = f"{outer_atoms[0]:2d} {inner_atom:2d} {outer_atoms[1]:2d}  {func}"
        angles_lines.append(angle_line)

    angles_lines = "\n".join(angles_lines)
    return angles_lines
