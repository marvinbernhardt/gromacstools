import configparser
from itertools import combinations
from itertools import permutations

"""Note:
Theese tools do reading, generating and writing tables from and to itp files.
However: NO PARAMETERS ARE SUPPORTED. I wrote this mainly to deal with opls
itp files."""


def get_bonds(itp_file):
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes=[';'])
    config.optionxform = str
    config.read(itp_file)

    try:
        bonds_keys = config[' bonds '].keys()
    except KeyError:
        bonds_keys = []

    bonds = set()
    for line in bonds_keys:
        # bond is a set of {atom0, atom1}
        bond = frozenset(map(int, line.strip().split()[0:2]))
        bonds.add(bond)

    return bonds


def get_pairs(itp_file):
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes=[';'])
    config.optionxform = str
    config.read(itp_file)

    try:
        pairs_keys = config[' pairs '].keys()
    except KeyError:
        pairs_keys = []

    pairs = set()
    for line in pairs_keys:
        # pair is a set of {atom0, atom1}
        pair = frozenset(map(int, line.strip().split()[0:2]))
        pairs.add(pair)

    return pairs


def get_angles(itp_file):
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes=[';'])
    config.optionxform = str
    config.read(itp_file)

    try:
        angles_keys = config[' angles '].keys()
    except KeyError:
        angles_keys = []

    angles = set()
    for line in angles_keys:
        atoms = list(map(int, line.strip().split()[0:3]))
        # angle is an tuple of (middle_atom, set of {outer_atom1, outer_atom2})
        angle = (atoms[1], frozenset({atoms[0], atoms[2]}))
        angles.add(angle)

    return angles


def get_dihedrals(itp_file):
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes=[';'])
    config.optionxform = str
    config.read(itp_file)

    try:
        dihedrals_keys = config[' dihedrals '].keys()
    except KeyError:
        dihedrals_keys = []

    dihedrals = set()
    for line in dihedrals_keys:
        atoms = list(map(int, line.strip().split()[0:4]))
        # dihedral is a set of {tuple of (inner_atom, outer_atom), tuple of (inner_atom, outer_atom)}
        dihedral = frozenset({(atoms[1], atoms[0]), (atoms[2], atoms[3])})
        dihedrals.add(dihedral)

    return dihedrals


def pairs_from_bonds(bonds, distance=3):
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


def dihedrals_from_bonds(bonds):
    dihedrals = set()
    # iterate over all combinations of three bonds
    for bond_combo in permutations(bonds, 3):
        bond0 = bond_combo[0]
        bond1 = bond_combo[1]
        bond2 = bond_combo[2]
        common_atoms = bond0.intersection(bond1) | bond1.intersection(bond2) | bond2.intersection(bond0)
        assert len(common_atoms) < 3  # cycle of 3 atoms or double bond?

        # if they have two common atoms
        if len(common_atoms) == 2:
            outer_bonds = [bond for bond in [bond0, bond1, bond2] if bond != common_atoms]

            inner_atom0 = next(iter(common_atoms & outer_bonds[0]))
            inner_atom1 = next(iter(common_atoms & outer_bonds[1]))
            outer_atom0 = next(iter(outer_bonds[0] - common_atoms))
            outer_atom1 = next(iter(outer_bonds[1] - common_atoms))
            # dihedral is a set of {tuple of (inner_atom, outer_atom), tuple of (inner_atom, outer_atom)}
            dihedral = frozenset({(inner_atom0, outer_atom0), (inner_atom1, outer_atom1)})
            dihedrals.add(dihedral)

    return dihedrals


def angles_lines_from_angles(angles):

    angles_lines = []

    for angle in angles:
        inner_atom = angle[0]
        outer_atoms = list(angle[1])
        angle_line = f"{outer_atoms[0]:2d} {inner_atom:2d} {outer_atoms[1]:2d} 1"
        angles_lines.append(angle_line)

    angles_lines = '\n'.join(angles_lines)
    return(angles_lines)


def set_bonds(itp_file, bonds, bond_type='1'):
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes=['#', ';'])
    config.optionxform = str
    config.read(itp_file)

    config[' bonds '] = {}
    for bond in bonds:
        atoms = list(bond)
        bond_line = f"{atoms[0]:<2d} {atoms[1]:<2d}  {bond_type}"
        config[' bonds '][bond_line] = None

    with open(itp_file, 'w') as fp:
        config.write(fp)


def set_pairs(itp_file, pairs, pair_type='1'):
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes=['#', ';'])
    config.optionxform = str
    config.read(itp_file)

    config[' pairs '] = {}
    for pair in pairs:
        atoms = list(pair)
        pair_line = f"{atoms[0]:<2d} {atoms[1]:<2d}  {pair_type}"
        config[' pairs '][pair_line] = None

    with open(itp_file, 'w') as fp:
        config.write(fp)


def set_angles(itp_file, angles, angle_type='1'):
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes=['#', ';'])
    config.optionxform = str
    config.read(itp_file)

    config[' angles '] = {}
    for angle in angles:
        inner_atom = angle[0]
        outer_atoms = list(angle[1])
        angle_line = f"{outer_atoms[0]:<2d} {inner_atom:<2d} {outer_atoms[1]:<2d} {angle_type}"
        config[' angles '][angle_line] = None

    with open(itp_file, 'w') as fp:
        config.write(fp)


def set_dihedrals(itp_file, dihedrals, dihedral_type='3'):
    config = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes=['#', ';'])
    config.optionxform = str
    config.read(itp_file)

    config[' dihedrals '] = {}
    for dihedral in dihedrals:
        # dihedral is a set of {tuple of (inner_atom, outer_atom), tuple of (inner_atom, outer_atom)}
        dihedral_list = list(dihedral)
        atoms = [dihedral_list[0][1], dihedral_list[0][0], dihedral_list[1][0], dihedral_list[1][1]]
        dihedral_line = f"{atoms[0]:<2d} {atoms[1]:<2d} {atoms[2]:<2d} {atoms[3]:<2d} {dihedral_type}"
        config[' dihedrals '][dihedral_line] = None

    with open(itp_file, 'w') as fp:
        config.write(fp)
