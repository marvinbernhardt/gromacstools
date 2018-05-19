import fileinput
import numpy as np
import sys
from .general import run_bash


def get_natoms(gro_file):
    """returns the number of atoms from a .gro file"""
    with open(gro_file, 'r') as file:
        for i, line in enumerate(file):
            if i == 1:
                natoms = int(line.strip())
                break
    return natoms


def get_box(gro_file):
    """returns the box vector from a .gro file"""
    natoms = get_natoms(gro_file)
    with open(gro_file, 'r') as file:
        for line_nr, line in enumerate(file):
            if line_nr == natoms + 2:
                box = list(map(float, line.split()))
    return np.array(box)


def set_box(gro_file, box_length):
    natoms = get_natoms(gro_file)
    for line_nr, line in enumerate(fileinput.input(gro_file, inplace=True)):
        if line_nr == natoms + 2:
            line = f"{box_length} {box_length} {box_length}\n"
        sys.stdout.write(line)


def mix_water(gro_file, nmols):
    """perutes the indices of water molecules in a .gro file.
Does only work on files, where first nmols molecules are water with three atoms
"""
    with open(gro_file, 'r') as f:
        content = f.readlines()
    header = content[0:2]
    footer = content[-1:]
    coordlines = content[2:-1]
    mols = np.reshape(coordlines[:nmols*3], (nmols, 3))
    mols = np.random.permutation(mols)
    with open(gro_file, 'w') as f:
        for line in header:
            f.write(line)
        for i, mol in enumerate(mols):
            for j, atomline in enumerate(mol):
                atomline = "{:>5}{}{:>5}{}".format(i + 1, atomline[5:15], 3 * i + j + 1, atomline[20:])
                f.write(atomline)
        for line in footer:
            f.write(line)

def read_moltypes(gro_file, atom_mass_dict, abc_indicators_dict, sigma_dict):
    """read gro file into moltypes"""
    # read atoms from file
    atoms = []
    with open(gro_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                pass
            elif i == 1:
                natoms = int(line.strip())
            elif i < natoms + 2:
                atom = {}
                atom['mol_nr'] = int(line[0:5])
                atom['mol_name'] = line[5:10].strip()
                atom['name'] = line[10:15].strip()
                #atom['atom_nr'] = int(line[15:20])
                atom['r'] = np.array([float(line[20:28]), float(line[28:36]), float(line[36:44])])
                # try to read velocities if present in file
                try:
                    atom['v'] = np.array([float(line[44:52]), float(line[52:60]), float(line[60:68])])
                except:
                    atom['v'] = np.full(3, np.nan)

                atoms.append(atom)

    # add atoms mass from dict
    for atom in atoms:
        atom['mass'] = atom_mass_dict[atom['name']]

    # create mols list
    mols = []
    mol_nr_old = 1
    mol = []
    for atom in atoms:
        if atom['mol_nr'] != mol_nr_old:
            mols.append(mol)
            mol = []
            mol_nr_old = atom['mol_nr']
        del atom['mol_nr']
        mol.append(atom)
    mols.append(mol)

    # create moltypes list
    moltypes = []
    mol_name_old = ""
    for mol in mols:
        # check for new moltype
        mol_name = mol[0]['mol_name']
        if mol_name_old != mol_name:
            moltype = {'name': mol_name, 'mols': [], 'abc_indicators': abc_indicators_dict[mol_name],
                       'sigma': sigma_dict[mol_name]}
            moltypes.append(moltype)
            mol_name_old = mol_name

        # count mols of that type
        if mol[0]['mol_name'] == mol_name_old:
            moltype['mols'].append(mol)
        for atom in mol:
            del atom['mol_name']

    return moltypes

def set_atomname(gro_file, atomlist, atomname):
    """modifies atoms in atomlist to have name atomname"""
    if len(atomname) > 5:
            raise Exception("atomname has {} characters, but can only have 5".format(len(atomname)))
    for i, line in enumerate(fileinput.input(gro_file, inplace=True)):
        if i - 2 in atomlist:
            line = line[0:10] + "{:>5}".format(atomname) + line[15:]
        # sys.stdout is redirected to the file
        sys.stdout.write(line)


def set_coordinate(gro_file, atom, coord, value):
    """modifies the coordinate of atom to have value"""
    for i, line in enumerate(fileinput.input(gro_file, inplace=True)):
        if i == atom + 2:
            first_char = 20 + 8 * coord
            last_char = first_char + 8
            line = line[0:first_char] + "{:8.3f}".format(value) + line[last_char:]
        # sys.stdout is redirected to the file
        sys.stdout.write(line)


def set_molname(gro_file, atomlist, molname):
    """set molecule name of atoms in atomlist to be molname"""
    if len(molname) > 5:
            raise Exception("molname has {} characters, but can only have 5".format(len(molname)))
    for i, line in enumerate(fileinput.input(gro_file, inplace=True)):
        if i - 2 in atomlist:
            line = line[0:5] + "{:<5}".format(molname) + line[10:]
        # sys.stdout is redirected to the file
        sys.stdout.write(line)


def translate_with_pbc(gro_file, vector, out_file="conf-shifted.gro"):
    """tanslates the molecules in a gro_file along vector and uses pbc to put them all into the box"""

    x, y, z = vector[0], vector[1], vector[2]
    run_bash(f"gmx editconf -f {gro_file} -o conf-temp-trans.gro -translate {x} {y} {z}")
    run_bash("gmx trjconv -f conf-temp-trans.gro -o {out_file} -pbc mol <<< 0")
    run_bash("rm conf-temp-trans.gro")


def generate_ffc_crystal_atomic(filename, n_atoms, box_length, mol_name, atom_name):
    """generates a face-centerd-cubic crystal of single atoms.

Works only if natoms/4 is a cubic number.

Parameters
----------
filename : string
    Name of output file.
n_atoms: int
    Number of atoms.
box_length: scalar
    Length of the cubic box
mol_name: string
    Name of moltype (max. length 5)
atom_name: string
    Name of atoms (max. length 5)
"""

    def is_perfect_cube(number):
        return number in (x**3 for x in range(20))

    if not n_atoms % 4 == 0:
        raise ValueError("number of atoms not divisible by 4 (number of points in fcc unit cell)")
    if not is_perfect_cube(n_atoms // 4):
        raise ValueError(f"number of atoms divided by 4 ({n_atoms // 4}) is not a cubic number")

    lattice_points=[(0., 0., 0.), (0., .5, .5), (.5, 0., .5), (.5, .5, 0.)]
    nx = round((n_atoms // 4)**(1/3))
    a = b = c = box_length / nx

    with open(filename, 'w') as f:
        atom_nr = 1
        f.write("argon\n")
        f.write("  {}\n".format(n_atoms))
        for x in np.linspace(0., box_length, num=nx, endpoint=False):
            for y in np.linspace(0., box_length, num=nx, endpoint=False):
                for z in np.linspace(0., box_length, num=nx, endpoint=False):
                    for lattice_point in lattice_points:
                        xpos = x + lattice_point[0] * a
                        ypos = y + lattice_point[1] * b
                        zpos = z + lattice_point[2] * c
                        f.write("{0: 5d}{1:>5s}{2:>5s}{0: 5d}{3:8.3f}{4:8.3f}{5:8.3f}\n".format(atom_nr, mol_name, atom_name, xpos, ypos, zpos))
                        atom_nr += 1

        f.write(" {0:10.5f}{0:10.5f}{0:10.5f}\n".format(box_length))


def generate_ffc_crystal_water(filename, n_mols, box_length):
    """generates a face-centerd-cubic crystal of SPC water.

All molecules are orientated in the same direction. Works only
if nmols/4 is a cubic number.

Parameters
----------
filename : string
    Name of output file.
n_mols: int
    Number of molecules.
box_length: scalar
    Length of the cubic box
"""

    def is_perfect_cube(number):
        return number in (x**3 for x in range(20))

    n_atoms = n_mols * 3
    if not n_mols % 4 == 0:
        raise ValueError("number of molecules not divisible by 4 (number of points in fcc unit cell)")
    if not is_perfect_cube(n_mols // 4):
        raise ValueError(f"number of molecules divided by 4 ({n_mols // 4}) is not a cubic number")

    nx = round((n_mols // 4)**(1/3))

    atoms = [("OW", [ 0., 0., 0.00646032]), ("HW1", [ 0.08164904, 0., -0.05127557]), ("HW2", [-0.08164904, 0., -0.05127557])]

    lattice_points=[(0., 0., 0.), (0., .5, .5), (.5, 0., .5), (.5, .5, 0.)]
    a = b = c = box_length / nx

    with open(filename, 'w') as f:
        atom_nr = 1
        mol_nr = 1
        f.write("water\n")
        f.write("  {}\n".format(n_atoms))
        for x in np.linspace(0., box_length, num=nx, endpoint=False):
            for y in np.linspace(0., box_length, num=nx, endpoint=False):
                for z in np.linspace(0., box_length, num=nx, endpoint=False):
                    for lattice_point in lattice_points:
                        for atom in atoms:
                            atom_name = atom[0]
                            atom_pos = atom[1]
                            xpos = x + lattice_point[0] * a + atom_pos[0]
                            ypos = y + lattice_point[1] * b + atom_pos[1]
                            zpos = z + lattice_point[2] * c + atom_pos[2]
                            f.write("{0: 5d}SOL    {1:3s}{2: 5d}{3:8.3f}{4:8.3f}{5:8.3f}\n"\
                                    .format(mol_nr, atom_name, atom_nr, xpos, ypos, zpos))
                            atom_nr += 1
                        mol_nr += 1

        f.write(" {0:10.5f}{0:10.5f}{0:10.5f}\n".format(box_length))
