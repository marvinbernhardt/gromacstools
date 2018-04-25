from copy import deepcopy
import numpy as np
from scipy import linalg

class Atom:
    """Atom class. Can contain positions and velocities.
    Modifications to objects ['name', 'mass', 'pos', 'vel'] will also modify
    the distinctive topology (__setattr__ method is altered)."""

    def __init__(self, dt_atom, index_in_top, index_in_moltype, index_in_mol):
        self.__dt_atom = dt_atom
        self.__dict__['name'] = dt_atom['name']
        self.__dict__['mass'] = dt_atom['mass']
        self.index_in_top = index_in_top
        self.index_in_moltype = index_in_moltype
        self.index_in_mol = index_in_mol
        try:
            self.__dict__['pos'] = dt_atom['pos']
        except:
            pass
        try:
            self.__dict__['vel'] = dt_atom['vel']
        except:
            pass

    def __setattr__(self, attribute, value):
        """additionally modify dt_atom in the distinctive topology"""
        if attribute in ['name', 'mass', 'pos', 'vel']:
            self.__dt_atom[attribute] = value
        self.__dict__[attribute] = value


class Molecule:
    """Molecule class. Contains atoms() function.
    It is not possible to chainge it's name (that comes from the moltype)
    nor its atoms (not implemented)."""

    def __init__(self, dt_mol, name, index_in_top, index_in_moltype,
                 firstatom, moltype_nr):
        self.name = name
        self.__dt_atoms = dt_mol['atoms']
        self.mass = sum([dt_atom['mass'] for dt_atom in self.__dt_atoms])
        self.natoms = len(self.__dt_atoms)
        self.index_in_top = index_in_top
        self.index_in_moltype = index_in_moltype
        self.firstatom = firstatom
        self.moltype_nr = moltype_nr

    def atoms(self):
        """Return list of Atom objects in this Molecule."""
        atoms = []
        for index_in_mol, dt_atom in enumerate(self.__dt_atoms):
            index_in_moltype = (self.index_in_moltype * self.natoms
                                + index_in_mol)
            index_in_top = self.firstatom + index_in_mol
            atoms.append(Atom(dt_atom, index_in_top, index_in_moltype,
                              index_in_mol))
        return atoms


class Moltype:
    """Moltype class. Contains mols() and atoms() function.
    Modifications to objects ['name', 'rot_treat', 'abc_indicators', 'sigma'] will also modify
    the distinctive topology (__setattr__ method is altered)."""

    def __init__(self, dt_moltype, index_in_top, firstmol, firstatom):
        self.__dt_moltype = dt_moltype
        self.__dict__['name'] = dt_moltype['name']
        self.dt_mols = dt_moltype['mols']
        self.__dict__['rot_treat'] = dt_moltype['rot_treat']
        self.__dict__['abc_indicators'] = dt_moltype['abc_indicators']
        self.__dict__['sigma'] = dt_moltype['sigma']
        self.nmols = len(self.dt_mols)
        self.natomtypes = len(self.atomtypes())
        self.natoms = self.nmols * self.natomtypes
        self.mass = self.nmols * sum([atomtype['mass'] for atomtype
                                      in self.atomtypes()])
        self.index_in_top = index_in_top
        self.firstmol = firstmol
        self.firstatom = firstatom

    def __setattr__(self, attribute, value):
        """additionally modify dt_moltype in the distinctive topology"""
        if attribute in ['name', 'rot_treat', 'abc_indicators', 'sigma']:
            self.__dt_moltype[attribute] = value
        self.__dict__[attribute] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def mols(self):
        """Return list of Molecule objects in this molecule type."""
        mols = []
        firstatom = self.firstatom
        for index_in_moltype, dt_mol in enumerate(self.dt_mols):
            index_in_top = self.firstmol + index_in_moltype
            mols.append(Molecule(dt_mol, self.name, index_in_top,
                                 index_in_moltype,
                                 firstatom, self.index_in_top))
            firstatom += len(self.atomtypes())
        return mols

    def atoms(self):
        """Return list of Atom objects in this molecule type."""
        atoms = []
        for mol in self.mols():
            for atom in mol.atoms():
                atoms.append(atom)
        return atoms

    def atomtypes(self):
        """Return list of atomtypes. Needed internally."""
        atomtypes = []
        for dt_atom in self.dt_mols[0]['atoms']:
            atomtype = {'name': dt_atom['name'], 'mass': dt_atom['mass']}
            atomtypes.append(atomtype)
        return atomtypes


class Topology:
    """Topology Class for MD simulations

    Uses a distinctive topology based on dicts and lists insinde.
    Gives lists of objects to the outside.

    Contains moltypes(), mols() and atoms() function."""

    def moltypes(self):
        moltypes = []
        firstmol = 0
        firstatom = 0
        for index_in_top, dt_moltype in enumerate(self.distinctive_top):
            moltype = Moltype(dt_moltype, index_in_top, firstmol, firstatom)
            moltypes.append(moltype)
            firstmol += moltype.nmols
            firstatom += moltype.natoms
        return moltypes

    def nmoltypes(self):
        return len(self.moltypes())

    def nmols(self):
        return sum([moltype.nmols for moltype in self.moltypes()])

    def natoms(self):
        return sum([moltype.natoms for moltype in self.moltypes()])

    def mass(self):
        return sum([moltype.mass for moltype in self.moltypes()])

    def mols(self):
        mols = []
        for moltype in self.moltypes():
            for mol in moltype.mols():
                mols.append(mol)

        return mols

    def atoms(self):
        atoms = []
        for mol in self.mols():
            for atom in mol.atoms():
                atoms.append(atom)

        return atoms

    def load_top_file(self, top_filename):
        """not implemented yes"""
        pass

    def load_simple_top(self, simple_top):
        """Loads a simple topology (moltypes with nmols) and converts
it to a distinctive topology."""
        self.distinctive_top = deepcopy(simple_top)

        # expand nmols to mols
        for dt_moltype in self.distinctive_top:
            mol = {'atoms': [atom.copy() for atom in dt_moltype['atoms']]}
            dt_moltype['mols'] = [deepcopy(mol) for molnr in range(dt_moltype['nmols'])]
            del dt_moltype['nmols']
            del dt_moltype['atoms']

    def load_gro_file_pos_vel(self, gro_filename):
        """Loads a positions and velocities from a gro file into an
        existing distinctive topology."""

        assert(hasattr(self, "distinctive_top"))

        atoms = self.atoms()
        with open(gro_filename, 'r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    pass
                elif i == 1:
                    pass
                elif i < self.natoms() + 2:
                    atom = atoms[i - 2]
                    atom.pos = np.array([float(line[20:28]),
                                         float(line[28:36]),
                                         float(line[36:44])])
                    # try to read velocities if present in file
                    try:
                        atom.vel = np.array([float(line[44:52]),
                                             float(line[52:60]),
                                             float(line[60:68])])
                    except IndexError:
                        atom.vel = np.full(3, np.nan)


    def load_gro_file_names(self, gro_filename):
        """Loads atom and mol names from a gro file into an
        existing distinctive topology."""

        assert(hasattr(self, "distinctive_top"))

        atoms = self.atoms()
        with open(gro_filename, 'r') as f:
            for i, line in enumerate(f):
                if i == 0:
                    pass
                elif i == 1:
                    pass
                elif i < self.natoms() + 2:
                    atom = atoms[i - 2]
                    atom.name = line[10:15].strip()


    def load_gro_file(self, gro_filename, atom_mass_dict, rot_treat_dict,
                      abc_indicators_dict, sigma_dict):
        """Loads a gro file into a distinctive topology with
positions and velocities."""
        atoms = []
        with open(gro_filename, 'r') as f:
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
                    atom['pos'] = np.array([float(line[20:28]),
                                            float(line[28:36]),
                                            float(line[36:44])])
                    # try to read velocities if present in file
                    try:
                        atom['vel'] = np.array([float(line[44:52]),
                                                float(line[52:60]),
                                                float(line[60:68])])
                    except IndexError:
                        atom['vel'] = np.full(3, np.nan)

                    atoms.append(atom)

        # create distinctive topology
        dt_moltypes = []
        dt_mols = []
        dt_atoms = []

        mol_nr_old = 1
        mol_name_old = atoms[0]['mol_name']

        dt_mol = {'atoms': []}
        dt_moltype = {'mols': []}

        # go throug atoms
        for atom_nr, atom in enumerate(atoms):
            # add new atom to dt_atoms
            dt_atom = {'name': atom['name'],
                       'pos': atom['pos'],
                       'vel': atom['vel'],
                       'mass': atom_mass_dict[atom['name']]}
            dt_atoms.append(dt_atom)

            # if last atom
            if atom_nr == natoms - 1:
                dt_mol = {'atoms': dt_atoms}
                dt_mols.append(dt_mol)
                mol_name = atom['mol_name']
                dt_moltype = {'name': mol_name,
                              'mols': dt_mols,
                              'rot_treat': rot_treat_dict[mol_name],
                              'abc_indicators': abc_indicators_dict[mol_name],
                              'sigma': sigma_dict[mol_name]}
                dt_moltypes.append(dt_moltype)
                break

            next_atom = atoms[atom_nr + 1]

            # if at end of molecule
            if next_atom['mol_nr'] != mol_nr_old:
                dt_mol = {'atoms': dt_atoms}
                dt_mols.append(dt_mol)
                dt_atoms = []
                mol_nr_old = next_atom['mol_nr']

            # if at end of moltype
            if next_atom['mol_name'] != mol_name_old:
                mol_name = atom['mol_name']
                dt_moltype = {'name': mol_name,
                              'mols': dt_mols,
                              'rot_treat': rot_treat_dict[mol_name],
                              'abc_indicators': abc_indicators_dict[mol_name],
                              'sigma': sigma_dict[mol_name]}
                dt_moltypes.append(dt_moltype)
                dt_mols = []
                mol_name_old = next_atom['mol_name']

        self.distinctive_top = dt_moltypes

    def save_gro_file(self, gro_filename, box):
        """Saves a distinctive topology in a gro file."""
        with open(gro_filename, 'w') as f:
            f.write("exported from SimTop\n")
            f.write(f"{self.natoms()}\n")
            for moltype in self.moltypes():
                for mol in moltype.mols():
                    for atom in mol.atoms():
                        f.write(f"{mol.index_in_top + 1:>5}")
                        f.write(f"{mol.name:>5}")
                        f.write(f"{atom.name:>5}")
                        f.write(f"{atom.index_in_top + 1:>5}")
                        try:
                            f.write(f"{atom.pos[0]:> 8.3f}")
                            f.write(f"{atom.pos[1]:> 8.3f}")
                            f.write(f"{atom.pos[2]:> 8.3f}")
                        except IndexError:
                            raise KeyError("""There are no positions in this
                                           topology. Not possible to write
                                           gro file!""")
                        try:
                            f.write(f"{atom.vel[0]:> 8.4f}")
                            f.write(f"{atom.vel[1]:> 8.4f}")
                            f.write(f"{atom.vel[2]:> 8.4f}")
                        except IndexError:
                            pass
                        f.write("\n")
            f.write(f"{box[0]} {box[1]} {box[2]}\n")

    def save_doscalc_parameters_file(self, nsamples, nblocks, nblocksteps,
                                     parameters_file="params.txt"):
        """
        Generate a parameters file for my dos-calc code.
        """

        with open(parameters_file, 'w') as f:

            f.write(str(nsamples) + "\n")
            f.write(str(nblocks) + "\n")
            f.write(str(nblocksteps) + "\n")
            f.write(str(self.nmoltypes()) + "\n")
            for moltype in self.moltypes():
                moltype_atommasses = (atomtype['mass'] for atomtype in moltype.atomtypes())
                moltype_string = f"{moltype.nmols} "\
                                 f"{moltype.natomtypes} "\
                                 f"{' '.join(map(str, moltype_atommasses))} "\
                                 f"{moltype.rot_treat} "\
                                 f"{' '.join(map(str, moltype.abc_indicators))}"

                f.write(moltype_string + '\n')


def com_atomlist(atomlist):
    com = np.zeros(3)
    mass = 0
    for atom in atomlist:
        com += atom.pos * atom.mass
        mass += atom.mass
    com /= mass
    return com


def move_atomlist(atomlist, vector):
    for atom in atomlist:
        atom.pos += vector


def decompose_velocities_of_molecule(molecule):
    """decompose the velocity of a molecule into translation,
    rotation and vibration parts."""

    positions = np.array([atom.pos for atom in molecule.atoms()])
    velocities = np.array([atom.vel for atom in molecule.atoms()])
    m_atommasses = np.array([atom.mass for atom in molecule.atoms()])

    if molecule.natoms == 1:
        return velocities, np.zeros((1, 3)), np.zeros((1, 3))

    center_of_mass = m_atommasses @ positions / molecule.mass
    mol_velocity_trn = m_atommasses @ velocities / molecule.mass
    velocities_trn = np.repeat([mol_velocity_trn], molecule.natoms, axis=0)
    positions_rel = positions - center_of_mass

    angular_momentum = m_atommasses @ np.cross(positions_rel, velocities)

    moi_tensor = np.zeros((3, 3))
    for i in range(len(positions_rel)):
        moi_tensor += m_atommasses[i] * ((positions_rel[i] @ positions_rel[i]) * np.identity(3)
                                         - np.tensordot(positions_rel[i], positions_rel[i], axes=0))

    if molecule.natoms == 2:
        raise Exception("linear molecules not implemented")
    else:
        angular_velocity = linalg.solve(moi_tensor, angular_momentum)

    velocities_rot = np.cross(angular_velocity, positions_rel)
    velocities_vib = velocities - velocities_trn - velocities_rot

    return velocities_trn, velocities_rot, velocities_vib

def calc_kinetic_energy_distribution(topology, kT):
    """calculate how much kinetic energy is in translation, rotation and vibration
    in units of kT"""
    e_kin_trn = 0
    e_kin_rot = 0
    e_kin_vib = 0
    e_kin_tot = 0

    for mol in topology.mols():
        vel_trn, vel_rot, vel_vib = decompose_velocities_of_molecule(mol)
        m_atommasses = np.array([atom.mass for atom in mol.atoms()])

        e_kin_trn += np.sum(m_atommasses @ vel_trn**2)
        e_kin_rot += np.sum(m_atommasses @ vel_rot**2)
        e_kin_vib += np.sum(m_atommasses @ vel_vib**2)

    e_kin_trn /= topology.nmols() * kT
    e_kin_rot /= topology.nmols() * kT
    e_kin_vib /= topology.nmols() * kT

    return e_kin_trn, e_kin_rot, e_kin_vib

def show_kinetic_energy_distribution(topology, kT):
    """prints out how much kinetic energy is in translation, rotation and vibration
    in units of kT"""

    e_kin_trn, e_kin_rot, e_kin_vib = calc_kinetic_energy_distribution(topology, kT)

    print(f"""kinetic energy distribution:
trn: {e_kin_trn:7.4f}
rot: {e_kin_rot:7.4f}
vib: {e_kin_vib:7.4f}""")
