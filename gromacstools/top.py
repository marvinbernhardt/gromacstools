import numpy as np


class Atom:
    """Atom class. Can contain positions and velocities."""

    def __init__(self, dt_atom, index_in_top, index_in_moltype, index_in_mol):
        self.name = dt_atom['name']
        self.mass = dt_atom['mass']
        self.index_in_top = index_in_top
        self.index_in_moltype = index_in_moltype
        self.index_in_mol = index_in_mol
        try:
            self.pos = dt_atom['pos']
        except:
            pass
        try:
            self.vel = dt_atom['vel']
        except:
            pass


class Molecule:
    """Molecule class. Contains atoms() function."""
    def __init__(self, dt_mol, name, index_in_top, index_in_moltype,
                 firstatom, moltype_nr):
        self.name = name
        self.dt_atoms = dt_mol['atoms']
        self.mass = sum([dt_atom['mass'] for dt_atom in self.dt_atoms])
        self.natoms = len(self.dt_atoms)
        self.index_in_top = index_in_top
        self.index_in_moltype = index_in_moltype
        self.firstatom = firstatom
        self.moltype_nr = moltype_nr

    def atoms(self):
        """Return list of Atom objects in this Molecule."""
        atoms = []
        for index_in_mol, dt_atom in enumerate(self.dt_atoms):
            index_in_moltype = (self.index_in_moltype * self.natoms
                                + index_in_mol)
            index_in_top = self.firstatom + index_in_mol
            atoms.append(Atom(dt_atom, index_in_top, index_in_moltype,
                              index_in_mol))
        return atoms


class Moltype:
    """Moltype class. Contains mols() and atoms() function."""
    def __init__(self, dt_moltype, index_in_top, firstmol, firstatom):
        self.name = dt_moltype['name']
        self.dt_mols = dt_moltype['mols']
        self.abc_indicators = dt_moltype['abc_indicators']
        self.nmols = len(self.dt_mols)
        self.natomtypes = len(self.atomtypes())
        self.natoms = self.nmols * self.natomtypes
        self.mass = self.nmols * sum([atomtype['mass'] for atomtype
                                      in self.atomtypes()])
        self.index_in_top = index_in_top
        self.firstmol = firstmol
        self.firstatom = firstatom

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
        self.distinctive_top = simple_top.copy()

        # expand nmols to mols
        for dt_moltype in self.distinctive_top:
            mol = {'atoms': dt_moltype['atomtypes']}
            mols = [mol for molnr in range(dt_moltype['nmols'])]
            dt_moltype['mols'] = mols
            del dt_moltype['nmols']
            del dt_moltype['atomtypes']

    def load_gro_file(self, gro_filename, atom_mass_dict, abc_indicators_dict,
                      sigma_dict):
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
                    except:
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
                        f.write(f"{mol.index_in_top:>5}")
                        f.write(f"{mol.name:>5}")
                        f.write(f"{atom.name:>5}")
                        f.write(f"{atom.index_in_top:>5}")
                        try:
                            f.write(f"{atom.pos[0]:> 8.3f}")
                            f.write(f"{atom.pos[1]:> 8.3f}")
                            f.write(f"{atom.pos[2]:> 8.3f}")
                        except:
                            raise KeyError("""There are no positions in this
                                           topology. Not possible to write
                                           gro file!""")
                        try:
                            f.write(f"{atom.vel[0]:> 8.4f}")
                            f.write(f"{atom.vel[1]:> 8.4f}")
                            f.write(f"{atom.vel[2]:> 8.4f}")
                        except:
                            pass
                        f.write("\n")
            f.write(f"{box[0]} {box[1]} {box[2]}\n")

    def save_2pt_parameters_file(self, nsamples, nblocks, nblocksteps,
                                 parameters_file="params.txt"):
        """
        Generate a parameters file for my 2PT code.
        """

        def fwrite_int(number):
            f.write(str(number) + "\n")

        def fwrite_list(_list):
            f.write(' '.join(map(str, _list)) + "\n")

        with open(parameters_file, 'w') as f:

            fwrite_int(nsamples)
            fwrite_int(nblocks)
            fwrite_int(nblocksteps)
            fwrite_int(self.natoms())
            fwrite_int(self.nmols())
            fwrite_int(self.nmoltypes())
            fwrite_list([moltype.firstmol for moltype in self.moltypes()])
            fwrite_list([moltype.firstatom for moltype in self.moltypes()])
            fwrite_list([moltype.nmols for moltype in self.moltypes()])
            fwrite_list([moltype.natomtypes for moltype in self.moltypes()])
            fwrite_list([' '.join(map(str, moltype.abc_indicators))
                         for moltype in self.moltypes()])
            fwrite_list([mol.firstatom for mol in self.mols()])
            fwrite_list([mol.natoms for mol in self.mols()])
            fwrite_list([mol.mass for mol in self.mols()])
            fwrite_list([mol.moltype_nr for mol in self.mols()])
            fwrite_list([atom.mass for atom in self.atoms()])


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
