def get_total_mass(moltypes):
    return sum((sum((at['mass'] for at in mt['atoms'])) * mt['nmols']
                for mt in moltypes))


def get_nmols(moltypes):
    return sum((mt['nmols'] for mt in moltypes))


def get_average_molar_mass(moltypes):
    return get_total_mass(moltypes) / get_nmols(moltypes)


def get_natoms(moltypes):
    return sum((len(mt['atoms']) * mt['nmols'] for mt in moltypes))


def count_atomname(moltypes, atomname):
    return sum((sum((1 for at in mt['atoms'] if at['name'] == atomname)) * mt['nmols']
               for mt in moltypes))


def get_mol_mass(moltypes, moltype_index, single_mol=False):
    mt = moltypes[moltype_index]
    return sum((at['mass'] for at in mt['atoms'])) * (1 if single_mol else mt['nmols'])
