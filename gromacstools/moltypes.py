import numpy as np
from copy import deepcopy


def generate_parameters_file(
    moltypes, nsamples, nblocks, nblocksteps, parameters_file="params.txt"
):
    """
    Generate a parameters file for my 2PT code.
    """

    moltypes = deepcopy(moltypes)

    nmoltypes = len(moltypes)
    nmols = sum([len(moltype["mols"]) for moltype in moltypes])
    natoms = sum([len(mol) for moltype in moltypes for mol in moltype["mols"]])

    # moltypes
    for moltype in moltypes:
        moltype["nmols"] = len(moltype["mols"])
        moltype["mass"] = sum([atom["mass"] for atom in moltype["mols"][0]])
        moltype["natomtypes"] = len(moltype["mols"][0])
        moltype["natoms"] = moltype["natomtypes"] * moltype["nmols"]
    for i in range(nmoltypes):
        moltypes[i]["firstmol"] = sum([moltype["nmols"] for moltype in moltypes[0:i]])
        moltypes[i]["firstatom"] = sum([moltype["natoms"] for moltype in moltypes[0:i]])

    # mols
    mols = []
    for h, moltype in enumerate(moltypes):
        for i, mol in enumerate(moltype["mols"]):
            mols.append(
                {
                    "name": moltype["name"],
                    "mass": moltype["mass"],
                    "natoms": len(mol),
                    "moltypenr": h,
                }
            )
    for i in range(nmols):
        mols[i]["firstatom"] = sum([mol["natoms"] for mol in mols[0:i]])

    # atoms
    atoms = []
    for moltype in moltypes:
        for mol in moltype["mols"]:
            for atom in mol:
                atoms.append(atom)

    with open(parameters_file, "w") as f:

        def fwrite(obj):
            f.write(obj + "\n")

        fwrite(str(nsamples))
        fwrite(str(nblocks))
        fwrite(str(nblocksteps))
        fwrite(str(natoms))
        fwrite(str(nmols))
        fwrite(str(nmoltypes))
        fwrite(" ".join(map(str, [moltype["firstmol"] for moltype in moltypes])))
        fwrite(" ".join(map(str, [moltype["firstatom"] for moltype in moltypes])))
        fwrite(" ".join(map(str, [moltype["nmols"] for moltype in moltypes])))
        fwrite(" ".join(map(str, [moltype["natomtypes"] for moltype in moltypes])))
        fwrite(
            " ".join(
                map(
                    str,
                    [
                        " ".join(map(str, moltype["abc_indicators"]))
                        for moltype in moltypes
                    ],
                )
            )
        )
        fwrite(" ".join(map(str, [mol["firstatom"] for mol in mols])))
        fwrite(" ".join(map(str, [mol["natoms"] for mol in mols])))
        fwrite(" ".join(map(str, [mol["mass"] for mol in mols])))
        fwrite(" ".join(map(str, [mol["moltypenr"] for mol in mols])))
        fwrite(" ".join(map(str, [atom["mass"] for atom in atoms])))

    return moltypes


def shrink_box_remove_mols(moltypes, box_center, box_vectors, padding):
    """
    Returns a new moltypess list and a new box; molecules are all within that box.

    padding space is added at the positive limits of the box for easier restarting of a
    simulation
    """
    box_center = np.array(box_center)
    box_vectors = np.array(box_vectors)

    box_vertex0 = box_center - box_vectors / 2
    box_vertex1 = box_center + box_vectors / 2

    newbox = box_vectors + padding
    newmoltypes = deepcopy(moltypes)

    # identifying mols that have a com out of the box
    for moltype in newmoltypes:
        mol_index_out_of_box = []
        for molnr, mol in enumerate(moltype["mols"]):
            com = mol_calc_com(mol)

            if any(com - box_vertex0 < 0):
                mol_index_out_of_box.append(molnr)
            elif any(com - box_vertex1 > 0):
                mol_index_out_of_box.append(molnr)

        # remove those mols
        newmols = [
            mol
            for i, mol in enumerate(moltype["mols"])
            if i not in mol_index_out_of_box
        ]
        print("deleted {} {}".format(len(mol_index_out_of_box), moltype["name"]))
        moltype["mols"] = newmols

        # shifting molecules
        for mol in moltype["mols"]:
            for atom in mol:
                atom["r"] -= box_vertex0

    return newmoltypes, newbox


def write_gro_file(moltypes, out_filename, box):
    """writes mols array to a gro file"""
    with open(out_filename, "w") as f:
        f.write("generated from python\n")
        f.write(
            str(
                len(
                    [
                        atom
                        for moltype in moltypes
                        for mol in moltype["mols"]
                        for atom in mol
                    ]
                )
            )
            + "\n"
        )

        atomnr = 1
        molnr = 1
        for moltype in moltypes:
            for mol in moltype["mols"]:
                for atom in mol:
                    string = (
                        "{:>5d}{:<5s}{:>5s}{:>5d}{:8.3f}{:8.3f}{:8.3f}"
                        "{:8.4f}{:8.4f}{:8.4f}\n".format(
                            molnr,
                            moltype["name"],
                            atom["name"],
                            atomnr,
                            atom["r"][0],
                            atom["r"][1],
                            atom["r"][2],
                            atom["v"][0],
                            atom["v"][1],
                            atom["v"][2],
                        )
                    )
                    f.write(string)
                    atomnr += 1
                molnr += 1

        f.write("{:10.5f}{:10.5f}{:10.5f}\n".format(box[0], box[1], box[2]))
