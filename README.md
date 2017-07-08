# gromacs-tools

## available functions:

general
run_bash(command, stdin=None, logging=False):

gro
get_box(gro_file)
get_natoms(gro_file)
mix_water(gro_file, nmols)
read_moltypes(gro_file, atom_mass_dict, abc_indicators_dict, sigma_dict)
set_atomname(gro_file, atomlist, atomname)
set_coordinate(gro_file, atom, coord, value)
set_molname(gro_file, atomlist, molname)
translate_with_pbc(gro_file, vector, out_file="conf-shifted.gro")
generate_ffc_crystal_atomic(file, natoms, boxlength, mol_name, atom_name)
generate_ffc_crystal_water(filename, n_mols, box_length)

mdp
get_parameter(mdp_file, key)
set_parameter(mdp_file, key, value)

mol
bounding_box(mol)
calc_com(mol)

moltypes
generate_parameters_file(moltypes, nsamples, nblocks, nblocksteps, parameters_file="params.txt")
shrink_box_remove_mols(moltypes, box_center, box_vectors, padding)
write_gro_file(moltypes, out_filename, box)

remote
check_slurm_job(jobid, remote_host)
pull_files(filelist, remote_host, remote_dir, exclude="")
push_files(filelist, remote_host, remote_dir, exclude="")
run_slurm_array(command, remote_host, remote_dir, array_start, array_end, array_step=1, name="name", time_minutes="600", mem_per_cpu=1750, cpus_per_task=None, dry_run=False)
run_slurm(command, remote_host, remote_dir, name="name", time_minutes="600", mem_per_cpu=1750, cpus_per_task=None, dry_run=False)

xvg
plot(filename)
load(filename)

