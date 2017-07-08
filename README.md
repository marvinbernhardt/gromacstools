# gromacs-tools

## available functions:

### general
```python
run_bash(command, stdin=None, logging=False):
```

### gro
```python
def get_box(gro_file)
def get_natoms(gro_file)
def mix_water(gro_file, nmols)
def read_moltypes(gro_file, atom_mass_dict, abc_indicators_dict, sigma_dict)
def set_atomname(gro_file, atomlist, atomname)
def set_coordinate(gro_file, atom, coord, value)
def set_molname(gro_file, atomlist, molname)
def translate_with_pbc(gro_file, vector, out_file="conf-shifted.gro")
def generate_ffc_crystal_atomic(file, natoms, boxlength, mol_name, atom_name)
def generate_ffc_crystal_water(filename, n_mols, box_length)
```

### mdp
```python
def get_parameter(mdp_file, key)
def set_parameter(mdp_file, key, value)
```

### mol
```python
def bounding_box(mol)
def calc_com(mol)
```

### moltypes
```python
def generate_parameters_file(moltypes, nsamples, nblocks, nblocksteps, parameters_file="params.txt")
def shrink_box_remove_mols(moltypes, box_center, box_vectors, padding)
def write_gro_file(moltypes, out_filename, box)
```

### remote
```python
def check_slurm_job(jobid, remote_host)
def pull_files(filelist, remote_host, remote_dir, exclude="")
def push_files(filelist, remote_host, remote_dir, exclude="")
def run_slurm_array(command, remote_host, remote_dir, array_start, array_end, array_step=1, name="name", time_minutes="600", mem_per_cpu=1750, cpus_per_task=None, dry_run=False)
def run_slurm(command, remote_host, remote_dir, name="name", time_minutes="600", mem_per_cpu=1750, cpus_per_task=None, dry_run=False)
```

### xvg
```python
def plot(filename)
def load(filename)
```
