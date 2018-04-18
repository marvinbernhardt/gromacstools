import subprocess
import shlex
from .general import run_bash


def check_slurm_job(jobid, remote_host):
    """Checks the status of a slurm job"""
    proc = subprocess.Popen(['ssh', remote_host, *shlex.split('"source /etc/profile; sacct -nPj {} -o State"'.format(jobid))], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    try:
        return stdout.decode().split()
    except IndexError:
        print("No output - wrong jobid?")
        return None


def check_pbs_job(jobid, remote_host):
    """Checks the status of a PBS job"""
    proc = subprocess.Popen(['ssh', remote_host, *shlex.split(f'"qstat -f {jobid}"')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    try:
        output = stdout.decode().strip()
        return next((line for line in output.splitlines() if 'job_state = ' in line)).split(' = ')[1:]
    except IndexError:
        print("No output - wrong jobid?")
        return None


def pull_files(filelist, remote_host, remote_dir, exclude=""):
    """Copies a list of files from a directory on a remote host to the currend directory."""
    filelist_string = ' :{}/'.format(remote_dir).join(filelist)
    if exclude == "":
        run_bash(f'rsync -az {remote_host}:{remote_dir}/{filelist_string} ./', logging=False)
    else:
        run_bash(f'rsync -az --exclude={exclude} {remote_host}:{remote_dir}/{filelist_string} ./', logging=False)


def push_files(filelist, remote_host, remote_dir, exclude=""):
    """Copies a list of files on a remote host into a specified directory."""
    filelist_string = ' '.join(filelist)
    if exclude == "":
        run_bash(f'rsync -az {filelist_string} {remote_host}:{remote_dir}/', logging=False)
    else:
        run_bash(f'rsync -az --exclude="{exclude}" {filelist_string} {remote_host}:{remote_dir}/', logging=False)


def run_slurm_array(command, remote_host, remote_dir, array_start, array_end, array_step=1, name="name", time_minutes="600", mem_per_cpu=1750, cpus_per_task=None, dry_run=False):
    """Runs a slurm job array on a remote host. Settings should be good for OpenMP."""
    if cpus_per_task is not None:
        cpus_per_task_line = f"#SBATCH --cpus-per-task={cpus_per_task}\n"
    else:
        cpus_per_task_line = ""
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --array={array_start}-{array_end}:{array_step}
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={time_minutes}
{cpus_per_task_line}
echo -n "started $SLURM_JOB_ID at "; date
echo -n "running on "; hostname

shopt -s expand_aliases
source $HOME/.bash_profile
mkdir -p $HPC_LOCAL
set -euo pipefail

{command}

echo -n "finished $SLURM_JOB_ID at "; date
"""
    # write and copy to remote
    with open("sbatch.sh", 'w') as f:
        f.write(sbatch_script)
    run_bash("rsync -az sbatch.sh {0}:{1}/".format(remote_host, remote_dir), logging=False)

    # run sbatch on remote and remote jobid
    if not dry_run:
        proc = subprocess.Popen(['ssh', remote_host, *shlex.split('"source /etc/profile; cd {0}; sbatch --parsable sbatch.sh"'.format(remote_dir))], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        stderr_string = stderr.decode().rstrip()
        if stderr_string.isdigit():
            return stderr_string
        else:
            return None
    else:
        return None


def run_slurm(command, remote_host, remote_dir, name="name", time_minutes="600", mem_per_cpu=1750, cpus_per_task=None, dry_run=False):
    """Runs a slurm job on a remote host. Settings should be good for OpenMP."""
    if cpus_per_task is not None:
        cpus_per_task_line = f"#SBATCH --cpus-per-task={cpus_per_task}\n"
    else:
        cpus_per_task_line = ""
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={time_minutes}
{cpus_per_task_line}
echo -n "started $SLURM_JOB_ID at "; date
echo -n "running on "; hostname

shopt -s expand_aliases
source $HOME/.bash_profile
mkdir -p $HPC_LOCAL
set -euo pipefail

{command}

echo -n "finished $SLURM_JOB_ID at "; date
"""
    # write and copy to remote
    with open("sbatch.sh", 'w') as f:
        f.write(sbatch_script)
    run_bash("rsync -az sbatch.sh {0}:{1}/".format(remote_host, remote_dir), logging=False)

    # run sbatch on remote and remote jobid
    if not dry_run:
        proc = subprocess.Popen(['ssh', remote_host, *shlex.split('"source /etc/profile; cd {0}; sbatch --parsable sbatch.sh"'.format(remote_dir))], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        stderr_string = stderr.decode().rstrip()
        if stderr_string.isdigit():
            return stderr_string
        else:
            return None
    else:
        return None


def run_pbs(command, remote_host, remote_dir, name="name",
            walltime="01:00:00", queue='chickencurry',
            nodes=1, dry_run=False):
    """Runs a pbs job on a remote host."""
    if queue == "chickencurry":
        nodes_ppn_line = "#PBS -l nodes=1:ppn=24"
    elif queue == "pizza":
        nodes = nodes * 24  # hack since enzo names cores nodes in this case
        nodes_ppn_line = f"#PBS -l nodes={nodes}"
    else:
        raise KeyError(f"unknown queue {queue}")
    qsub_script = f"""#!/bin/bash
#PBS -j oe
#PBS -q {queue}
#PBS -l walltime={walltime}
{nodes_ppn_line}
#PBS -N {name}

if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi

echo -n "started $PBS_JOBID at "; date
echo -n "running on "; hostname

shopt -s expand_aliases
source $HOME/.bash_profile
set -euo pipefail

{command}

echo -n "finished $PBS_JOBID at "; date
"""
    # write and copy to remote
    with open("qsub.sh", 'w') as f:
        f.write(qsub_script)
    run_bash(f"rsync -az qsub.sh {remote_host}:{remote_dir}/", logging=False)

    # run qsub on remote and remote jobid
    if not dry_run:
        proc = subprocess.Popen(['ssh', remote_host, *shlex.split(f'"source /etc/profile; cd {remote_dir}; qsub qsub.sh"')], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        stdout_string = stdout.decode().rstrip()
        jobid = stdout_string.split('.')[0]
        if jobid.isdigit():
            return jobid
        else:
            return None
    else:
        return None
