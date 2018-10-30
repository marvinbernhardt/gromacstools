import subprocess
import os


def run_bash(command, stdin=None, print_stdout=False, print_stderr=False):
    """
    Runs a command in a sh (not bash :/) shell and optionally prints output.

    Returns stdout as string.

    Raises an exception if there is a non-zero return code and prints returncode,
    stdout and stderr in that case.

    'set -euo pipefail' is run before the command, in order to stop early if
    there is an error.
    """

    # bash strict mode, always!
    command_full = "set -euo pipefail\n" + command

    try:
        cp = subprocess.run(command_full, shell=True, check=True,
                            capture_output=True, text=True, input=stdin)
    except subprocess.CalledProcessError as e:
        print("#=#=# The command #=#=#")
        print(e.cmd)
        print("#=#=# returned #=#=#")
        print(e.returncode)
        print("#=#=# stdout is #=#=#")
        print(e.stdout)
        print("#=#=# stderr is #=#=#")
        print(e.stderr)
        raise

    # print
    if print_stdout:
        print(cp.stdout)
    if print_stderr:
        print(cp.stderr)

    return cp.stdout


class WorkingDir():
    """
    Context Manager to execute code in a certain directory.

    path is a directory.
    Usage:
    >>> with WorkingDir(path):
    ...     do something here
    """
    def __init__(self, path):
        self.path = os.path.expanduser(path)
        self.old_directory = os.getcwd()
        os.makedirs(self.path, exist_ok=True)

    def __enter__(self):
        os.chdir(self.path)

    def __exit__(self, *args):
        os.chdir(self.old_directory)
