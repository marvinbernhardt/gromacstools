import subprocess
import os


def run_bash(command, stdin=None, logging=False, print_stdout=False,
             print_stderr=False):
    """
    Runs a command in a bash shell and optionally prints output or logs the
    output in log/other_command.log.

    Returns stdout as string.

    Raises an exception if there is a non-zero return code and prints stdout
    and stderr in that case.

    'set -euo pipefail' is put before the command, in order to stop early if
    there is an error.
    """

    # bash strict mode, always!
    command_full = "set -euo pipefail\n" + command

    cp = subprocess.run(command_full, shell=True, check=True,
                        capture_output=True, text=True, input=stdin)
    cp.check_returncode()

    # write to logfile
    if logging:
        # create log directory if needed
        if not os.path.exists("log"):
            os.makedirs("log")

        # check if gmx command
        commands = command.strip().split()
        logfile_name = f"log/{commands[0]}.log"

        # print log to file
        with open(logfile_name, "w+") as logfile:
            logfile.write(cp.stdout)
            logfile.write(cp.stderr)

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
