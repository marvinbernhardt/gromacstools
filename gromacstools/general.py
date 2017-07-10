import subprocess
import os


def run_bash(command, stdin=None, logging=False):
    """
    Runs a command in a bash shell and optionally logs the output in log/other_command.log.

    Raises an exception if there is a non-zero return code and prints stdout and stderr in that case.
    Gromacs commands like 'gmx ...' get an own log file.
    """
    proc = subprocess.Popen(['bash', '-c', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    # write to logfile
    if logging:
        # create log directory if needed
        if not os.path.exists("log"):
            os.makedirs("log")

        # check if gmx command
        commands = command.split()
        if commands[0] == "gmx" and len(commands) > 1:
            logfile_name = "log/gmx_" + commands[1] + ".log"
        else:
            logfile_name = "log/other_command.log"

        # print log to file
        with open(logfile_name, "w+") as logfile:
            logfile.write(stdout.decode('utf-8'))
            logfile.write(stderr.decode('utf-8'))

    # raise error if returncode != 0
    if proc.returncode:
        print(stdout.decode())
        print(stderr.decode())
        raise Exception("There was a non-zero return code.")

    return stdout.decode()


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
