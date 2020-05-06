import re


def get_parameter(mdp_file, key):
    """retruns the value of a parameter in a .mdp file"""

    search_string = r"{0}\s*=\s*(.*)$".format(key)
    with open(mdp_file, "r") as file:
        for line in file:
            if line.startswith(key):
                return re.match(search_string, line).group(1)
                break


def set_parameter(mdp_file, key, value):
    """sets the value of a parameter in a .mdp file"""

    search_string = r"{0}\s*=.*".format(key)
    search_pattern = re.compile(search_string)
    replacement_string = r"{0} = {1}".format(key, value)
    changed = False
    with open(mdp_file, "r") as file:
        newlines = []
        for line in file:
            if re.search(search_pattern, line):
                newlines.append(re.sub(search_string, replacement_string, line))
                changed = True
            else:
                newlines.append(line)

    if not changed:
        newlines.append(replacement_string + '\n')

    with open(mdp_file, "w") as file:
        for line in newlines:
            file.write(line)
