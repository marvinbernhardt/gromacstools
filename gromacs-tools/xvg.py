import shlex
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

def plot(filename):
    """plots a grace file. Does not work with all grace features!"""
    data, header = load(filename)
    data = np.array(data)

    if header['type'].lower() == 'xy':
        plt.figure(figsize=(10,3))
        for column in range(len(data[0]) - 1):
            plt.plot(data[:,0], data[:,column + 1], label=header['legend'][column + 1])
        plt.xlabel(header['xlabel'])
        plt.ylabel(header['ylabel'])
        plt.legend()
        plt.show()

    elif header['type'].lower() == 'xydy':
        plt.figure(figsize=(10,3))
        plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], label=header['ylabel'])
        plt.xlabel(header['xlabel'])
        plt.ylabel(header['ylabel'])
        plt.legend()
        plt.show()


def load(filename):
    """parses a grace file and returns its data and header info. Does not work with all grace features!"""
    grace_types_legend = {
        'XY': ['x', 'y'],
        'XYDX': ['x', 'y', 'dx'],
        'XYDY': ['x', 'y', 'dy'],
        'XYDXDX': ['x', 'y', '-dx', 'dx'],
        'XYDYDY': ['x', 'y', '-dy', 'dy'],
        'XYDXDY': ['x', 'y', 'dx', 'dy'],
        'XYDXDXDYDY': ['x', 'y', '-dx', 'dx', '-dy', 'dy'],
        'XYZ': ['x', 'y', 'z'],
        'XYR': ['x', 'y', 'r'],
    }

    def convert_xvg_string(string):
        replace_list = [
            (r'\xl\f{}', 'λ'),
            (r'\xD\f{}', 'Δ')
        ]
        for pattern in replace_list:
            string = string.replace(pattern[0], pattern[1])
        return string

    header = {}
    pure_data = StringIO()
    with open(filename, "r") as f:
        for line in f:
            if line.startswith('#'):
                pass
            elif line.startswith('@'):
                line = convert_xvg_string(line)
                parts = shlex.split(line)
                if parts[0] == "@TYPE":
                    header['type'] = parts[1]
                    header['legend'] = grace_types_legend[header['type'].upper()]
                elif parts[1] == "title" and parts[2].startswith('"'):
                    header['title'] = ' '.join(parts[2:])
                elif parts[1:3] == ["xaxis", "label"]:
                    header['xlabel'] = ' '.join(parts[3:])
                elif parts[1:3] == ["yaxis", "label"]:
                    header['ylabel'] = ' '.join(parts[3:])
                elif parts[1][0] == 's' and parts[2] == "legend":
                    column = int(parts[1][1:]) + 1 # grace does not count x column
                    if column < len(header['legend']):
                        header['legend'][column] = ' '.join(parts[3:])
                    else:
                        header['legend'].append(' '.join(parts[3:]))
            else:
                pure_data.write(line.strip() + "\n")

    # header post-processing
    header['legend'][0] = header['xlabel']

    pure_data.seek(0)
    dataFrame = pd.read_csv(pure_data, delim_whitespace=True, header=None)
    dataFrame.columns = header['legend']
    return dataFrame, header
