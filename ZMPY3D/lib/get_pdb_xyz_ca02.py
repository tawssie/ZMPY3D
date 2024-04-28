# MIT License
#
# Copyright (c) 2024 Jhih-Siang Lai
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import math
import numpy as np

def get_pdb_xyz_ca02(file_name):

    with open(file_name, 'r') as file:
        lines = file.readlines()

    atom_lines = [line for line in lines if line.startswith('ATOM')]

    ca_lines = [line for line in atom_lines if line[13:15] == 'CA']

    xyz = []
    aa_names = []
    for line in ca_lines:
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        aa_name = line[17:20].strip()

        xyz.append((x, y, z))
        aa_names.append(aa_name)

        if any(map(math.isnan, [x, y, z])):
            raise ValueError('has nan in XYZ')

    xyz_matrix = np.array(xyz)

    return xyz_matrix, aa_names

