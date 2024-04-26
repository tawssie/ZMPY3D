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

from .read_file_as_string_list import read_file_as_string_list
from .write_string_list_as_file import write_string_list_as_file

import numpy as np

def set_pdb_xyz_rot_m_01(input_file_name, rot_m, output_file_name):
    def sprintf2(format_spec, numbers):
        return [format_spec % num for num in numbers]

    contents = read_file_as_string_list(input_file_name)

    x = [float(s[30:38]) for s in contents if s.startswith('ATOM')]
    y = [float(s[38:46]) for s in contents if s.startswith('ATOM')]
    z = [float(s[46:54]) for s in contents if s.startswith('ATOM')]
    
    target_xyz = np.array([x, y, z]).T
    
    temp_xyz = rot_m @ np.hstack((target_xyz, np.ones((target_xyz.shape[0], 1)))).T
    temp_xyz = np.round(temp_xyz, decimals=3)
    
    x_formatted = sprintf2('%8.3f', temp_xyz[0, :])
    y_formatted = sprintf2('%8.3f', temp_xyz[1, :])
    z_formatted = sprintf2('%8.3f', temp_xyz[2, :])
    

    atom_count = 0
    for i, s in enumerate(contents):
        if s.startswith("ATOM"):
            contents[i] = s[:30] + x_formatted[atom_count] + y_formatted[atom_count] + z_formatted[atom_count] + s[54:]
            atom_count += 1
    
    write_string_list_as_file(contents, output_file_name)

