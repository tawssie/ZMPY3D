# GPL3.0 License, volume
# Copyright (C) 2024  Jhih-Siang Lai
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import numpy as np

def get_total_residue_weight(aa_name_list, residue_weight_map):

    weight_multiplier = 1
    
    total_residue_weight = 0
    for aa_name in aa_name_list:
        if aa_name not in residue_weight_map:
            aa_name = 'ASP'
        total_residue_weight += residue_weight_map[aa_name] * weight_multiplier

    return total_residue_weight

