# GPL3.0 License, volume
#
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


from .get_residue_weight_map01 import *
from .get_residue_radius_map01 import *
import math

def get_global_parameter02():
    param = {}

    param['doable_grid_width_list'] = [0.25, 0.5, 2, 4, 8, 16]
    param['default_radius_multiplier'] = 1.8

    sd_cutoff = 3.0
    three_over_2pi_32 = (3 / (2 * math.pi)) ** (3 / 2)
    density_multiplier = 100


    residue_name = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'MSE', 'PHE', 'PRO',
        'SER', 'THR', 'TRP', 'TYR', 'VAL', 'A', 'T', 'C', 'G', 'U', 'I',
        'DA', 'DT', 'DG', 'DC', 'DU', 'DI'
    ]


    residue_weight_map = get_residue_weight_map01()
    residue_radius_map = get_residue_radius_map01()

    param['sd_cutoff'] = sd_cutoff
    param['three_over_2pi_32'] = three_over_2pi_32
    param['density_multiplier'] = density_multiplier
    param['residue_name'] = residue_name
    param['residue_weight_map'] = residue_weight_map
    param['residue_radius_map'] = residue_radius_map

    return param

