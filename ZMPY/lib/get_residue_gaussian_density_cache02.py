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

from .calculate_box_by_grid_width import *

def get_residue_gaussian_density_cache02(param):
    residue_name = param['residue_name']
    doable_grid_width_list = param['doable_grid_width_list']

    aa_box_list_gw_1 = calculate_box_by_grid_width(param, 1.00)

    total_gaussian_density_gw_1 = [np.sum(box) for box in aa_box_list_gw_1]


    aa_box_map_ref = dict(zip(residue_name, aa_box_list_gw_1))

    grid_widths = [1] + doable_grid_width_list

    aa_box_map_list = [aa_box_map_ref]

    for gw in grid_widths[1:]:

        temp_aa_box_list = calculate_box_by_grid_width(param, gw)

        temp_total_gaussian_density = [np.sum(box) for box in temp_aa_box_list]

        density_scalar = [td_gw1 / (td * gw**3) for td_gw1, td in zip(total_gaussian_density_gw_1, temp_total_gaussian_density)]

        adjusted_temp_aa_box_list = [box * scalar for box, scalar in zip(temp_aa_box_list, density_scalar)]

        aa_box_map_list.append(dict(zip(residue_name, adjusted_temp_aa_box_list)))

    residue_box = dict(zip(grid_widths, aa_box_map_list))

    return residue_box
