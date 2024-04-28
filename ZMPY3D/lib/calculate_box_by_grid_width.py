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


import numpy as np

def calculate_box_by_grid_width(param, grid_width):
    residue_name = param['residue_name']
    residue_weight_map = param['residue_weight_map']
    residue_radius_map = param['residue_radius_map']

    sd_cutoff = param['sd_cutoff']
    three_over_2pi_cube_root = param['three_over_2pi_32']
    density_multiplier = param['density_multiplier']

    aa_box_list = []

    for aa in residue_name:
        weight = residue_weight_map[aa]
        radius = residue_radius_map[aa]

        sigma = np.sqrt(radius ** 2 / 5.0)
        box_edge = int(np.ceil(2 * sd_cutoff * sigma / grid_width))

        if box_edge % 2 == 0:
            box_edge += 1

        center = box_edge // 2
        sqr_radius = center ** 2

        x, y, z = np.mgrid[0:box_edge, 0:box_edge, 0:box_edge] - center

        is_within_radius = x**2 + y**2 + z**2 <= sqr_radius

        x_sx = (x**2 + y**2 + z**2) * grid_width**2 / (radius**2 / 5.0)
        gaus_val = (weight * three_over_2pi_cube_root / (radius**2 / 5.0 * sigma) *
                    np.exp(-0.5 * x_sx) * density_multiplier)

        residue_unit_box = np.zeros((box_edge, box_edge, box_edge))
        residue_unit_box[is_within_radius] = gaus_val[is_within_radius]

        aa_box_list.append(residue_unit_box)

    return aa_box_list

