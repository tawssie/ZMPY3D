# GPL 3.0 License, volume
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

def fill_voxel_by_weight_density04(xyz, aa_name_list, residue_weight_map, grid_width, residue_box):
    min_bbox_point = np.min(xyz, axis=0)
    max_bbox_point = np.max(xyz, axis=0)
    dimension_bbox_unscaled = max_bbox_point - min_bbox_point

    max_box_edge = max([box.shape[0] for box in residue_box.values()])

    dimension_bbox_scaled = np.ceil((dimension_bbox_unscaled / grid_width) + max_box_edge).astype(int)
    corner_xyz = min_bbox_point - max_box_edge * grid_width / 2

    weight_multiplier = 1
    num_of_atom = xyz.shape[0]

    voxel3d = np.zeros(dimension_bbox_scaled)

    for i in range(num_of_atom):
        aa_name = aa_name_list[i]
        if aa_name not in residue_weight_map:
            aa_name = 'ASP'

        coord = xyz[i, :]
        aa_box = residue_box[aa_name]
        box_edge = aa_box.shape[0]

        coord_box_corner = np.fix(np.round((coord - corner_xyz) / grid_width - box_edge / 2)).astype(int)

        start = coord_box_corner
        end = coord_box_corner + box_edge

        voxel3d[start[0]:end[0], start[1]:end[1], start[2]:end[2]] += aa_box * weight_multiplier

    return voxel3d, corner_xyz
