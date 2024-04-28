# GPL3.0 License, JSL from ZM, this is my originality, to calculate 3D bbox moment by tensordot
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

def calculate_bbox_moment06(voxel3d, max_order, xyz_sample_struct):
    extend_voxel3d = np.zeros(np.array(voxel3d.shape) + 1)
    extend_voxel3d[:-1, :-1, :-1] = voxel3d

    diff_extend_voxel3d = np.diff(np.diff(np.diff(extend_voxel3d, axis=0), axis=1), axis=2)

    x_power = np.power(xyz_sample_struct['X_sample'][1:, np.newaxis], np.arange(1, max_order + 2))
    y_power = np.power(xyz_sample_struct['Y_sample'][1:, np.newaxis], np.arange(1, max_order + 2))
    z_power = np.power(xyz_sample_struct['Z_sample'][1:, np.newaxis], np.arange(1, max_order + 2))

    bbox_moment = np.tensordot(
        z_power, 
        np.tensordot(
            y_power, 
            np.tensordot(x_power, diff_extend_voxel3d, axes=([0], [0])), 
            axes=([0], [1])
        ), 
        axes=([0], [2])
    )

    p, q, r = np.meshgrid(np.arange(1, max_order + 2), np.arange(1, max_order + 2), np.arange(1, max_order + 2), indexing='ij')
    bbox_moment = -np.transpose(bbox_moment, (2, 1, 0)) / p / q / r

    volume_mass = bbox_moment[0, 0, 0]
    center = [bbox_moment[1, 0, 0], bbox_moment[0, 1, 0], bbox_moment[0, 0, 1]] / volume_mass

    return volume_mass, center, bbox_moment
