# GPL3.0 License, JSL
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

def calculate_molecular_radius03(voxel3d, center, volume_mass, default_radius_multiplier):
    has_weight = voxel3d > 0

    voxel_list = voxel3d[has_weight]

    x_coord, y_coord, z_coord = np.where(has_weight)
    
    voxel_list_xyz = np.stack([x_coord, y_coord, z_coord], axis=1)
    center = np.array(center)

    voxel_dist2center_squared = np.sum((voxel_list_xyz - center) ** 2, axis=1)

    average_voxel_mass2center_squared = np.sum(voxel_dist2center_squared * voxel_list) / volume_mass

    average_voxel_dist2center = np.sqrt(average_voxel_mass2center_squared) * default_radius_multiplier
    max_voxel_dist2center = np.sqrt(np.max(voxel_dist2center_squared))

    return average_voxel_dist2center, max_voxel_dist2center

