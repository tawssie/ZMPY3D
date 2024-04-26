# GPL3.0, JSL from ZM
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

def calculate_bbox_moment_2_zm05(max_order, g_cache_complex, g_cache_pqr_linear, g_cache_complex_index, clm_cache3d, bbox_moment):

    def complex_nan():
        c = np.zeros(1, dtype=np.complex128)
        c[0] = np.nan
        return c

    max_n = max_order + 1
    
    bbox_moment = np.reshape(np.transpose(bbox_moment, (2, 1, 0)), -1)
    
    zm_geo = g_cache_complex * bbox_moment[g_cache_pqr_linear - 1] 
    
    zm_geo_sum = np.zeros(max_n * max_n * max_n, dtype=np.complex128)
    
    np.add.at(zm_geo_sum, g_cache_complex_index - 1, zm_geo)


    zm_geo_sum[zm_geo_sum == 0.0] = complex_nan()


    z_moment_raw = zm_geo_sum * (3.0 / (4.0 * np.pi))
    z_moment_raw = z_moment_raw.reshape((max_n, max_n, max_n))
    z_moment_raw = np.transpose(z_moment_raw, (2, 1, 0))
    z_moment_scaled = z_moment_raw * clm_cache3d  # CLMCache3D is a 3D matrix, so operations are direct

    return z_moment_scaled, z_moment_raw


    
