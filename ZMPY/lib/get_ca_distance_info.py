# GPL3.0 License, JSLai
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

def get_ca_distance_info(xyz):

    xyz_center = np.mean(xyz, axis=0)
    xyz_dist2center = np.sqrt(np.sum((xyz - xyz_center) ** 2, axis=1))

    percentiles_for_geom = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
    percentile_list = np.percentile(xyz_dist2center, percentiles_for_geom)
    percentile_list = percentile_list.reshape(-1, 1)

    std_xyz_dist2center = np.std(xyz_dist2center, ddof=1)


    n = len(xyz_dist2center)
    mean_distance = np.mean(xyz_dist2center)
    s = (n / ((n - 1) * (n - 2))) * np.sum(((xyz_dist2center - mean_distance) / std_xyz_dist2center) ** 3)

    fourth_moment = np.sum(((xyz_dist2center - mean_distance) / std_xyz_dist2center) ** 4)
    k = (n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * fourth_moment -
         3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))

    return percentile_list, std_xyz_dist2center, s, k


