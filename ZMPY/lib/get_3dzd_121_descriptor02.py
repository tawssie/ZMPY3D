# GPL 3.0 License, JSL from BIOZ
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

def get_3dzd_121_descriptor02(z_moment_scaled):
    z_moment_scaled[np.isnan(z_moment_scaled)] = 0
    z_moment_scaled_norm = np.abs(z_moment_scaled) ** 2

    z_moment_scaled_norm_positive = np.sum(z_moment_scaled_norm, axis=2)

    z_moment_scaled_norm[:, :, 0] = 0
    z_moment_scaled_norm_negative = np.sum(z_moment_scaled_norm, axis=2)

    zm_3dzd_invariant = np.sqrt(z_moment_scaled_norm_positive + z_moment_scaled_norm_negative)
    zm_3dzd_invariant[zm_3dzd_invariant < 1e-20] = np.nan

    return zm_3dzd_invariant


