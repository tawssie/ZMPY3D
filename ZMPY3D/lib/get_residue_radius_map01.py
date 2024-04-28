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

import math


def get_residue_radius_map01():
    coef = math.sqrt(5.0 / 3.0)

    residue_radii = {
        'ALA': 1.963913 * coef,
        'ARG': 3.374007 * coef,
        'ASN': 2.695111 * coef,
        'ASP': 2.525241 * coef,
        'CYS': 2.413249 * coef,
        'GLN': 3.088783 * coef,
        'GLU': 2.883527 * coef,
        'GLY': 1.841949 * coef,
        'HIS': 2.652737 * coef,
        'ILE': 2.575828 * coef,
        'LEU': 2.736953 * coef,
        'LYS': 3.177825 * coef,
        'MET': 2.959014 * coef,
        'MSE': 2.959014 * coef,
        'PHE': 2.979213 * coef,
        'PRO': 2.266054 * coef,
        'SER': 2.184637 * coef,
        'THR': 2.366486 * coef,
        'TRP': 3.248871 * coef,
        'TYR': 3.217711 * coef,
        'VAL': 2.351359 * coef,
        'A': 4.333750 * coef,
        'T': 3.700942 * coef,
        'G': 4.443546 * coef,
        'C': 3.954067 * coef,
        'U': 3.964129 * coef,
        'I': 4.0 * coef,
        'DA': 4.333750 * coef,
        'DT': 3.700942 * coef,
        'DG': 4.443546 * coef,
        'DC': 3.954067 * coef,
        'DU': 3.964129 * coef,
        'DI': 4.0 * coef
    }

    return residue_radii

