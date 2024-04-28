# GPL3.0, JSL, ZM, bioz
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


from .eigen_root import *
import numpy as np

def calculate_ab_rotation_02(z_moment_raw, target_order2_norm_rotate):
    ind_real = 2

    if target_order2_norm_rotate % 2 == 0:
        abconj_coef = [
            z_moment_raw[target_order2_norm_rotate, 2, 2],
            -z_moment_raw[target_order2_norm_rotate, 2, 1],
            z_moment_raw[target_order2_norm_rotate, 2, 0],
            np.conj(z_moment_raw[target_order2_norm_rotate, 2, 1]),
            np.conj(z_moment_raw[target_order2_norm_rotate, 2, 2])
        ]
        abconj_sol = eigen_root(abconj_coef)
        n_abconj = 4
    else:
        abconj_coef = [
            z_moment_raw[target_order2_norm_rotate, 1, 1],
            -z_moment_raw[target_order2_norm_rotate, 1, 0],
            -np.conj(z_moment_raw[target_order2_norm_rotate, 1, 1])
        ]
        abconj_sol = eigen_root(abconj_coef)
        n_abconj = 2

    k_re = np.real(abconj_sol)
    k_im = np.imag(abconj_sol)
    k_im2 = k_im ** 2
    k_re2 = k_re ** 2
    k_im3 = k_im ** 3
    k_im4 = k_im ** 4
    k_re4 = k_re ** 4

    f20 = np.real(z_moment_raw[ind_real, 2, 0])
    f21 = z_moment_raw[ind_real, 2, 1]
    f22 = z_moment_raw[ind_real, 2, 2]

    f21_im = np.imag(f21)
    f21_re = np.real(f21)
    f22_im = np.imag(f22)
    f22_re = np.real(f22)

    coef4 = 4 * f22_re * k_im * (-1 + k_im2 - 3 * k_re2) - 4 * f22_im * k_re * (1 - 3 * k_im2 + k_re2) - 2 * f21_re * k_im * k_re * (-3 + k_im2 + k_re2) + 2 * f20 * k_im * (-1 + k_im2 + k_re2) + f21_im * (1 - 6 * k_im2 + k_im2 ** 2 - k_re2 ** 2)
    coef3 = 2 * (-4 * f22_im * (k_im + k_im3 - 3 * k_im * k_re2) + f21_re * (-1 + k_im4 + 6 * k_re2 - k_re4) + 2 * k_re * (f22_re * (2 + 6 * k_im2 - 2 * k_re2) + f21_im * k_im * (-3 + k_im2 + k_re2) + f20 * (-1 + k_im2 + k_re2)))

    bimbre_coef = np.array([coef4, coef3, np.zeros_like(coef4), coef3, -coef4]).T

    bimbre_sol_real = [np.real(eigen_root(bc)) for bc in bimbre_coef]

    is_abs_bimre_good = [np.abs(x) > 1e-7 for x in bimbre_sol_real]

    a_list = []
    b_list = []

    for i in range(len(bimbre_sol_real)):
        bre = 1 / np.sqrt((1 + k_im2[i] + k_re2[i]) * (1 + np.power(bimbre_sol_real[i], 2)))
        bim = bimbre_sol_real[i] * bre
        b = np.vectorize(complex)(bre, bim)
        a = abconj_sol[i] * np.conj(b)

        a_list.append(a[is_abs_bimre_good[i]])
        b_list.append(b[is_abs_bimre_good[i]])

    ab_list = np.concatenate((np.concatenate(a_list).reshape(-1, 1), np.concatenate(b_list).reshape(-1, 1)), axis=1)

    return ab_list




