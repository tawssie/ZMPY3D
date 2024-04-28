# GPL3.0 License, JSL, ZM
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

def calculate_zm_by_ab_rotation01(z_moment_raw, binomial_cache, ab_list, max_order, clm_cache, s_id, n, l, m, mu, k, is_nlm_value):
    zm_rotated_list = [None] * len(ab_list)

    a = ab_list[:, 0]
    b = ab_list[:, 1]
    a = a.flatten()
    b = b.flatten()

    aac = np.real(a * np.conj(a)).astype(np.complex128)
    bbc = np.real(b * np.conj(b)).astype(np.complex128)
    bbcaac = -bbc / aac

    abc = -(a / np.conj(b))
    ab = a / b

    bbcaac_pow_k_list = np.log(bbcaac)[:, None] * np.arange(max_order + 1)
    aac_pow_l_list = np.log(aac)[:, None] * np.arange(max_order + 1)
    ab_pow_m_list = np.log(ab)[:, None] * np.arange(max_order + 1)
    abc_pow_mu_list = np.log(abc)[:, None] * np.arange(-max_order, max_order + 1)

    f_exp = np.zeros(len(s_id), dtype=np.complex128)
    f_exp[mu >= 0] = z_moment_raw[n[mu >= 0], l[mu >= 0], mu[mu >= 0]]
    f_exp[(mu < 0) & (mu % 2 == 0)] = np.conj(z_moment_raw[n[(mu < 0) & (mu % 2 == 0)], l[(mu < 0) & (mu % 2 == 0)], -mu[(mu < 0) & (mu % 2 == 0)]])
    f_exp[(mu < 0) & (mu % 2 != 0)] = -np.conj(z_moment_raw[n[(mu < 0) & (mu % 2 != 0)], l[(mu < 0) & (mu % 2 != 0)], -mu[(mu < 0) & (mu % 2 != 0)]])

    f_exp = np.log(f_exp)

    max_n = max_order + 1
    clm = clm_cache[l * max_n + m].astype(np.complex128)
    clm = clm.flatten()

    bin = binomial_cache[l - mu, k - mu].astype(np.complex128) + binomial_cache[l + mu, k - m].astype(np.complex128)

    for zm_i in range(len(zm_rotated_list)):
        al = aac_pow_l_list[zm_i, l]
        al = al.flatten()

        abpm = ab_pow_m_list[zm_i, m]
        abpm = abpm.flatten()

        amu = abc_pow_mu_list[zm_i, max_order + mu]
        amu = amu.flatten()

        bbk = bbcaac_pow_k_list[zm_i, k]
        bbk = bbk.flatten()

        nlm = f_exp + al + clm + abpm + amu + bbk + bin

        z_nlm = np.zeros(is_nlm_value.shape, dtype=np.complex128)
        np.add.at(z_nlm, s_id, np.exp(nlm))

        zm = np.full((np.prod(z_moment_raw.shape), ), np.nan, dtype=np.complex128)
        zm[is_nlm_value] = z_nlm
        zm = zm.reshape(z_moment_raw.shape)
        zm = np.transpose(zm, (2, 1, 0))

        zm_rotated_list[zm_i] = zm

    return zm_rotated_list

