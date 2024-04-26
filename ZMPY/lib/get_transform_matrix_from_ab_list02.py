# MIT License, BIOZ
#
# Copyright (c) 2024 Jhih-Siang Lai
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import numpy as np

def get_transform_matrix_from_ab_list02(a, b, center_scaled):
    a2pb2 = a ** 2 + b ** 2
    a2mb2 = a ** 2 - b ** 2

    m33_linear = np.array([
        np.real(a2pb2),
        -np.imag(a2mb2),
        2 * np.imag(a * b),
        np.imag(a2pb2),
        np.real(a2mb2),
        -2 * np.real(a * b),
        2 * np.imag(a * np.conj(b)),
        2 * np.real(a * np.conj(b)),
        np.real(a * np.conj(a)) - np.real(b * np.conj(b))
    ])
    
    scale = 1.0
    
    m33 = m33_linear.reshape([3, 3])
    m44 = np.zeros((4, 4))
    m44[0:3, 0:3] = m33 * scale
    m44[0:3, 3] = center_scaled.flatten()
    m44[3, 3] = 1
    transform = np.linalg.inv(m44)
    
    return transform



