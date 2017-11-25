# -*- coding: utf-8 -*-
"""
@author: Michael Yeh

C.-C. M. Yeh, N. Kavantzas, and E. Keogh, "Matrix Profile VI: Meaningful
Multidimensional Motif Discovery," IEEE ICDM 2017.
https://sites.google.com/view/mstamp/
http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
"""

import scipy.io as sio
import matplotlib.pyplot as plt
from mstamp_stomp import mstamp as mstamp_stomp
from mstamp_stamp import mstamp as mstamp_stamp


def plot_motifs(matrix_profile, dimensionality=1):
    motif_at = matrix_profile[dimensionality - 1, :].argsort()[:2]

    plt.figure(figsize=(14, 7))
    for i in range(3):
        plt.subplot(4, 1, i + 1)
        plt.plot(data.T[i, :])
        plt.title('$T_{}$'.format(i + 1))
        for m in motif_at:
            plt.plot(range(m, m + sub_len), data.T[i, :][m:m + sub_len], c='r')
        plt.xlim((0, matrix_profile.shape[1]))

    plt.subplot(414)
    plt.title('{}-dimensional Matrix Profile'.format(dimensionality))
    plt.plot(matrix_profile[dimensionality - 1, :])
    for m in motif_at:
        plt.axvline(m, c='r')
    plt.xlim((0, matrix_profile.shape[1]))
    plt.tight_layout()


if __name__ == '__main__':
    mat = sio.loadmat('toy_data.mat')
    data = mat['data']
    sub_len = mat['sub_len'][0][0]

    # using the stomp based method to compute the multidimensional matrix
    # profile
    mat_pro_1, pro_idx_1 = mstamp_stomp(data.T, sub_len,
                                        return_dimension=False)

    # plot the matrix profile as image
    plt.figure()
    plt.title('Matrix Profile (STOMP)')
    plt.imshow(mat_pro_1, extent=[0, 1, 0, 1])

    # using the stamp based method to compute the multidimensional matrix
    # profile
    mat_pro_2, pro_idx_2 = mstamp_stamp(data.T, sub_len,
                                        return_dimension=False)

    # plot the matrix profile as image
    plt.figure()
    plt.title('Matrix Profile (STAMP)')
    plt.imshow(mat_pro_2, extent=[0, 1, 0, 1])

    plot_motifs(mat_pro_2)

    plt.show()
