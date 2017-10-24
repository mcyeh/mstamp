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

if __name__ == '__main__':
    mat = sio.loadmat('toy_data.mat')
    data = mat['data']
    sub_len = mat['sub_len'][0][0]

    # using the stomp based method to compute the multidimensional matrix
    # profile
    [mat_pro_1, pro_idx_1, pro_dim_1] = mstamp_stomp(data.T, sub_len)

    # plot the matrix profile as image
    plt.figure()
    plt.imshow(mat_pro_1, extent=[0, 1, 0, 1])

    # using the stamp based method to compute the multidimensional matrix
    # profile
    [mat_pro_2, pro_idx_2, pro_dim_2] = mstamp_stamp(data.T, sub_len)

    # plot the matrix profile as image
    plt.figure()
    plt.imshow(mat_pro_2, extent=[0, 1, 0, 1])
