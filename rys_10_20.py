# -*- coding: utf-8 -*-
from __future__ import division, print_function
import matplotlib.pyplot as plt
import scipy.signal as ss
import numpy as np


def rys(tytul =  ''): 
    lab =  ['Fp1', 'Fpz', 'Fp2','F7', 'F3', 'Fz', 'F4','F8', 'T3' ,'C3', 'Cz', 'C4', 'T4', 'T5', 'P3', 'Pz','P4', 'T6', 'O1', 'O2', 'Oz', 'legend']
    #plt.ion()
    f = plt.figure()
    plt.suptitle(tytul, fontsize = 15)
    axis =  []
    index = 0
    axis.append(f.add_subplot(5, 5, 2))
    plt.title(lab[0])
    axis.append(f.add_subplot(5, 5, 3))
    plt.title(lab[1])
    axis.append(f.add_subplot(5, 5, 4))
    plt.title(lab[2])
    axis.append(f.add_subplot(5, 5, 6))
    plt.title(lab[3])
    axis.append(f.add_subplot(5, 5, 7))
    plt.title(lab[4])
    axis.append(f.add_subplot(5, 5, 8))
    plt.title(lab[5])
    axis.append(f.add_subplot(5, 5, 9))
    plt.title(lab[6])
    axis.append(f.add_subplot(5, 5, 10))
    plt.title(lab[7])
    axis.append(f.add_subplot(5, 5, 11))
    plt.title(lab[8])
    axis.append(f.add_subplot(5, 5, 12))
    plt.title(lab[9])
    axis.append(f.add_subplot(5, 5, 13))
    plt.title(lab[10])
    axis.append(f.add_subplot(5, 5, 14))
    plt.title(lab[11])
    axis.append(f.add_subplot(5, 5, 15))
    plt.title(lab[12])
    axis.append(f.add_subplot(5, 5, 16))
    plt.title(lab[13])
    axis.append(f.add_subplot(5, 5, 17))
    plt.title(lab[14])
    axis.append(f.add_subplot(5, 5, 18))
    plt.title(lab[15])
    axis.append(f.add_subplot(5, 5, 19))
    plt.title(lab[16])
    axis.append(f.add_subplot(5, 5, 20))
    plt.title(lab[17])
    axis.append(f.add_subplot(5, 5, 22))
    plt.title(lab[18])
    axis.append(f.add_subplot(5, 5, 24))
    plt.title(lab[19])
    axis.append(f.add_subplot(5, 5, 23))
    plt.title(lab[20])
    axis.append(f.add_subplot(5, 5, 21))
    axis[-1].axis('off')
    return {lab_: axis[i] for i, lab_ in enumerate(lab)}
