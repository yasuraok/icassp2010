# -*- coding: utf-8 -*-
"""
$Date::                           $
$Rev::                            $
$Author::                         $

gauss窓を追加
"""

import numpy as np
from scipy import hanning, bartlett, blackman, hamming#, kaiser

def gaussian(length, var=None):
    if var == None:
        halfptr = length / 4.
        # solve exp(- (halfptr**2 / 2 * var)) = 0.5 for variance
        var = -halfptr ** 2 / (2 * np.log(0.5))
    ptr = np.arange((1 - length) / 2., (length + 1) / 2.)
    Y = np.exp(-ptr ** 2 / (2 * var))
    #Y = ((len+1)/2.) * (Y / sum(Y));
    return Y


if __name__ == '__main__':
    from pylab import *
    fftLen = 1024
    std = fftLen / 2
        
    Y1 = gaussian(fftLen)
    Y2 = hanning(fftLen)
    Y3 = hamming(fftLen)
    Y4 = blackman(fftLen)
    Y5 = bartlett(fftLen)
#    Y6 = kaiser(fftLen)
    plot(Y1, label="gaussian")
    plot(Y2, label="hanning")
    plot(Y3, label="hamming")
    plot(Y4, label="blackman")
    plot(Y5, label="bartlett")
#    plot(Y6, label="kaiser")
    ylim(0, 1)
    legend()
    show()
