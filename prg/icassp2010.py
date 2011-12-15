# -*- coding: utf-8 -*-
"""
$Date:: 2012-03-15 00:34:05 +0900#$
$Rev:: 147                        $
$Author:: yasuraoka               $

Example of ICASSP 2010
"""

import numpy as np
import os.path as osp
import pylab
from PyQt4.QtGui import QApplication as QA

from audio import wavio
from signalproc import higmm, dereverberate
from stft import specgram, ispecgram

def exec_estim(conf, wavpath):
    pylab.ion()
    fig = pylab.figure()
    fig.hold(False)
#    QA.processEvents()


    fftlength = conf['fftLen']
    fs = conf['fs']
    overlap = conf['overlap']
    gaussian = conf['window'](fftlength)
    
    # STFT analysis
    (x, fs, bits) = wavio.audioreadmono(wavpath)
    (X, F, T) = specgram(x[:882000], fftlength, fs, gaussian, overlap)
#    X[:, :5] *= conf["nfloor"] # kill 50Hz or 60Hz noise
    PX = np.abs(X)**2
    Y  = X.copy()

    F, B = X.shape
    srcEstimator = higmm.HIGMM()
    srcEstimator.initParam(F, B, conf, higmm.HIGMMParam(numS = 8, numH = 150))
    revEstimator = dereverberate.EstimRev()
    G = revEstimator.initParam(F, B, conf, dereverberate.DerevParam(p_step = 8, derevorder = 80))
    
    
    # observation
    ax  = fig.add_subplot(2,2,1)
    ax.hold(False)
    ax.imshow(10 * np.log10(PX).T, None, None, 'auto', 'nearest', None, -70, 40, 'lower')
    ax.figure.canvas.draw()
    QA.processEvents()

    for i in xrange(3):
        srcEstimator.estim(np.abs(Y)**2)
        PS = srcEstimator.synth()
        
        # source model spectra
        ax  = fig.add_subplot(2,2,2)
        ax.imshow(10 * np.log10(abs(PS)).T, None, None, 'auto', 'nearest', None, -70, 40, 'lower')
        ax.figure.canvas.draw()
        QA.processEvents()
        
        G, Y = revEstimator.estimAndDerev(X, PS, G)
    
        # reverberation 
        ax  = fig.add_subplot(2,2,3)
        ax.imshow(20 * np.log10(abs(X - Y)).T, None, None, 'auto', 'nearest', None, -70, 40, 'lower')
        ax.figure.canvas.draw()
        QA.processEvents()

        # dereverberated
        ax  = fig.add_subplot(2,2,4)
        ax.imshow(20 * np.log10(abs(Y)).T, None, None, 'auto', 'nearest', None, -70, 40, 'lower')
        ax.figure.canvas.draw()
        QA.processEvents()


    y = ispecgram(Y, fftlength, fs, gaussian, overlap)
    outpath = osp.splitext(wavpath)[0]
    wavio.audiowrite(y, fs, bits, outpath + "_derev.wav")

    Yw = revEstimator.derev_wiener(X, G, Y)
    y = ispecgram(Yw, fftlength, fs, gaussian, overlap)
    wavio.audiowrite(y, fs, bits, outpath + "_derev_wiener.wav")

    print "Completed. Press any key to exit."
    raw_input()
    return


if __name__ == "__main__":
    from window import gaussian

    conf = { \
      'fs'            : 44100,
      'bits'          : 16,
      'fftLen'        : 2048,
      'overlap'       : 2048 - 512, # fftLen / shift must be integer
      'window'        : gaussian,
      'srcreviter'    : 3,
      'srciter'       : 1,
      'nfloor'        : 10 ** (-120 / 20.), # -120dB
    }
    
    import argparse    
    parser = argparse.ArgumentParser(description='ICASSP 2010 smallest example')
    parser.add_argument('wavpath', help="Input wave file path")
    parser.add_argument('--version', action='version', version='%(prog)s 0.01') # version
    args = parser.parse_args()
        
    exec_estim(conf, args.wavpath)


