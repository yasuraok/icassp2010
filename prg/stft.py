# -*- coding: utf-8 -*-
"""
$Date::                           $
$Rev::                            $
$Author::                         $

至って普通のスペクトログラム算出
"""

from numpy import empty, linspace, r_, zeros, real, conj, hstack, vstack
from scipy.fftpack import fft, ifft

def sec2frm(sec, fftlength, fs, overlap):
    hop = fftlength - overlap
    tmp = sec * fs - fftlength
    return tmp / hop

def sgn2frm(sgnLen, fftlength, overlap):
    hop = fftlength - overlap
    tmp = sgnLen - fftlength
    frmLen = ((tmp / hop + 1) if tmp >= 0 else 1)
    return frmLen

def frm2sgn(frmLen, fftlength, overlap):
    hop = fftlength - overlap
    frmLen = frmLen
    sgnLen = hop * (frmLen - 1) + fftlength
    return sgnLen

#def specgram(i_x, fftlength, fs, window, overlap):
#    # X is assumed to be complex spectrogram obtained by function "specgram"
#    if i_x.ndim == 1:
#        x = i_x.reshape(len(i_x), 1)
#    elif i_x.ndim == 2:
#        x = i_x
#    else:
#        assert False
#    
#    hop    = fftlength - overlap
#    tmp    = x.shape[0] - fftlength
#    #frmLen = ((tmp/hop + 1) + int(bool(tmp%hop)) if tmp >= 0 else 1)
#    frmLen = ((tmp/hop + 1) if tmp >= 0 else 1)
#    sgnLen = hop * (frmLen-1) + fftlength
#    binLen = fftlength / 2.0 + 1
#    
#    #x = hstack([x, tile(0.0, (sgnLen - x.size))])
#    Y = empty((binLen, frmLen), dtype=complex)
#    for f in range(0, frmLen):
#        offset = f*hop
#        Y[:,f] = fft(x[offset:offset+fftlength,:].T * window)[0,0:binLen]
#    
#    F = linspace(0, fs/2, binLen)
#    T = r_[0 : float(tmp)/fs : float(hop)/fs]
#    return Y, F, T

def specgram(i_x, fftlength, fs, window, overlap):
    # X is assumed to be complex spectrogram obtained by function "specgram"
    if i_x.ndim == 1:
        x = i_x.reshape(len(i_x), 1)
    elif i_x.ndim == 2:
        x = i_x
    else:
        assert False
    
    hop = fftlength - overlap
    tmp = x.shape[0] - fftlength
    #frmLen = ((tmp/hop + 1) + int(bool(tmp%hop)) if tmp >= 0 else 1)
    frmLen = ((tmp / hop + 1) if tmp >= 0 else 1)
    sgnLen = hop * (frmLen - 1) + fftlength
    binLen = fftlength / 2 + 1
    
    #x = hstack([x, tile(0.0, (sgnLen - x.size))])
    Y = empty((frmLen, binLen), dtype=complex)
    for f in range(0, frmLen):
        offset = f * hop
        Y[f, :] = fft(x[offset:offset + fftlength, :].T * window)[0, 0:binLen]
    
    F = linspace(0, fs / 2, binLen)
    T = r_[0 : float(tmp) / fs : float(hop) / fs]
    return Y, F, T

def overlapdecomp(i_x, fftlength, fs, window, overlap):
    if i_x.ndim == 1:
        x = i_x.reshape(len(i_x), 1)
    elif i_x.ndim == 2:
        x = i_x
    else:
        assert False
    
    hop = fftlength - overlap
    tmp = x.shape[0] - fftlength
    frmLen = ((tmp / hop + 1) if tmp >= 0 else 1)
    #sgnLen = hop * (frmLen-1) + fftlength
    #x = vstack([x, zeros((sgnLen-x.shape[0], x.shape[1]), dtype=x.dtype)])
    
    Y = empty((frmLen, fftlength), dtype=x.dtype)
    for f in range(0, frmLen):
        offset = f * hop
        #print offset+fftlength, x.shape[0]
        Y[f, :] = x[offset:offset + fftlength, 0] * window
    
    return Y
    

def overlapsynth(X, fftlength, fs, window, overlap):
    hop = fftlength - overlap
    frmLen = X.shape[0]
    sgnLen = hop * (frmLen - 1) + fftlength
    
    y = zeros((sgnLen), dtype=float)
    w2 = zeros((sgnLen), dtype=float)
    window2 = window ** 2
    for f in xrange(0, frmLen):
        offset = f * hop
        w2[offset:offset + fftlength] += window2
    (w2[0], w2[-1]) = (1, 1)
        
    for f in range(0, frmLen):
        offset = f * hop
        y[offset:offset + fftlength] += X[f, :] * window
#    print "1", y.shape
    y = y / w2
#    print "2", y.shape
#    print "3", y.reshape(sgnLen, 1).shape
    return y.reshape(sgnLen, 1)


def ispecgram(X, fftlength, fs, window, overlap):
    # X is assumed to be complex spectrogram obtained by function "specgram"
    # X[frmLen, binLen]
    
    hop = fftlength - overlap
    frmLen = X.shape[0]
    sgnLen = hop * (frmLen - 1) + fftlength
    
    [y, w2] = [zeros((sgnLen, 1)) for i in xrange(2)]
    window2 = window ** 2
    for f in xrange(0, frmLen):
        offset = f * hop
        w2[offset:offset + fftlength, 0] += window2
    (w2[0:fftlength], w2[-fftlength:-1]) = (1, 1)
        
    for f in range(0, frmLen):
        Xf = X[f, :]
        offset = f * hop
        y[offset:offset + fftlength, 0] += real(ifft(hstack([Xf, conj(Xf[-2:0:-1])]))) * window
        
    #for n in xrange(2000):
    #    print n, y[n], (y / w2)[n]
    
    y = y / w2
    #y = y * (hop / sum(window))
    
    return y

def ispecgram_simple(X, fftlength, fs, window, overlap):
    # X is assumed to be complex spectrogram obtained by function "specgram"
    
    hop = fftlength - overlap
    frmLen = X.shape[1]
    sgnLen = hop * (frmLen - 1) + fftlength
    
    [y, w] = [zeros((sgnLen, 1)) for i in xrange(2)]
    for f in xrange(0, frmLen):
        offset = f * hop
        w[offset:offset + fftlength, 0] += window
    (w[0:fftlength], w[-fftlength:-1]) = (1, 1)
    
    for f in range(0, frmLen):
        Xf = X[:, f]
        offset = f * hop
        y[offset:offset + fftlength, 0] += real(ifft(hstack([Xf, conj(Xf[-2:0:-1])])))
    
    y = y / w
    #y = y * (hop / sum(window))
    
    return y
