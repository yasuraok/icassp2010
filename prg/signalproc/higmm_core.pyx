# -*- coding: utf-8 -*-
"""
$Date::                           $
$Rev::                            $
$Author::                         $

"""
import numpy as np

cimport cython
cimport numpy as np; np.import_array(); np.import_ufunc()

cdef double nfloor = 1e-06 # -120dB noise floor
cdef inline float float_max(float a, float b) nogil:
    return a if a >= b else b

import pylab
from PyQt4.QtGui import QApplication as QA


cpdef estim1frame(freqs, f0f, ampf, var, AXf, param):
    cdef unsigned int i
    for i in xrange(10):
        gmm, distIdx = calcHGMM(freqs, f0f, var, param.numH)
        model   = applyAmp2HGMM(gmm, ampf, distIdx)
#        jensen1  = AXf[:] * model[:,:] / model.sum(0)[None,:]
        jensen  = (AXf[None,:] / model.sum(0)[None,:]) * model[:,:]
#        print jensen.sum(0)
#        print AXf
#        
#        
#        modelsum = model.sum(0)
#        print (AXf * np.log(AXf / modelsum) - (AXf - modelsum)).sum()
##        print distIdx
#        
#        pylab.plot(modelsum[:300])
#        pylab.ylim(0, 15)
#        pylab.hold(True)
##        pylab.plot(model[1,:300])
##        pylab.plot(distIdx[0,:300] / 2.)
##        pylab.plot(jensen[1,:300])
#        pylab.plot(AXf[:300])
#        pylab.hold(False)
#        QA.processEvents()
#        raw_input()
        
        # Update amplitude of each harmonic peak
        ampf[:] = sumupAmp(ampf[:], distIdx, jensen, model, param.numH)
        # Update fundamental frequency
        f0f[:]  = updateF0(f0f, freqs, distIdx, jensen)
        
#        print f0f


cpdef synth1frame(freqs, f0f, ampf, var, param):
    gmm, distIdx = calcHGMM(freqs, f0f, var, param.numH)
    model        = applyAmp2HGMM(gmm, ampf, distIdx)
    return model.sum(0)


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calcHGMM(
        np.ndarray[np.float_t, ndim=1, mode="c"] freqs,
        np.ndarray[np.float_t, ndim=1, mode="c"] f0f,
        np.ndarray[np.float_t, ndim=1, mode="c"] var,
        np.uint_t numH):
    """
    Calculate harmonic gaussian array
    """
    
    distIdx = np.minimum((freqs[None,:] / f0f[:,None] - 0.5).astype(np.uint) + 1, numH - 1) #[S,B]
    mean    = distIdx * f0f[:,None]
    gmm     = nfloor + np.exp(-(freqs[None,:] - mean)**2 / (2 * var[:,None])) #[S,B]
    return gmm, distIdx


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef applyAmp2HGMM(
        np.ndarray[np.float_t, ndim=2, mode="c"] gmm,      #[S,B]
        np.ndarray[np.float_t, ndim=2, mode="c"] amp,      #[S,H]
        np.ndarray[np.uint_t,  ndim=2, mode="c"] distIdx): #[S,B]
    """
    Multiply amplitude with harmonic gmm
    """

    cdef unsigned int s, b, S, H, B
    (S, H, B) = gmm.shape[0], amp.shape[1], gmm.shape[1]
    cdef np.ndarray[np.float_t, ndim=2, mode="c"] ret = np.empty((S,B), dtype=float)
    for s in xrange(S):
        for b in xrange(B):
            ret[s,b] = gmm[s,b] * amp[s,distIdx[s,b]]
    return ret


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sumupAmp(
        np.ndarray[np.float_t, ndim=2, mode="c"] amp,
        np.ndarray[np.uint_t,  ndim=2, mode="c"] distIdx,
        np.ndarray[np.float_t, ndim=2, mode="c"] value,
        np.ndarray[np.float_t, ndim=2, mode="c"] model,
        unsigned int H):
    """
    """

    cdef unsigned int s, b, h, S, B
    S, B = value.shape[0], value.shape[1]
    cdef np.ndarray[np.float_t, ndim=1, mode="c"] nume = np.zeros((H), dtype=float)
    cdef np.ndarray[np.float_t, ndim=1, mode="c"] deno = np.zeros((H), dtype=float)
    
    for s in xrange(S):
        nume[:] = deno[:] = nfloor
        for b in xrange(B):
            h = distIdx[s,b]
            nume[h] += value[s,b]
            deno[h] += model[s,b]
        for h in xrange(H):
            amp[s,h] *= nume[h] / deno[h]
    return amp


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef updateF0(
        np.ndarray[np.float_t, ndim=1, mode="c"] F0f,     # [S]
        np.ndarray[np.float_t, ndim=1, mode="c"] freqs,   # [B]
        np.ndarray[np.uint_t,  ndim=2, mode="c"] distIdx, # [S,B]
        np.ndarray[np.float_t, ndim=2, mode="c"] value):  # [S,B]
    """
    """
    cdef unsigned int s, b, S, B
    cdef double nume, deno
    S, B = value.shape[0], value.shape[1]

    for s in xrange(S):
        nume = deno = nfloor
        for b in xrange(B):
            nume  += value[s,b] * distIdx[s,b] * freqs[b]
            deno  += value[s,b] * (distIdx[s,b] * distIdx[s,b])
#            vnume += AUXHk * bias[s,b]
#            vdeno += AUXHk
        F0f[s] = nume / deno
    return F0f
