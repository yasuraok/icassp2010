# -*- coding: utf-8 -*-
"""
$Date::                           $
$Rev::                            $
$Author::                         $

Blind dereververation based on complex linear inverse filtering
and IS-divergence minimization
http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?arnumber=5496223
"""

import sys, numpy as np
from scipy.linalg import solve
from scipy.fftpack import fft, ifft
from fftfilt import fftfilt

################################################################################

class DerevParam(object):
    __slots__ = (
        "p_step",
        "derevorder",
    )
    
    __defaults__ = {
        "p_step": 10,
        "derevorder": 80,
    }
    
    def __init__(self, **kwargs):
        for k, v in self.__defaults__.items():
            self.__setattr__(k, kwargs.get(k, v))


class EstimRev(object):
    __slots__ = (
        "fltLen",
        "zfltLen", # array containing zero for zero padding
        "param",
    )

    def initParam(self, F, B, conf, param = None):
        self.param   = param if param != None else DerevParam()
        self.fltLen  = self.param.p_step + self.param.derevorder
        self.zfltLen = np.zeros(self.fltLen, dtype=complex)

        # initialize 
        G = np.zeros((self.fltLen, B), dtype=complex)
        G[0, :] = 1
        return G
        
    def estimAndDerev(self, X, srcP, G):
        p_step, derevorder, fltLen = self.param.p_step, self.param.derevorder, self.fltLen
        nfloor = 10 ** ((20 * np.log10(srcP.max()) - 150) / 20)
        (F, B) = X.shape
        
        isrcP = (1.0 / np.maximum(srcP, nfloor))
        Y = np.empty((F, B), dtype=complex) # dereverberated signal
        
        print "Estimating reverb filter".ljust(40),
                
        for b in xrange(B):
            Xb = X[:, b]
            iSb = isrcP[:, b]
            
            R, r = self.calcRr(Xb, iSb, Xb)
            G[p_step:, b] = -solve(R, r, sym_pos=True, overwrite_a=True, overwrite_b=True)

            if int(b / (B/10.)) != int((b+1) / (B/10.)):
                print "*",
                sys.stdout.flush()
            
            Y[:, b] = fftfilt(G[:, b], Xb)[0:F]
        print
            
        return G, Y

    def calcRr(self, _Xb, _iSb, _Yb):
        """
        Calculate modified correlation matrix R and r using FFT
        Xb:  observed signal for calculating R
        iSb: inverce source power spectrum
        Xb:  observed signal for calculating r
        """
        p_step, derevorder, fltLen = self.param.p_step, self.param.derevorder, self.fltLen
        
        Xb  = np.r_[self.zfltLen, _Xb, self.zfltLen]
        iSb = np.r_[_iSb, self.zfltLen]
        Yb  = np.r_[_Yb, self.zfltLen]

        fftXbc = fft(Xb[fltLen:]).conj()
        r = ifft(fft(Yb * iSb, overwrite_x=True) * fftXbc, overwrite_x=True)[p_step:fltLen]
        R = np.empty((derevorder, derevorder), dtype=complex)
        for i in xrange(p_step, fltLen):
            fftXbiSb = fft(Xb[fltLen - i:-i] * iSb, overwrite_x=True)
            R[:, i - p_step] = ifft(fftXbiSb * fftXbc, overwrite_x=True)[p_step:fltLen]
        
        return R, r
    
    
    def derev_wiener(self, X, G, Y=None):
        """
        Power-domain filtering based on wiener filter
        This suppresses reverberation stronger but makes non-linear noises.
        """
        if type(Y) == type(None):
            (F, B) = X.shape
            Y = np.empty((F, B), dtype=complex) # dereverberated signal
            for b in range(B):
                Y[:, b] = fftfilt(G[:, b], X[:, b])[0:F]
        filt = (Y.real ** 2 + Y.imag ** 2) / (X.real ** 2 + X.imag ** 2)
        return filt * X
    

#    def calcRr(self, Xb, iSb, Yb):
#        derevorder, p_step, revorder, fltLen = self.derevorder, self.p_step, self.revorder, self.fltLen
#        F = len(Xb)
#        
#        z   = np.zeros(fltLen, dtype=complex)
#        tmp = np.zeros(F+fltLen, dtype=complex)
#        ii  = r_[0:derevorder]
#        R = empty((derevorder, derevorder), dtype=complex)
#        r = empty(derevorder, dtype=complex)
#    
#        conjX = Xb.conj();
#        tmp1 = fft(r_[iSb, z])
#        for j in xrange(0, p_step):
#            tmp[0:j] = 0
#            tmp[j:F] = conjX[j:] * Xb[0:F-j]
#            tmp2 = ifft(tmp1 * fft(tmp).conj())
#            i = ii[:derevorder-j]
#            R[i, i+j] = tmp2[i+p_step]
#        for j in xrange(p_step, derevorder):
#            tmp[0:j] = 0
#            tmp[j:F] = conjX[j:] * Xb[0:F-j]
#            tmp2 = ifft(tmp1 * fft(tmp).conj())
#            i = ii[:derevorder-j]
#            R[i, i+j] = tmp2[i+p_step]
#            r[j-p_step] = tmp2[0]
#        for j in xrange(derevorder, fltLen):
#            tmp[0:j] = 0
#            tmp[j:F] = conjX[j:] * Xb[0:F-j]
#            r[j-p_step] = correlate(iSb, tmp[0:F].conj(), 0)[0]
#            
#        for i in xrange(1, derevorder):
#            j = r_[0:i]
#            R[i, j] = R[j, i].conj();
#        return R, r


#    def setCorr_weave(self, p_step, fltLen, corr, R):
#        code_corr = """
#        for(int i=p_step; i<fltLen; ++i){
#            j = i - p_step;
#            for(int j=0; j<i; ++j){
#                R(i,j) = corr(i-j, p_step + j);
#            }
#        }
#        """
#        from scipy import weave
#        weave.inline(code_corr,['p_step','fltLen','corr','R'], \
#                     type_converters = weave.converters.blitz, \
#                     support_code = "#include <math.h>", libraries = ['m'])
#        return R


#    def calcRr_fast(self, b, _iSb, _Yb):
#        """
#        Calculate modified correlation matrix R using precalculation of X
#        """
#        derevorder, p_step, revorder, fltLen = self.derevorder, self.p_step, self.revorder, self.fltLen
#        
#        idx = r_[0:fltLen]
#        
#        #iSb = r_[_iSb, self.zF]
#        iSb = r_[_iSb, self.zfltLen]
#        fftiSb = fft(iSb)
#        R = empty((derevorder, derevorder), dtype=complex)
##        R = rand(derevorder, derevorder*2).view(complex)   
#
#        corr = ifft(fftiSb * self.XXfftc[b])
#        for i in xrange(p_step, fltLen):
#            j = (i - p_step)
#            R[idx[:-i] + j, idx[:-i]] = corr[j, p_step : fltLen - j]
#            R[idx[:-i], idx[:-i] + j] = corr[j, p_step : fltLen - j].conj()
#
#        Yb = r_[_Yb, self.zfltLen]
#        r = ifft(fft(Yb * iSb) * self.Xfftc[b])[p_step:fltLen]
#        
#        return R, r
#        
#
#    def setXfft(self, X):
#        derevorder, p_step, revorder, fltLen = self.derevorder, self.p_step, self.revorder, self.fltLen
#        (F, B) = X.shape
#        
#        self.XXfftc = [empty((derevorder, F + fltLen), dtype=complex) for b in xrange(B)]
#        self.Xfftc = [None for b in xrange(B)]
#        print (" calculating Xfft [fft length=" + str(F + fltLen) + "]").ljust(40),
#
#        c = Counter(0, B, 10)
#        for b in xrange(B):
#            c.progress(b)
#            Xb = r_[self.zfltLen, X[:, b], self.zfltLen]
#            l = len(Xb)
#            for i in xrange(derevorder):
#                self.XXfftc[b][i, :] = fft(Xb[fltLen:].conj() * Xb[fltLen - i:l - i]).conj()
#            self.Xfftc[b] = fft(Xb[fltLen:]).conj()
#        print " finished[", str(c.time()), "]"
#        return
#    
#    def derev(self, X, G):
#        (F, B) = X.shape
#        Y = empty((F, B), dtype=complex) # dereverberated signal
#        for b in range(B):
#            Y[:, b] = fftfilt(G[:, b], X[:, b])[0:F]
#        return Y
