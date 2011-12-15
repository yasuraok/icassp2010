# -*- coding: utf-8 -*-
"""
$Date::                           $
$Rev::                            $
$Author::                         $

Source estimator by using harmonic and inharmonic GMM
"""

import numpy as np
from itertools import izip

setup_args={'options': { 'build_ext': {
    'compiler': 'mingw32',
    'include_dirs': [np.get_include()],
}}}
import pyximport; pyximport.install(setup_args=setup_args)
from higmm_core import estim1frame, synth1frame

nfloor = 1e-06 # -120dB noise floor

#import pylab
#pylab.hold(False)
#pylab.ion()

class HIGMMParam(object):
    __slots__ = (
        "numS", # number of harmonic components in each time frame
        "numH", # maximum number of overtones
    )
    
    __defaults__ = {
        "numS": 8,
        "numH": 40,
    }
    
    def __init__(self, **kwargs):
        for k, v in self.__defaults__.items():
            self.__setattr__(k, kwargs.get(k, v))

class HIGMM(object):
    """
    Harmonic-Inharmonic Gaussian Mixture Model
    F: time Frame length
    B: number of frequency Bins
    S: number of Source model in each time frame
    H: number of Harmonics in each source
    """
    
    __slots__ = (
        "f0",  # array[F, S]    Fundamental Frequency
        "amp", # array[F, S, H] amplitude of each harmonic peak
        "var", # array[S]
        "freqs",
        "param",
    )

    def initParam(self, F, B, conf, param = None):
        self.param = param if param != None else HIGMMParam()
        self.freqs = np.linspace(0, conf['fs']/2., B)

        S, H = param.numS, param.numH
        self.f0  = np.maximum(np.random.rand(F, S) * 1200, 50)
        self.amp = np.tile((1.0 / np.arange(1, H + 1) ** 2)[None,None,:], (F, S, 1))
        self.amp[:,:,0] = 0
#        self.amp = np.ones((F, S, H), dtype=float)
        self.var = np.ones((S),       dtype=float) * 450

        # for inharmonic component
        self.f0[:,0] = 600
        self.var[0]    = 30000
                

    def estim(self, X):
        """
        X: spectrogram numpy[F, B]
        """
        AX = np.abs(X)
        F, B = X.shape
        
#        fig = pylab.figure()
#        ax  = fig.add_subplot(111)
        
        print "Estimating HIGMM... ".ljust(40),
        for f in xrange(F):
            f0f, ampf, var = self.f0[f,:], self.amp[f,:,:], self.var
            if int(f / (F/10.)) != int((f+1) / (F/10.)):
                print "*",
            estim1frame(self.freqs, f0f, ampf, var, AX[f,:], self.param)
        print 

    def synth(self):
        F, B = len(self.f0), len(self.freqs)
        AY = np.empty((F, B), dtype=float)
                
        for f in xrange(F):
            f0f, ampf, var = self.f0[f,:], self.amp[f,:,:], self.var
            AY[f,:] = synth1frame(self.freqs, f0f, ampf, var, self.param)

        return AY
