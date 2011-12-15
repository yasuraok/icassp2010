# -*- coding: utf-8 -*-
"""
$Date::                           $
$Rev::                            $
$Author::                         $

fftを利用した畳み込みの高速計算
scipyにはfftconvolveがあるが，
・2のn乗以外の点数で異常に遅い (おそらく切り分け方法が悪い)
・複素数に対する共役の定義がmatlabと異なる (確か)
という理由により自分で作成
"""

from numpy import ceil, log2, r_, zeros, hstack
from scipy.fftpack import ifft, fft
#from util.whos import whos
MINLEN = 64

def fftfilt_simple(i_ir, i_data, mode='full'):
    """
    区間分割を行わない簡単なFFT畳み込み
    使用メモリが大きい
    """
    type = i_ir.dtype
    minLen = min(len(i_ir), len(i_data))
    maxLen = max(len(i_ir), len(i_data))
    fftlen = 2 ** ceil(log2(maxLen))
    ir = r_[i_ir, zeros(fftlen * 2 - len(i_ir))]
    data = r_[i_data, zeros(fftlen * 2 - len(i_data))]
    if mode == 'same':
        return type.type(ifft(fft(ir) * fft(data))[(minLen - 1) / 2:(minLen - 1) / 2 + maxLen])
    if mode == 'valid':
        return type.type(ifft(fft(ir) * fft(data))[minLen - 1:maxLen])
    else: # mode = 'full'
        return type.type(ifft(fft(ir) * fft(data))[0:minLen + maxLen - 1])

def fftfilt(ii_x, ii_y, mode='full'):
    """
    区間分割を行うFFT畳み込み
    """
    type = ii_x.dtype
    type = ii_y.dtype
    minLen = min(len(ii_x), len(ii_y))
    maxLen = max(len(ii_x), len(ii_y))
    i_x = hstack([ii_x, zeros(MINLEN - len(ii_x), type)]) if len(ii_x) < MINLEN else ii_x
    i_y = hstack([ii_y, zeros(MINLEN - len(ii_y), type)]) if len(ii_y) < MINLEN else ii_y
    fftLen = int(2 ** ceil(log2(min(len(i_x), len(i_y)))))
    if len(i_x) < len(i_y):
        x = hstack([i_x, zeros(fftLen * 2 - len(i_x), type)])
        batchnum = (len(i_y) / fftLen) + int(bool(len(i_y) % fftLen))
        sgnLen = (batchnum + 1) * fftLen
        y = hstack([i_y, zeros(sgnLen - len(i_y), type)])
    else:
        x = hstack([i_y, zeros(fftLen * 2 - len(i_y), type)])
        batchnum = (len(i_x) / fftLen) + int(bool(len(i_x) % fftLen))
        sgnLen = (batchnum + 1) * fftLen
        y = hstack([i_x, zeros(sgnLen - len(i_x), type)])
    z = zeros((sgnLen), dtype=type)
    yy = zeros((fftLen * 2), dtype=type)
    for i in range(batchnum):
        (bgn, batch, end) = (i * fftLen, (i + 1) * fftLen, (i + 2) * fftLen)
        yy[:fftLen] = y[bgn:batch]
        z[bgn:end] += ifft(fft(x) * fft(yy))
    #if mode == 'same':
    #    return type.type()[1:maxLen+1])
    #if mode == 'valid':
    #    return type.type(ifft(fft(ir) * fft(data))[minLen-1:maxLen])
    #else: # mode = 'full'
    if mode == 'same':
        return z[(minLen - 1) / 2:(minLen - 1) / 2 + maxLen]
    if mode == 'valid':
        return z[minLen - 1:maxLen]
    else: # mode = 'full'
        return z[0:minLen + maxLen - 1]

def fftfilt_simple_ar(i_ir, i_ar, i_data, mode='full'):
    type = i_ir.dtype
    minLen = min(len(i_ir), len(i_data))
    maxLen = max(len(i_ir), len(i_data))
    fftlen = 2 ** ceil(log2(maxLen))
    ir = r_[i_ir, zeros(fftlen * 2 - len(i_ir))]
    ar = r_[i_ar[1:], zeros(fftlen * 2 - len(i_ar[1:]))]
    data = r_[i_data, zeros(fftlen * 2 - len(i_data))]
    if mode == 'same':
        return type.type(ifft(fft(ir) / fft(ar) * fft(data))[minLen / 2 - 1:minLen / 2 + maxLen - 1])
    if mode == 'valid':
        return type.type(ifft(fft(ir) / fft(ar) * fft(data))[minLen - 1:maxLen])
    else: # mode = 'full'
        return type.type(ifft(fft(ir) / fft(ar) * fft(data))[0:minLen + maxLen - 1])

class seqfilt:
    default_require = 16
    def __init__(self, filt):
        self.filt = filt
        self.filtlen = len(filt)
        self.buf = zeros((self.filtlen,), dtype=filt.dtype)
        
    def calc(self, time, val):
        pos = time % self.filtlen
        ret = 0
        self.buf[pos] = val
        for p in xrange(self.filtlen):
            d = (pos - p) % self.filtlen
            ret += self.filt[p] * self.buf[d]
        return ret


class seqfiltr:
    default_require = 16
    def __init__(self, filt):
        self.filt = filt
        self.filtlen = len(filt)
        self.mask = (int)(2 ** ceil(log2(len(filt)))) - 1
#        print "mask: ", self.mask
        self.buf = zeros((self.mask + 1,), dtype=filt.dtype)
        
    def calc(self, time, val):
        pos = time & self.mask
#        print pos
        ret = 0
        self.buf[pos] = val
        for p in xrange(self.filtlen):
            d = (pos - p) & self.mask
#            print d
            ret += self.filt[p] * self.buf[d]
        return ret
#
#
#  capacity = Pow2((uint)capacity);
#  this.data = new T[capacity];
#  this.top = this.bottom = 0;
#  this.mask = capacity - 1;


if __name__ == '__main__':
    from scipy import rand, convolve, array
#    from util.testcode import testcode
    
#    xLen = max(int(rand()*10), 1)
#    yLen = max(int(rand()*1000), 1)
#    print xLen, yLen
#    x = rand(xLen)
#    y = rand(yLen)
#    testcode((convolve, fftfilt_simple, fftfilt), x, y, 'same')
#    testcode((convolve, fftfilt_simple, fftfilt), x, y, 'valid')
#    testcode((convolve, fftfilt_simple, fftfilt), x, y)

    x = rand(100)
    y = rand(1000)
    print convolve(x, y)
    for i in xrange(10):
        sf = seqfilt(x)
        [sf.calc(f, v) for f, v in enumerate(y)]
    print "mod ringbuffer"
    for i in xrange(10):
        sf = seqfiltr(x)
        [sf.calc(f, v) for f, v in enumerate(y)]
    print "and ringbuffer"

