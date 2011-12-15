# -*- coding: utf-8 -*-
"""
$Date::                           $
$Rev::                            $
$Author::                         $

wavファイルとflacファイルを読み書きする
"""

import sys, os, wave, subprocess, os.path as osp
from cStringIO import StringIO
from scipy import array, float64, int16, maximum, minimum, empty, int8, int32, frombuffer

#FLACPATH = "flac" # please specify if you need to change
FLACPATH = "C:/Program Files (x86)/FLAC/flac.exe"

typeinfo = {
  1 : int8,
  2 : int16,
  4 : int32,
}

class Error(Exception):
    pass

def audioreadmono(wavfile, normalize=False):
    print "Reading wave file [", wavfile, "] ...",
    (x, fs, bits) = audioread(wavfile, normalize=normalize);
    x = (x.sum(1) / x.shape[1])
    print "done."
    return (x, fs, bits)

def audioread(filename, normalize=False, dtype=float64):
    ext = osp.splitext(filename)[1]
    if   ext.upper() == '.WAV':
        return wavread(filename, normalize, dtype)
    elif ext.upper() == '.FLAC':
        return flacread(filename, normalize, dtype)
    else:
        raise Error, 'not a supported audio file format: ' + filename

def audiowrite(x, fs, bits, filename, normalize=False, compress=True):
    ext = osp.splitext(filename)[1]
    if   ext.upper() == '.WAV':
        return wavwrite(x, fs, bits, filename, normalize, compress)
    elif ext.upper() == '.FLAC':
        return flacwrite(x, fs, bits, filename, normalize, compress)
    else:
        raise Error, 'not a supported audio file format: ' + filename

def wavread(filename, normalize=False, dtype=float64):
    win = wave.open(filename, 'r')
    length = win.getnframes() * win.getsampwidth() * win.getnchannels()
    nframes = win.getnframes()
    x = frombuffer(win.readframes(length), dtype='c').view(typeinfo[win.getsampwidth()]).reshape(nframes, win.getnchannels()).astype(dtype)
    win.close()
    if normalize:
        x /= abs(x).max()
    else:
        x = x / (2 ** (win.getsampwidth() * 8 - 1))
#    print "==============================================12341234==="
    return x, win.getframerate(), win.getsampwidth() * 8
       

def wavwrite(i_x, fs, bits, filename, normalize=False, compress=True):
    if i_x.ndim == 1:
        x = i_x.reshape(len(i_x), 1)
    elif i_x.ndim == 2:
        x = i_x
    else:
        assert False
    
    win = wave.open(filename, 'w')
    win.setnchannels(int16(x.shape[1]))
    win.setsampwidth(int16(bits / 8))
    win.setframerate(fs)
    if compress:
        positivelimit = ((2 ** (bits - 1) - 1) / float(2 ** (bits - 1)))
        x = maximum(-1.0, minimum(x, positivelimit))
#    if compress:
#        x = maximum(-1.0, minimum(x, 1.0))
    if normalize:
        x /= abs(x).max()
    x = x * (2 ** (bits - 1));
    win.writeframes(array(x, int16).tostring())
    win.close()
    return

class wavwritestream():
    def __init__(self, channel, fs, bits, filename, normalize=False, compress=True):
        self.filename = filename
        self.win = wave.open(filename, 'w')
        self.channel = channel
        self.win.setnchannels(int16(channel))
        self.bits = bits
        self.win.setsampwidth(int16(bits / 8))
        self.fs = fs
        self.win.setframerate(fs)
        self.normalize = normalize
        self.compress = compress
        return
    def close(self):
        self.win.close()
        return
    def write(self, i_x):
        if i_x.ndim == 1:
            x = i_x.reshape(len(i_x), 1)
        elif i_x.ndim == 2:
            x = i_x
        else:
            assert False
        if self.compress:
            positivelimit = ((2 ** (self.bits - 1) - 1) / float(2 ** (self.bits - 1)))
            x = maximum(-1.0, minimum(x, positivelimit))
        if self.normalize:
            x /= max(abs(x))
        x = x * (2 ** (self.bits - 1))
        self.win.writeframes(array(x, int16).tostring())
        return
    def path(self):
        return self.filename
    def __del(self):
        try:
            self.close()
        except Exception:
            print Exception 
        return

class flacwritestream(wavwritestream):
    def __init__(self, channel, fs, bits, filename, normalize=False, compress=True):
        self.wavfile = chext(filename, "wav")
        self.flacfile = filename
        open(self.flacfile, 'w').close() # lock file
        self.win = wave.open(self.wavfile, 'w')
        self.channel = channel
        self.win.setnchannels(int16(channel))
        self.bits = bits
        self.win.setsampwidth(int16(bits / 8))
        self.fs = fs
        self.win.setframerate(fs)
        self.normalize = normalize
        self.compress = compress
        return
    def close(self):
        self.win.close()
        command = FLACPATH + " --delete-input-file " + self.wavfile + " -f -o " + self.flacfile
        #print command
        ret = perform(command)
        return
    def path(self):
        return self.flacfile
    def __del__(self):
        try:
            self.close()
        except Exception:
            print Exception
        return
        
def checkexists(path):
    if not osp.exists(path):
        raise Error, 'File does not exist: ' + path


def flacread(flacfile, normalize=False, dtype=float64):
    """
    まずwavで吐く→flacをコマンド呼び出しして処理
    """
    wavfile = chext(flacfile, "wav")
    command = [FLACPATH, "-d", flacfile, "-f", "-o", wavfile]
    perform(command)
    (x, fs, bits) = wavread(wavfile, normalize, dtype)
    try: os.remove(wavfile)
    except OSError: pass# [Errno 17] File exists:
    return (x, fs, bits)

def flacwrite(x, fs, bits, flacfile, normalize=False, compress=True):
    """
    まずwavで吐く→flacをコマンド呼び出しして処理
    """
    open(flacfile, 'w').close() # lock file
    wavfile = chext(flacfile, "wav")
    wavwrite(x, fs, bits, wavfile, normalize, compress)
    command = [FLACPATH, "--delete-input-file", wavfile, "-f", "-o", flacfile]
    #print command
    ret = perform(command)
    print ret
    return

def flacreadUsingPipe(flacfile, normalize=False, dtype=float64):
    """
    flacのデコード結果をstdoutに吐いてパイプから得る方法
    unixでしかうまくいかなかったので注意
    """
    checkexists(flacfile)
#    command = "flac -d "+flacfile+" --silent --stdout"
#    wav = (os.popen(command).read() if sys.platform == 'win32' else commands.getoutput(command))
    p = subprocess.Popen([FLACPATH , "-d", flacfile, "--silent", "--stdout"],
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
#                        stderr = subprocess.STDOUT,
#                        close_fds=True)
    #wav = p.wait()
    (out, stdin) = (p.stdout, p.stdin)
#    return wavio.wavread(StringIO(wav))
    return wavread(out, normalize=normalize, dtype=dtype)

def flacwriteUsingPipe(i_x, fs, bits, flacfile, normalize=False, compress=True):
    """
    flacのエンコード結果をstdinに吐いてパイプから得る方法
    unixでしかうまくいかなかったので注意
    """
    wav = StringIO()
    wavwrite(i_x, fs, bits, wav, normalize=normalize, compress=compress)
    
    p = subprocess.Popen([FLACPATH, "-o", flacfile, "-f", "-"],
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
#                        stderr = subprocess.PIPE,
#                        close_fds=True)
#    p.stdin.write(wav,get)
    (out, err) = p.communicate(input=wav.getvalue())
#    p.terminate()
    return out

                                      
def chext(str, ext):
    return osp.splitext(str)[0] + "." + ext

def flacencode(wavfile, delete=False):
    flacfile = chext(wavfile, "flac")
    if delete == True:
        command = [FLACPATH, "--delete-input-file", wavfile, "-f", "-o", flacfile]
    else:
        command = [FLACPATH, wavfile, "-f", "-o", flacfile]
    perform(command)
    return flacfile

def flacdecode(flacfile, delete=False):
    wavfile = chext(flacfile, "wav")
    command = [FLACPATH, "-d", flacfile, "-f", "-o", wavfile]
    perform(command)
    if delete:
        os.remove(flacfile)
    return wavfile

def perform(command, debug = True):
    """
    command: コマンドをスペースで区切るところで切り分けたリスト 例： ["C:/hoge foo/bar.exe", "-v"]
    debug:   Trueならコマンドとコマンド実行結果の標準出力をprint
    この関数を使うメリットは
    1) 標準出力が取り出せる
    2) 1を実現するPython標準関数commands.getoutputがwin32では動作しないを回避
    3) ファイルパスを渡す時のスペース回避を自動処理 (Subprocessモジュールの機能)
    4) unicode文字列を自動処理
    5) windowsなら不要なDOS窓がでない？
    """
    if debug: print " ".join(command)
    isWin = sys.platform == "win32"
    encoding = "cp932" if isWin else "UTF-8"
    command = [unicode(c).encode(encoding) for c in command]
    p = subprocess.Popen(command, shell=isWin, stdin=subprocess.PIPE , stdout=subprocess.PIPE)
    ret = p.stdout.read()
    if debug: print ret
    return ret

if __name__ == '__main__':
    args = [
    "C:/Users/yasuraok/Slides/ppt_sound/orig/RM-J022_dry.wav",
    "C:/Users/yasuraok/Slides/ppt_sound/orig/RM-J022_wet.wav",
    "C:/Users/yasuraok/Slides/ppt_sound/orig/RM-J022_wet_derev_i.wav",
    "C:/Users/yasuraok/Slides/ppt_sound/orig/RM-J022_wet_derev_is.wav",
    ]
    thresh = 0.05
    
    def argthresh(seq, thresh):
        for i, s in enumerate(seq):
            if s > thresh:
                return i
        return len(seq)

    for wav in args:
        (x, fs, bits) = wavread(wav)
        x = x / x.max() # normalize
        trim = argthresh(x.sum(1), thresh)
        wavwrite(x [trim:, :], fs, bits, wav)
        print wav, " is normalized and trimmed."

    
"""
def allocate(shape, fs, bits, filename):
    win = wave.open(filename, 'w')
    win.setnchannels(int16(shape[1]))
    win.setsampwidth(int16(bits / 8))
    win.setframerate(fs)
    buf = empty(shape, int16).tostring()
    win.writeframes(buf)
    win.close()

def rewrite(filename)
"""

    #print "done"
    #x = array(y, copy=False).view(int16).reshape(win.getnframes(), win.getnchannels())
    #x = array(list(buf)).view(int16).reshape(win.getnframes(), win.getnchannels())
    #print "#####===============================================",memory(),"=============="
    #struct.unpack('%dh' % (len(buf) / win.getsampwidth()), buf)
    #x = array(buf, dtype)
    #print "#####===============================================",memory(),"=============="
    #del buf, y
    #print "#####===============================================",memory(),"=============="
    #x = reshape(x, [win.getnframes(), win.getnchannels()])
    #whos(locals())




# #   raw_input()
# #   print "call"
#    win    = wave.open(filename, 'r')
# #   print "read"
# #   print win.getnframes(), win.getsampwidth(), win.getnchannels()
##    raw_input()
#    length = win.getnframes() * win.getsampwidth() * win.getnchannels()
#    buf = win.readframes(length)
# #   print len(buf)
#    if len(buf) < length:
#        nframes = len(buf) / (win.getsampwidth() * win.getnchannels())
#    else:
#        nframes = win.getnframes()
##    raw_input()
#    win.close()
#    x = empty((nframes * win.getsampwidth() * win.getnchannels(), ), dtype='c')
#    x[:] = buf[:nframes * win.getsampwidth() * win.getnchannels()]
#    x = array(x.view(typeinfo[win.getsampwidth()]).reshape(nframes, win.getnchannels()), dtype)
#    if normalize:
#        x /= abs(x).max()
#    else:
#        x = x / (2 ** (win.getsampwidth() * 8 - 1))
#    return x, win.getframerate(), win.getsampwidth() * 8
