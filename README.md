Example program of ICASSP2010 proceedings: MUSIC DEREVERBERATION USING HARMONIC STRUCTURE SOURCE MODEL AND WIENER FILTER.

# setup
The author is using Mac OS 10.11, Python2.7 on Anaconda. Cython is also used for fast numerical computation.

1. Install homebrew (see: http://brew.sh/)
2. Install pyenv

   ```
   $ brew install pyenv
   $ echo 'export PYENV_ROOT="${HOME}/.pyenv"' >> ~/.bash_profile
   $ echo 'export PATH="${PYENV_ROOT}/bin:$PATH"' >> ~/.bash_profile
   $ echo 'eval "$(pyenv init -)"' >> ~/.bash_profile
   ```

3. Install Anaconda 2.7

   ```
   pyenv install anaconda2-4.1.1
   pyenv global anaconda2-4.1.1
   ```

4. Install CMake (see: https://cmake.org/download/)

# usage

1. Build Cython C++ extension.

   ```
   cmake ./
   make
   ```

   (If you see errors, please check how CMakeLists.txt works)

2. Prepare sample \*.wav file and execute entry python script.

   ```
   python prg/icassp2010.py sample.wav
   ```
   
   Two output files will be generated. Both are given by inverse STFT of dereverberated spectrogram.
   
   - sample_derev: dereverberated signal by using inverse linear filtering, which is mentioned at the last of section 2.2 in the proceedings
   - sample_derev_wiener: dereverberated signal by using the proposed Wiener filtering, which is mentioned at section 3.1 in the proceedings
