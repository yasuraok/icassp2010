Example program of ICASSP2010 proceeding: MUSIC DEREVERBERATION USING HARMONIC STRUCTURE SOURCE MODEL AND WIENER FILTER.

# setup
The author is using Mac OS 10.11, Python2.7 on Anaconda. Cython is also used for fast numerical computation.

1. install homebrew (see: http://brew.sh/)
2. install pyenv

   ```
   $ brew install pyenv
   $ echo 'export PYENV_ROOT="${HOME}/.pyenv"' >> ~/.bash_profile
   $ echo 'export PATH="${PYENV_ROOT}/bin:$PATH"' >> ~/.bash_profile
   $ echo 'eval "$(pyenv init -)"' >> ~/.bash_profile
   ```

3. install Anaconda 2.7

   ```
   pyenv install anaconda2-4.1.1
   pyenv global anaconda2-4.1.1
   ```

4. install CMake (see: https://cmake.org/download/)

# usage

1. build Cython C++ extension

   ```
   cmake ./
   ```

   (If you see errors, please check how CMakeLists.txt works)

2. prepare sample \*.wav file and execute entry python script

   ```
   python prg/icassp2010.py sample.wav
   ```
