# installation for velocyto on macOS 10.14.6

conda create -n velocyto python=3.7 numpy scipy cython numba matplotlib scikit-learn h5py click
conda activate velocyto
brew install gcc
# figure out path to homebrew gcc

CC=/usr/local/Cellar/gcc/11.2.0_3/bin/gcc-11 pip install velocyto


