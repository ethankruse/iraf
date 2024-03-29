#!/bin/bash -x

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh;
else
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh;
fi
# accept all defaults without human input
bash ~/miniconda.sh -b -p "$HOME"/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

# no human input needed anymore
conda config --set always_yes yes
# quiet update mode
conda update -q conda
# list the current conda environment info
conda info -a
conda create --yes -n test python="$PYTHON_VERSION"
source activate test
conda install -q numpy scipy matplotlib h5py setuptools pytest pytest-cov pip astropy
# pip install coveralls

# Build the extension
python setup.py develop
