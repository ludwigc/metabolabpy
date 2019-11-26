#!/bin/bash

conda create -n pypi
conda activate pypi

PACKAGE="metabolabpy"
VERSION="0.0.1"
echo $PACKAGE

git clone https://github.com/ludwigc/metabolabpy
cd $PACKAGE
python setup.py test
python setup.py sdist

conda install -c conda-forge twine

twine upload --repository-url https://test.pypi.org/legacy dist/*
pip install --index-url https://test.pypi.org/simple --extra-url https://pypi.org/simple $PACKAGE

twine upload dist/*
pip install -U $PACKAGE

conda deactivate pypi
conda remove --name pypi --all
