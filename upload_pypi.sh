#!/bin/bash

conda activate pypi
python setup.py sdist

twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#pip install --index-url https://test.pypi.org/simple --extra-index-url https://pypi.org/simple metabolabpy==0.6.1

conda deactivate 
#comment
