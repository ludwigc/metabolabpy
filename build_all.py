#!python

import metabolabpy.__init__ as ml_version
import os

os.system('python setup.py sdist')
system_cmd = 'twine upload dist/metabolabpy-' + ml_version.__version__ + '.tar.gz'
os.system(system_cmd)

os.chdir('../qtmetabolabpy')
os.system('python setup.py sdist')
system_cmd = 'twine upload dist/qtmetabolabpy-' + ml_version.__version__ + '.tar.gz'
os.system(system_cmd)

os.chdir('../metabolabpy')
