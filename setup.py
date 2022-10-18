#!/usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools
import sys
import metabolabpy
#from setuptools.command.install import install
#import os
#import shutil
#
#class CustomInstall(install):
#    def run(self):
#        install.run(self)
#        if sys.platform == 'darwin':
#          homeDir = os.path.expanduser("~")
#          fName  = homeDir + os.path.sep + 'appify.metabolabpy'
#          fName2 = homeDir + os.path.sep + 'createStarter.metabolabpy'
#          f = open(fName, 'w')
#          f.write('#!/usr/bin/env bash\n')
#          f.write('\n')
##          f.write('APPNAME=${2:-$(basename "$1" ".sh")}\n')
#          f.write('DIR="$APPNAME.app/Contents/MacOS"\n')
#          f.write('\n')
#          f.write('if [ -a "$APPNAME.app" ]; then\n')
#          f.write('  echo "$PWD/$APPNAME.app already exists :("\n')
#          f.write('  exit 1\n')
#          f.write('fi\n')
#          f.write('\n')
#          f.write('mkdir -p "$DIR"\n')
#          f.write('cp "$1" "$DIR/$APPNAME"\n')
#          f.write('chmod +x "$DIR/$APPNAME"\n')
#          f.write('\n')
#          f.write('echo "$PWD/$APPNAME.app"\n')
#          f.close()
#          os.chmod(fName, 0o777)
#          import metabolabpy
#          mlStarterPath = os.path.dirname(metabolabpy.__file__) + os.path.sep + 'mlStarter'
#          f = open(fName2, 'w')
#          f.write('#!/usr/bin/env bash\n')
#          f.write('cd ' + homeDir + '\n')
#          f.write(fName + ' $(which metabolabpy) MetaboLabPy\n')
#          f.write('cp -r ' + mlStarterPath + os.path.sep + 'Contents/* ' + homeDir + os.path.sep + 'MetaboLabPy.app' + os.path.sep + 'Contents\n')
#          f.write('cp ' + mlStarterPath + os.path.sep + 'Icon ' + homeDir + os.path.sep + 'MetaboLabPy.app' + os.path.sep + "Icon$'\\r'\n")
#          f.write('\n')
#          f.write('mv ' + homeDir + os.path.sep + 'MetaboLabPy.app /Applications\n')
#          f.close()
#          os.chmod(fName2, 0o777)
#          if os.path.isdir('/Applications/MetaboLabPy.app'):
#            shutil.rmtree('/Applications/MetaboLabPy.app')
#            
#          os.system(homeDir + os.path.sep + 'createStarter.metabolabpy')
#          os.remove(homeDir + os.path.sep + 'appify.metabolabpy')
#          os.remove(homeDir + os.path.sep + 'createStarter.metabolabpy')
#          
#


def main():

    if not ((sys.version_info[0] != 3) or (sys.version_info[1] >= 6)):
        sys.exit("Python >=3.6 is required ")

    # read the contents of your README file
    with open('README.rst', encoding='utf-8') as f:
        long_description = f.read()

    setuptools.setup(name="metabolabpy",
        version=metabolabpy.__version__,
        description="Python package for data processing of NMR 1D and 2D metabolomics and metabolism tracing data",
        long_description=long_description,
        author="Christian Ludwig",
        author_email="C.Ludwig@bham.ac.uk ",
        url="https://github.com/ludwigc/metabolabpy",
        license="GPLv3",
        platforms=['MacOS, Windows, UNIX'],
        keywords=['Metabolomics', 'Tracing', 'NMR', 'Data Processing', '13C'],
        packages=setuptools.find_packages(),
        test_suite='tests.suite',
        install_requires=open('requirements.txt').read().splitlines(),
        include_package_data=True,
        classifiers=[
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.9",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
        ],
        #entry_points={
        # 'console_scripts': [
        #     'metabolabpy = metabolabpy.__main__:main'
        # ]
        #}
    )


if __name__ == "__main__":
    main()
  
