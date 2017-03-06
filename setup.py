#!/usr/bin/env python2

from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES

# Install data files into the same directory as source files.
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

with open('README.rst') as file:
    readme = file.read()

setup(
    name='show_my_designs',
    version='1.2.0',
    author='Kale Kundert',
    author_email='kale.kundert@ucsf.edu',
    url='https://github.com/Kortemme-Lab/show_my_designs',
    license='GPLv3',
    description=" A GUI for visualizing and interacting with score vs RMSD plots.",
    long_description=readme,
    py_modules=[
        'show_my_designs',
    ],
    data_files=[
        'show_my_designs.png',
    ],
    install_requires=[
        #'pygtk',
        'docopt',
        'pyyaml',
        'matplotlib',
        'numpy',
        'scipy',
        'pandas',
        'numexpr',
    ],
    entry_points = {
        'console_scripts': ['show_my_designs=show_my_designs:main'],
    },
)
