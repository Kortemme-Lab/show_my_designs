#!/usr/bin/env python2

from setuptools import setup

setup(
    name='show_my_designs',
    version='1.2.2',
    author='Kale Kundert',
    author_email='kale.kundert@ucsf.edu',
    url='https://github.com/Kortemme-Lab/show_my_designs',
    license='GPLv3',
    description=" A GUI for visualizing and interacting with score vs RMSD plots.",
    long_description=open('README.rst').read(),
    packages=[
        'show_my_designs',
    ],
    package_data={
        'show_my_designs': ['*.png'],
    },
    install_requires=[
        #'pygtk' is required, but pip usually can't install it.
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
