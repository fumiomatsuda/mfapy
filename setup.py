# -*- coding: utf-8 -*-
# Learn more: https://github.com/fumiomatsuda/mfapy
#

from setuptools import setup, find_packages
import sys
import os, os.path
pyversion = sys.version_info[0]
pkgdir = 'python%s' % pyversion

def load_requires_from_file(fname):
    if not os.path.exists(fname):
        raise IOError(fname)
    return [pkg.strip() for pkg in open(fname, 'r')]

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='mfapy',
    version='0.6.3',
    package_dir = {'mfapy': ""},
    py_modules =  ["metabolicmodel", "carbonsource", "mfapyio", "mdv", "optimize"],
    description = "mfapy: a python package for the 13C metabolic flux analysis",
    platforms = ["Windows", "Linux", "Unix"],
    long_description=readme,
    author='Fumio Matsuda',
    author_email='fmatsuda@osaka-u.ac.jp',
    install_requires=load_requires_from_file('requirements.txt'),
    url='https://github.com/fumiomatsuda/mfapy',
    license=license,

    test_suite='tests'
)
