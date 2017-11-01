#! /usr/bin/env python


# --- import -------------------------------------------------------------------------------------


import os

from setuptools import setup, find_packages


# --- define -------------------------------------------------------------------------------------


here = os.path.abspath(os.path.dirname(__file__))


extra_files = []
extra_files.append(os.path.join(here, 'CONTRIBUTORS'))
extra_files.append(os.path.join(here, 'LICENSE'))
extra_files.append(os.path.join(here, 'README.md'))
extra_files.append(os.path.join(here, 'VERSION'))


# --- setup --------------------------------------------------------------------------------------


with open(os.path.join(here, 'requirements.txt')) as f:
    required = f.read().splitlines()


with open(os.path.join(here, 'VERSION')) as version_file:
    version = version_file.read().strip()


setup(
    name='semiconductor_photophysics',
    version=version,
    packages=find_packages(),
    package_data={'': extra_files},
    install_requires=required,
    author='Darien Morrow',
    author_email='darienmorrow@gmail.com',
    license='MIT',
    url='https://github.com/wright-group/semiconductor_photophysics',
    keywords='semiconductor photophysics spectroscopy science',
    classifiers=['Development Status :: 0 - Planning',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering',
                 'Natural Language :: English',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5']
)
