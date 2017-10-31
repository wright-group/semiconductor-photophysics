# --- import -------------------------------------------------------------------------------------


import os

from setuptools import setup, find_packages


# --- define -------------------------------------------------------------------------------------


here = os.path.abspath(os.path.dirname(__file__))


# --- setup --------------------------------------------------------------------------------------


extra_files = []
extra_files.append(os.path.join(here, 'LICENSE'))
extra_files.append(os.path.join(here, 'VERSION'))


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
    license='MIT'
)
