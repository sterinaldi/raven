import numpy
from setuptools import setup
scripts = ['raven=raven.raven:main']
pymodules = ['raven/raven']

with open("requirements.txt") as requires_file:
    requirements = requires_file.read().split("\n")

with open("README.md") as readme_file:
    long_description = readme_file.read()

setup(
    name = 'raven',
    description = 'RaVeN: Radial Velocities, Nonparametrics',
    author = 'Ste Rinaldi, María Claudia Ramírez-Tannus',
    author_email = 'stefano.rinaldi@uni-heidelberg.de, ramirez@mpia.de',
    url = 'https://github.com/sterinaldi/raven',
    python_requires = '>=3.9',
    packages = ['raven'],
    py_modules = pymodules,
    install_requires=requirements,
    include_dirs = ['figaro', numpy.get_include()],
    entry_points = {
        'console_scripts': scripts,
        },
    version='1.0.0',
    long_description=long_description,
    long_description_content_type='text/markdown',
    )
