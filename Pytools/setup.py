from setuptools import setup
from setuptools import find_packages

setup(name="west",
      version="4.1.0",
      packages=find_packages(),
      description='installation script for WEST',
      url='https://west-code.org',
      author='Marco Govoni',
      author_email='mgovoni@anl.gov',
      license='GPLv3',
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'pyyaml',
          'datetime',
          'setuptools'
      ],
      python_requires='>=3.6, <4',
      zip_safe=True)
