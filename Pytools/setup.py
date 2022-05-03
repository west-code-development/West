from setuptools import setup
from setuptools import find_packages
import json

with open('../VERSION.json',"r") as file:
    data = json.load(file)

setup(name=data["name"],
      version=data["version"],
      packages=find_packages(),
      description='installation script for WEST',
      url=data["url"],
      author='Marco Govoni',
      author_email='mgovoni@anl.gov',
      license=data["license"],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'pyyaml',
          'datetime',
          'setuptools',
          'sphinx_rtd_theme'
      ],
      python_requires='>=3.6, <4',
      zip_safe=True)

