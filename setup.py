#!/usr/bin/env python3
from setuptools import setup

setup(name = 'BRB',
      version = '0.0.1',
      description = 'dixitque PI fiat computationēs et facta est computationēs et vidit PI computationēs quod esset bona',
      author = "Devon P. Ryan",
      author_email = "ryan@ie-freiburg.mpg.de",
      scripts = ['bin/BigRedButton.py'],
      packages = ['BRB'],
      include_package_data = False,
      install_requires = ['configparser',
                          'bioblend'])
