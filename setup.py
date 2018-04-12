#!/usr/bin/env python3
from setuptools import setup

setup(name = 'BRB',
      version = '0.0.1',
      description = 'Push button, get results, ???, publish',
      author = "Devon P. Ryan",
      author_email = "ryan@ie-freiburg.mpg.de",
      scripts = ['bin/BigRedButton.py'],
      packages = ['BRB'],
      include_package_data = False,
      install_requires = ['configparser',
                          'bioblend'])
