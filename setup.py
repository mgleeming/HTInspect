#!/bin/env python

from setuptools import setup, find_packages

setup(
    name = 'viewHTResults',
    version = '1.0.0',
    packages = find_packages(),
    include_package_data = True,
    entry_points = {
        'gui_scripts': ['viewHTResults = viewHTResults.viewHTResults:main']
    },
)
