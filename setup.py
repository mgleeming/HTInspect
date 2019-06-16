#!/bin/env python

from setuptools import setup, find_packages

setup(
    name = 'HTInspect',
    version = '1.0.0',
    packages = find_packages(),
    include_package_data = True,
    entry_points = {
        'gui_scripts': ['HTInspect = HTInspect.HTInspect:main']
    },
)
