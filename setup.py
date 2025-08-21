# -*- coding: utf-8 -*-
# SPDX-License-Identifier: LGPLv3
"""
A minimal setup script for Anthropogenic mw.

All the remaining configuration is in pyproject.toml.
"""
from setuptools import setup, find_packages


def current_version():
    changelog = open('CHANGELOG.md', 'r')
    changelog_items = changelog.readlines()
    version = "0.0.0"
    version_date = "2025-01-01"
    for item in changelog_items:
        if item.startswith("##"):
            pos = item.split();
            version = pos[1][1:]
            version_date = pos[3]
    return version, version_date


Version, VersionDate = current_version()


setup(
    name='ha3py',
    version=Version,
    date=VersionDate,
    description='HyPy - Python package for assessing earthquake recurrence parameters from Incomplete Data Files.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Jan Wiszniowski, Andrzej Kijko',
    author_email='jwisz@igf.edu.pl, andrzej.kijko@up.ac.za',
    url='https://github.com/JanWiszniowski/ha3py',
    packages=find_packages(),
    install_requires=[
        'matplotlib>=3.9.2',
        'numpy>=1.2.0',
        'scipy>=1.0.0',
        # 'obspy>=1.2.0',
    ],
    license='LGPLv3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU Lesser General Public License, Version 3',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ],
    entry_points={
        'console_scripts': ['ha3 = ha3py.ha3:main',
                            'ha_m_max = ha3py.compute_m_max:main',
                            'ha_config = ha3py.configuration:main',
                            'ha_import = ha3py.import_to_hapy:main',
                            'ha_compute = ha3py.compute:main',
                            'ha_print = ha3py.print_results:main',
                            'ha_plot = ha3py.plot_results:main',
                            'ha_simulate = ha3py.synthetic_data_generation:main',
                            ],
    },
    keywords='seismology, seismic hazard, maximum magnitude, incomplete data files',
)
