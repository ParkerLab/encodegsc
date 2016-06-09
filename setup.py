#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


requirements = [
    'matplotlib',
    'pylab',
    'rpy2'
]

setup(
    name='encode_gsc',
    version='0.9.0',
    description='ENCODE Genome Structure Correction Tools',
    long_description="""Assessing the significance of observations within large scale genomic studies using random subsampled genomic region is a difficult problem because there often exists a complex dependency structure between observations. GSC is a data subsampling approach based on a block stationary model for genomic features to alleviate the hidden dependencies. This model is motivated by earlier studies of DNA sequences, which show that there are global shifts in base composition, but that certain sequence characteristics are locally unchanging.""",
    author="ENCODE",
    url='https://www.encodeproject.org/software/gsc/',
    packages=[
        'encode_gsc',
    ],
    scripts=[
        'block_bootstrap.py',
        'segmentation.py',
    ],
    package_dir={'encode_gsc': 'encode_gsc'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='ENCODE Genome Structure Correction',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
