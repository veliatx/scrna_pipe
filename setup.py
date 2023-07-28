# -*- coding: utf-8 -*-

from os.path import abspath, dirname
from sys import path
from setuptools import setup, find_packages

setup(
    name='scrna_pipe',
    version='0.0.1',
    description="""Velia scRNA pipelines""",
    author='Stephen Federowicz',
    author_email='steve@veliatx.com',
    classifiers=[
        'Programming Language :: Python :: 3.7',
    ],
    keywords='microproteins',
    packages=find_packages(),
    install_requires=[
        'SQLAlchemy>=1.3.10,<2.0',
        'numpy>=1.17.2',
        'anndata2ri',
        'scanpy',
        'anndata',
        'sc_toolbox',
        'scipy>=1.3.1',
        'pytest>=4.6.6',
        'six>=1.12.0',
        'tornado>=4.5.3',
        'configparser>=4.0.2',
        'click',
    ],
    entry_points = {
        'console_scripts': ['scrna_preprocess=scrna_pipe.preprocess:run_scale',
                            'scrna_qc=scrna_pipe.preprocess:run_qc',
                            'scrna_annotation=scrna_pipe.annotation:main',
                            'scrna_differential=scrna_pipe.differential:main'],
    }
)