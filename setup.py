"""
BCR and TCR modelling tool
---------

..

"""

import pip
from setuptools import setup, find_packages

setup(
    name='bcr_models',
    version='0.0.0',
    url=None,
    license='Proprietary',
    author='Martin Jespersen, Mads Valdemar Anderson, Michael Schantz Klausen',
    author_email='',
    description='Something something IgG.class',
    long_description=__doc__,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    setup_requires=[

        'numpy>=1.6.0'

        ],
    install_requires=[
        'requests',
        'biopython>=1.63',

        'scikit-learn==0.15.2',
        'SciPy'
    ],
    entry_points = {'console_scripts': [
        'lyra_model = bcr_models.__main__:entry',
        'lyra_pdbdir2tpldb = bcr_models.scripts.pdbdir2tpldb:main',
        'lyra_bench_db = bcr_models.scripts.benchmark_db:main',
        'lyra_bench_dir = bcr_models.scripts.benchmark_dir:main',
    ]},
    package_data={'bcr_models': ['data/','data/pdb','data/rf3/*','data/rf2/*']} ,
)
