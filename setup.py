from setuptools import setup, find_packages
from MetaCanSNPer.CommandLineParser import __version__

setup(
    name="MetaCanSNPer",
    version=__version__,
    url=None, #"https://git-int.foi.se/bioinfo/metacansnper",
    description="MetaCanSNPer: A toolkit for SNP-typing genomes.",

    # Author details
    author="Fredrik Sörensen",
    author_email="fredrik.sorensen@foi.se",

    license='GNU GENERAL PUBLIC LICENSE version 3',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    python_requires=">=3.6",

    keywords="Bioinformatics SNP-typing sequence-data",

    install_requires=["PseudoPathy"],
    packages=find_packages(exclued=['contrib', 'docs', 'test*']),
    entry_points={"console_scripts": [
                    "MetaCanSNPer=MetaCanSNPer.CommandLineParser:main"
    ]})
