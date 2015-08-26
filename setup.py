#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = ['matplotlib==1.3.1', 'pysam==0.8.3']

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='riboplot',
    version='0.1.0',
    description="Plot read counts of RiboSeq data from BAM format alignment files",
    long_description=readme + '\n\n' + history,
    author="Vimalkumar Velayudhan",
    author_email='vimalkumarvelayudhan@gmail.com',
    url='https://github.com/vimalkumarvelayudhan/riboplot',
    packages=[
        'riboplot',
    ],
    package_dir={'riboplot':
                 'riboplot'},
    install_requires=requirements,
    license="GPL",
    zip_safe=False,
    keywords='riboplot',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    entry_points={
        'console_scripts': [
            'riboplot = riboplot.riboplot:run',
            'ribocount = riboplot.ribocount:run']
    },
    package_data={'riboplot': ['data/*.html', 'data/css/*.gif', 'data/css/*.css',
                               'data/js/*.js', 'data/js/*.map', 'scripts/*.sh']},
    exclude_package_data={'tests': ['data/*.bam', 'data/*.bai', 'data/*.fna', 'data/*.fai']}
)
