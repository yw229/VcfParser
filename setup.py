from setuptools import setup
import sys

requires = []

# python 2.6 does not have argparse
try:
    import argparse
except ImportError:
    requires.append('argparse')


setup(
    name='PyVCF',
    py_modules=['vcf', 'vcf_filter'],
    scripts=['vcf_melt', 'vcf_filter.py'],
    author='James Casbon and @jdoughertyii',
    author_email='casbon@gmail.com',
    description='Variant Call Format (VCF) parser for python',
    test_suite='test.test_vcf.suite',
    requires=requires,
    entry_points = {
        'vcf.filters': [
            'site_quality = vcf_filter:SiteQuality',
            'vgq = vcf_filter:VariantGenotypeQuality',
        ]
    },
    url='https://github.com/jamescasbon/PyVCF',
    version='0.2.2'
)
