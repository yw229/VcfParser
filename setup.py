from setuptools import setup
import sys

requires = []

# install ordereddict for python < 2.7
py_version = sys.version_info
if py_version.major == 2 and py_version.minor < 7:
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
    version='0.2.1'
)
