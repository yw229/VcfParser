from setuptools import setup

requires = []

# python 2.6 does not have argparse
try:
    import argparse
except ImportError:
    requires.append('argparse')

# get the version without an import
VERSION = "Undefined"
DOC = ""
inside_doc = False
for line in open('vcf/__init__.py'):
    if "'''" in line:
        inside_doc = not inside_doc
    if inside_doc:
        DOC += line.replace("'''", "")

    if (line.startswith('VERSION')):
        exec(line.strip())


setup(
    name='PyVCF',
    packages=['vcf', 'vcf.test'],
    scripts=['scripts/vcf_melt', 'scripts/vcf_filter.py'],
    author='James Casbon and @jdoughertyii',
    author_email='casbon@gmail.com',
    description='Variant Call Format (VCF) parser for Python',
    long_description=DOC,
    test_suite='vcf.test.test_vcf.suite',
    install_requires=['distribute'],
    requires=requires,
    entry_points = {
        'vcf.filters': [
            'site_quality = vcf.filters:SiteQuality',
            'vgq = vcf.filters:VariantGenotypeQuality',
        ]
    },
    url='https://github.com/jamescasbon/PyVCF',
    version=VERSION,
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
      ],
    keywords='bioinformatics',
    use_2to3=True,
    include_package_data=True,
    package_data = {
        '': ['*.vcf', '*.gz', '*.tbi'],
        }
)
