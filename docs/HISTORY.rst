Changes
=======

Pending
-------

* Fix setup.py for python < 2.7
* Add ``__eq__`` to ``_Record``
* Add ``is_het`` and ``is_variant`` to ``_Call``

0.2.2 Release
-------------

Documentation release

0.2.1 Release
-------------

* Add shebang to vcf_filter.py

0.2 Release 
-----------

* Replace genotype dictionary with a ``Call`` object
* Methods on ``Record`` and ``Call`` (thanks @arq5x)
* Shortcut parse_sample when genotype is None

0.1 Release 
-----------

* Added test code
* Added Writer class
* Allow negative number in ``INFO`` and ``FORMAT`` fields (thanks @martijnvermaat)
* Prefer ``vcf.Reader`` to ``vcf.VCFReader``
* Support compressed files with guessing where filename is available on fsock
* Allow opening by filename as well as filesocket
* Support fetching rows for tabixed indexed files
* Performance improvements (see ``test/prof.py``)
* Added extensible filter script (see FILTERS.md), vcf_filter.py 

Contributions
-------------

Project started by @jdoughertyii and taken over by @jamescasbon on 12th January 2011.


