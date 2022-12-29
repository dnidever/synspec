.. synspec documentation master file, created by
   sphinx-quickstart on Tue Feb 16 13:03:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*******
Synspec
*******

Introduction
============
|synspec| is a generical stellar spectral synthesis package that can be run from python.  It's essentially
a redistribution of the `Synspec <http://tlusty.oca.eu/Synspec49/synspec.html>`_
synthesis software by Ivan Hubeny and Thierry Lanz from Carlos Allende-Prieto's `synple <https://github.com/callendeprieto/synple>`_
package.  There is also a Python wrapper/driver based on the ``synple`` package and some code from Jon Holtzman
in the `APOGEE package <https://github.com/sdss/apogee>`_.  The setup.py file has also been modified to
automatically compile the Fortran code and copy them to the user's python scripts directory.  This reused code written by
Andy Casey in his `moog package <https://github.com/andycasey/moog>`_.

.. toctree::
   :maxdepth: 1

   install
   modules
	      

Description
===========
To run |synspec| you need 1) a model atmosphere, 2) a linelist (or multiple), and 3) the set of stellar parameters
and elemental abundances that you want to run.

1) Model Atmospheres

   Synspec can read TLUSTY, Kurucz/ATLAS, MARCS and PHOENIX model atmospheres.  See page 10 in the `Synspec manual <_static/syn43guide.pdf>`_ for the format and examples.

2) Linelists

   Synspec requires a specific linelist format.  See pages 9-10 in the `Synspec manual <_static/syn43guide.pdf>`_ for the format.
   
3) Stellar parameters and elemental abundances.

   The main stellar parameters are Teff, logg, [M/H], and [alpha/M].  These are the first four parameters in the
   main ``synthesis.synthesize()`` function.  The individual elements abundances are given in the ``elems`` parameters
   as a list of [element name, abundance] pairs, where abundance in the in [X/M] format relative to the overall metallicity.

More details on how the Fortran MOOG code works can be found in the `syn43guide.pdf <_static/syn43guide.pdf>`_ documentation.
   

Examples
========

.. toctree::
    :maxdepth: 1

    examples

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
