************
Installation
************


Important packages
==================
`synspec` is a package to run the `Synspec <http://tlusty.oca.eu/Synspec49/synspec.html>`_
synthesis software by Ivan Hubeny and Thierry Lanz from Carlos Allende-Prieto's `synple <https://github.com/callendeprieto/synple>`_
package. There is also a Python wrapper/driver based on the ``synple`` package and some code from Jon Holtzman in the
in the `APOGEE package <https://github.com/sdss/apogee>`_.

Installing Synspec
==================

The easiest way to install the code is with pip.  This will both compile and install the Fortran code as
well as the Python wrapper code.

.. code-block:: bash

    pip install synspec
 
*NOTE* that this will take a few minutes because it will download large linelist files and convert them from gzipped ASCII to binary format.
    
Fortran code
------------
    
The pip install will attempt to automatically compile the Fortran code and copy the binaries to your
Python scripts directory (which should be in your path).  It also attemps to download linelists.
If this fails for some reason, then you'll need to compile it yourself.  You'll likely want to do a
full git clone of the repository for this.
To compile the code you need the GNU Fortran compiler (``gfortran``).
The Fortran code lives in the `src/` directory.  All you should need to do is to cd into that
directory and type ``make``.
Copy the binaries ``synspec54``, ``rotin``, and ``list2bin`` to a directory in your path (e.g., ~/bin/ or /usr/local/bin/).  

Dependencies
============

- numpy
- scipy
- astropy
- matplotlib
- `dlnpyutils <https://github.com/dnidever/dlnpyutils>`_
