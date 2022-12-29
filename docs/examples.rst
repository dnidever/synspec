********
Examples
********


Running Synspec
===============
You need to supply |synspec| with 1) a model atmosphere, 2) a linelist (or multiple linelists), and 3) the set of stellar
parameters and abundances.

Here's a simple example:

    >>> from synspec import synthesis,utils,models
    >>> linelists = [utils.linelistsdir()+'gfATO.19.11',utils.linelistsdir()+'gfMOLsun.20.11']
    >>> atmod = utils.modelsdir()+'ksun.mod'
    >>> flux,cont,wave = synthesis.synthesize(5777,4.44,0.00,linelists=linelists,atmod=atmod,wrange=[6700,6800])

Now plot the spectrum:

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(wave,flux)

It should look like this.

.. image:: spectrum_example.png
  :width: 600
  :alt: Example Synspec synthetic spectrum

Abundances
----------
	
You can modify the global alpha abundance with `am` or individual abundances with `elems`.  The `elems` parameter
takes a list of [element name, abundance] pairs, where the abundance should be in the form [X/M], where M is the
overall metallicity that is used to scale the individual abundances.  For example, ``elems=[['Mg',0.55],['Ba',-0.15]]``
means a relative Magnesium abundance of +0.55 and a relative Barium abundance of -0.15.

Let's try it out:

    >>> flux2,cont2,wave2 = synthesis.synthesize(5777,4.44,0.00,linelists=linelists,atmod=atmod,wrange=[6700,6800],elems=[['Mg',0.55],['Ba',-0.15]])
    >>> plt.plot(wave,flux)
    >>> plt.plot(wave2,flux2)
    >>> plt.xlim(6700,6800)
