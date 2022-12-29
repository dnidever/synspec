import os
import numpy as np
import shutil
import subprocess
import tempfile
import time
from . import utils, atomic, atmos


def synthesize(teff,logg,mh=0.0,am=0.0,cm=0.0,nm=0.0,vmicro=2.0,elems=None,
               wrange=[15000.0,17000.0],dw=0.1,atmod=None,atmos_type='kurucz',
               dospherical=True,linelists=None,solarisotopes=False,workdir=None,
               save=False,verbose=False):
    """
    Code to synthesize a spectrum with MOOG.
    
    Parameters
    ----------
    teff : float
       Effective temperature in K.
    logg : float
       Surface gravity.
    mh : float, optional
       Metallicity, [M/H].  Deftauls is 0.0 (solar).
    am : float, optional
       Alpha abundance, [alpha/M].  Default is 0.0 (solar).
    cm : float, optional
       Carbon abundance, [C/M].  Default is 0.0 (solar).
    nm : float, optional
       Nitrogen abundance, [N/M].  Default is 0.0 (solar).
    vmicro : float, optional
       Microturbulence in km/s.  Default is 2 km/s.
    solarisotopes : bool, optional
       Use solar isotope ratios, else "giant" isotope ratios ( default False ).
    elems : list, optional
       List of [element name, abundance] pairs.
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step.  Default is 0.1 A.
    atmod : str, optional
       Name of atmosphere model (default=None, model is determined from input parameters).
    atmos_type : str, optional
       Type of model atmosphere file.  Default is 'kurucz'.
    dospherical : bool, optional
       Perform spherically-symmetric calculations (otherwise plane-parallel).  Default is True.
    linelists : list
       List of linelist file names.
    save : bool, optional
       Save temporary directory and files for synthesis.  Default=False.
    workdir : str, optional
       Directory to perform the work in.  By default a temporary directory is
         created and deleted after the work is done (unless save=True).
    verbose : bool, optional
       Verbose output to the screen.

    Returns
    -------
    flux : numpy array
       The fluxed synthetic spectrum.
    continuum : numpy array
       The continuum of the spectrum.
    wave : numpy array
       Wavelength array in A.

    Example
    -------

    flux,cont,wave = synthesize(5000.0,2.5,-1.0)

    """

    t0 = time.time()
    
    # Default abundances
    abundances = atomic.solar()
    abundances[2:] += mh
    abundances[6-1] += cm
    abundances[7-1] += nm
    for i in [8,10,12,14,16,18,20,22]:
        abundances[i-1] += am
    # Abundance overrides from els, given as [X/M]
    if elems is not None :
        for el in elems:
            atomic_num = atomic.periodic(el[0])
            abundances[atomic_num-1] = atomic.solar(el[0]) + mh + el[1]
    # Cap low abundances at -5.0
    #   that's what MOOG uses internally for the solar abundances
    for i in range(len(abundances)):
        if abundances[i]<-5:
            abundances[i] = -5.0
    # leave off the last one, 99, b/c for moog this means to scale ALL abundances by this value
    abundances = np.delete(abundances,98)
    
    # Change to temporary directory
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix='moog')
    cwd = os.getcwd()
    os.chdir(workdir)

    # Create the root name from the input parameters
    root = (atmos_type+'_t{:04d}g{:s}m{:s}a{:s}c{:s}n{:s}v{:s}').format(int(teff), atmos.cval(logg), 
                      atmos.cval(mh), atmos.cval(am), atmos.cval(cm), atmos.cval(nm),atmos.cval(vmicro))

    # Check that linelists and model atmosphere files exit
    if type(linelists) is str:
        linelists = [linelists]
    for l in linelists:
        if os.path.exists(l)==False:
            raise FileNotFoundError(l)
    if os.path.exists(atmod)==False:
        raise FileNotFoundError(atmod)

    if dospherical and ('marcs' in atmos_type) and logg <= 3.001:
        spherical= True
    else:
        spherical = False
    flux,cont,wave = do_moog(root,atmod,linelists,mh,am,abundances,wrange,dw,
                             save=save,solarisotopes=solarisotopes)

    os.chdir(cwd)
    if not save:
        shutil.rmtree(workdir)
    else:
        print('Saving temporary directory '+workdir)
        
    if verbose:
        print('dt = {:.3f}s'.format(time.time()-t0))
        
    return flux,cont,wave

    
def do_moog(root,atmod,linelists,mh,am,abundances,wrange,dw,
            solarisotopes=False,spherical=True,vmicro=2.0,
            molecules=None,save=False,verbose=False):
    """
    Runs MOOG for specified input parameters.

    Parameters
    ----------
    root : str
       Root of filenames to use for this MOOG run.
    atmod : str, optional
       Name of atmosphere model (default=None, model is determined from input parameters).
    linelists : list
       List of linelist file names.
    mh : float, optional
       Metallicity, [M/H].  Default is 0.0 (solar).
    am : float, optional
       Alpha abundance, [alpha/M].  Default is 0.0 (solar).
    abundances : list
       List of abundances in log epsilon format, log eps(X) = log(N(X)/N(H)) + 12.0.
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step.  Default is 0.1 A.
    solarisotopes : bool, optional
       Use solar isotope ratios, else "giant" isotope ratios.  Default is False.
    spherical : bool, optional
       Spherical atmosphere.  Default is True.
    vmicro : float, optional
       Microturbulent velocity in km/s.  Default is 2.0 km/s.
    molecules : list, optional
       List of molecules and ions to include in the molecular equilibrium calculations.
         By default, these ones are included:
         [606.0,106.0,607.0,608.0,107.0,108.0,112.0,707.0,708.0,
          808.0,12.1,60808.0,10108.0,101.0,60606.0,839.0,840.0,
          822.0,22.1,6.1,7.1,8.1,40.1,39.1]
    save : bool, optional
       Save temporary directory and files for synthesis.  Default=False.
    verbose : bool, optional
       Verbose output to the screen.

    Returns
    -------
    flux : numpy array
       The fluxed synthetic spectrum.
    continuum : numpy array
       The continuum of the spectrum.
    wave : numpy array
       Wavelength array in A.

    Example
    -------

    flux,cont,wave = do_moog(root,atmod,linefile,-0.1,0.2,abund,wrange=[15000.0,17000.0],dw=0.1)

    """

    # MOOG setup
    shutil.copy(atmod,'./'+os.path.basename(atmod))
    atmosfile = os.path.basename(atmod)
    alines = utils.readlines(atmosfile)
    
    # Prepare the end of the MOOG model atmosphere file 
    #NATOMS        3    0.10 
    #       6.0      7.70       8.0      8.30       7.0      7.50 
    #NMOL         21 
    #     606.0     106.0     607.0     608.0     107.0     108.0     112.0 
    #707.0 
    #     708.0     808.0      12.1   60808.0   10108.0     101.0   60606.0 
    #839.0 
    #     840.0     822.0      22.1      40.1      39.1 

    # MOOG internally stores the Asplund+2009 solar abundances
    # see Batom.f
    
    # Create MOOG control file
    addon = []
    natoms = len(abundances)
    addon.append('NATOMS        {0:d}    {1:.3f}'.format(natoms,mh))
     
    # CONVERT TO LOG EPSILON FORMAT 
    # log eps(X) = log(N(X)/N(H)) + 12.0 
    # [X/H] = log(n(X)/n(H)) - log(n(X)/n(H))_solar 
    # log(n(X)/n(H)) = [X/H] + log(n(X)/n(H))_solar 
    # log eps(X) = [X/H] + log(n(X)/n(H))_solar + 12.0 
    # log eps(X) = [X/Fe] + [Fe/H] + (log(n(X)/n(H))_solar + 12.0) 
    #INPUT ABUNDANCES: (log10 number densities, log H=12) 
    #      Default solar abundances: Anders and Grevesse 1989 
    # This is the (log(n(X)/n(H))_solar + 12.0) value 
    # from Batom.f 
    # Mg(12)= 7.58 
    # logeps_mg = alpha + metal + 7.58 
    # addon.append('      {0:4.1f}      {1:.3f}'.format(12.0,logeps_mg))  # Mg
    # The abundances have to be input as log epsilon
    for iel,abun in enumerate(abundances):
        addon.append("      {0:4.1f}      {1:8.3f}".format(iel+1,abun))
    # Molecules and ions to be included in the molecular equilibrium calculation
    #   the atoms will be automatically added
    # 6-C, 7-N, 8-O
    # 606 - C2
    # 106 - CH
    # 607 - CN
    # 608 - CO
    # 107 - NH
    # 108 - OH
    # 112 - MgH
    # 707 - N2
    # 708 - NH
    # 808 - O2
    # 21.1 - Sc II
    # 60808.0 - CO2
    # 10108.0 - H20
    # 101.0 - H2
    # 60606.0 - C3
    # 839.0 - YO
    # 840.0 - ZrO
    # 822.0 - TiO
    # 22.1 - Ti II
    # 40.1 - Zr II
    # 39.1 - Y II
    if molecules is not None:
        mols = molecules
    else:
        mols = [606.0,106.0,607.0,608.0,107.0,108.0,112.0,707.0,708.0,
                808.0,12.1,60808.0,10108.0,101.0,60606.0,839.0,840.0,
                822.0,22.1,6.1,7.1,8.1,40.1,39.1]
    addon.append('NMOL         '+str(len(mols)))
    mollines = ''
    for l in mols:
        mollines += '{:10.1f}'.format(float(l))
    # Put 7 per line
    for i in range(int(np.ceil(len(mols)/7))):
        addon.append(mollines[i*70:i*70+70])
    newalines = alines + addon
    utils.writelines(atmosfile,newalines,overwrite=True)
    
     
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 
    # RUN MOOG 
     
    # Wavelength parameters 
    w0 = wrange[0]
    w1 = wrange[1]
    fwhm = 0.01  # Gaussian broadening 
     
    # Read in the linelists
    lines = []
    for i in range(len(linelists)):
        lines1 = utils.readlines(linelists[i],noblank=True)  #,comment='#')
        lines += lines1
    nlinelist = len(lines)
    # Get wavelengths
    lwave = np.array([float(l.split()[0]) for l in lines]).astype(float)
    # Sort the lines if more than one linelist was input
    if len(linelists)>1:
        si = np.argsort(lwave)
        lwave = lwave[si]
        lines = np.char.array(lines)[si]
    else:
        lines = np.char.array(lines)        
    wavemin = np.min(lwave) 
    wavemax = np.max(lwave) 
    
    # Make temporary linelist file 
    tid,templist = tempfile.mkstemp(prefix='line')
    templist = os.path.basename(templist)
    gd, = np.where((lwave >= (w0-1)) & (lwave <= (w1+1)))  # allow 1A buffer on ends
    tlines = np.char.array(lines)[gd]
    utils.writelines(templist,tlines)
         
    # This is what the MOOG input file looks like 
    #terminal       'xterm' 
    #model_in       'atmos' 
    #lines_in       'vald.dat' 
    #atmosphere    1 
    #molecules     2 
    #lines         1 
    #flux/int      0 
    #damping       0 
    #abundances    1   1 
    #        3     -0.87 
    #synlimits 
    #  5085.0  5320.0    0.02    1.00 
    #plotpars      1 
    #  5085.0  5320.0    0.00    1.00 
    #   0.0       0.0    0.00   1.0 
    #    g      1.0      0.0    0.00     0.0     0.0 
    #opacit        0
    
    # Set up the parameter file 
    params = []
    params.append("synth")
    params.append("terminal       'xterm'")
    params.append("standard_out   out1")
    params.append("summary_out    out2")
    params.append("smoothed_out   out3")
    params.append("model_in       '{:s}'".format(atmosfile))
    params.append("lines_in       '{:s}'".format(templist))
    params.append("atmosphere    1")
    params.append("molecules     2")
    params.append("lines         1")
    params.append("flux/int      0")
    params.append("damping       0")
    #params.append("abundances    1   1" 
    #params.append("        3     -0.87" 
    #***!!!!! NOT REQUIRED IF INPUT IN MODEL ATMOSPHERE!!!!! *** 
    # ADD ABUNDANCES OF ALPHA ELEMENTS 
    # For spectrum use these alpha elements: Mg, Si, S, Ar, Ca and Ti. 
    # Mg-12, Si-14, S-16, Ar-18, Ca-20, Ti-22 
    #if (float(alpha) ne 0.0) then begin 
    #  PUSH,tpaarms,'abundances    6   1' 
    #  params.append('       12     '+strtrim(string(logeps_mg,format='(F8.2)'),2) ; Mg 
    #  params.append('       14     '+strtrim(string(logeps_si,format='(F8.2)'),2) ; Si 
    #  params.append('       16     '+strtrim(string(logeps_s,format='(F8.2)'),2) ; S 
    #  params.append('       18     '+strtrim(string(logeps_ar,format='(F8.2)'),2) ; Ar 
    #  params.append('       20     '+strtrim(string(logeps_ca,format='(F8.2)'),2) ; Ca 
    #  params.append('       22     '+strtrim(string(logeps_ti,format='(F8.2)'),2) ; Ti 
    #end 
    params.append("synlimits")
    params.append("  {0:10.3f}  {1:10.3f}  {2:10.3f}  1.00".format(w0,w1,dw))
    params.append("plotpars      1")
    params.append("  {0:10.3f}  {1:10.3f}   0.0   1.00".format(w0,w1))
    params.append("   0.0       0.0    0.00   1.0")
    params.append("    g     {:.3f}      0.0    0.00     0.0     0.0   ".format(fwhm))
    params.append("opacit        0")
    # MOOGSILENT is hardwired to use "batch.par"
    utils.writelines('batch.par',params)
    
    # Run MOOGSILENT
    ret = subprocess.check_output(['MOOGSILENT'],stderr=subprocess.STDOUT)    
    # Save the log file
    if type(ret) is bytes:
        ret = ret.decode()
    with open(root+'_MOOG.log','w') as f:
        f.write(ret)

    # Load the spectrum
    wave,flux = utils.read_synthfile('out2')
    cont = np.zeros(flux.shape,float)
    
    return flux,cont,wave
