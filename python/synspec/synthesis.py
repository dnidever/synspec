import os
import numpy as np
import shutil
import subprocess
import tempfile
import time
from . import utils, atomic, atmos, models

clight = 299792.458
epsilon = 0.6 #clv coeff.
bolk = 1.38054e-16  # erg/ K
zero = " 0 "
one =  " 1 "
two =  " 2 "

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

    # Default linelists
    if linelists is None:
        linelistdir = utils.linelistsdir()
        linelists = ['gfATO.19.11','gfMOLsun.20.11','gfTiO.20.11','H2O-8.20.11']
        linelists = [l+linelistdir for l in linelists]
        
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
    # Synspec expects N(X)/N(H), while this is log eps(X) = log(N(X)/N(H)) + 12.0.
    abundances = 10**(abundances-12.0)
    
    # Change to temporary directory
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix='synspec')
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
    flux,cont,wave = do_synspec(root,atmod,linelists,mh,am,abundances,wrange,dw,
                                solarisotopes=solarisotopes)

    os.chdir(cwd)
    if not save:
        shutil.rmtree(workdir)
    else:
        print('Saving temporary directory '+workdir)
        
    if verbose:
        print('dt = {:.3f}s'.format(time.time()-t0))
        
    return flux,cont,wave


def do_synspec(root,atmod,linelists,mh,am,abundances,wrange,dw=None,
               solarisotopes=False,spherical=True,vmicro=2.0,vrot=0.0,
               fwhm=0.0,vmacro=0.0,atom='ap18',strength=1e-4,
               lte=None,verbose=False):
    """
    Runs Synspec for specified input parameters.

    Interface to the fortran codes synspec/rotin that only requires two mandatory inputs: 
    a model atmosphere (modelfile) and the limits of the spectral range (wrange). The code 
    recognizes Kurucz, MARCS and Phoenix LTE model atmospheres. The sampling of the frequency 
    grid is chosen internally, but can also be set by adding a constant wavelength step (dw).
    The abundances and microturbulence velocity can be set through the abu and vmicro 
    parameters, but default values will be taken from the model atmosphere. Rotational and 
    Gaussian broadening can be introduced (vrot and fwhm parameters), as well as radial-tangential
    macroturbulence.

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
    abundances : list/array, optional
       Chemical abundances relative to hydrogen (N(X)/N(H)). 99 values.
       (default taken from input model atmosphere)
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step for the output fluxes, this will trigger interpolation at the end.
        A negative value will lead to interpolation to a uniform step in ln(lambda)
        (default is None for automatic frequency selection).
    solarisotopes : bool, optional
       Use solar isotope ratios, else "giant" isotope ratios.  Default is False.
    spherical : bool, optional
       Spherical atmosphere.  Default is True.
    vmicro : float, optional
       Microturbulent velocity in km/s.  Default is 2.0 km/s.
    vrot: float, optional
       Projected rotational velocity (km/s).  Default 0.
    fwhm: float, optional
       Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
         Ddefault 0.
    vmacro : float, optional
       Radial-tangential macroturbulence (km/s).  Default 0.
    atom : str, optional
       Type of opacities to use:
        'ap18' -- generic opacities used in Allende Prieto+ 2018
        'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
        'hhm' -- continuum opacity is simplified to H and H-
       Default is 'ap18'.
    strength: float, optional
       Threshold in the line-to-continuum opacity ratio for 
         selecting lines (default is 1e-4).
    lte: bool, optional
       This flag can be set to True to enforce LTE in NLTE models. MARCS, Kurucz, the 
        class of Phoenix models used here are always LTE models. Tlusty models
        can be LTE or NLTE, and this keyword will ignore the populations and compute
        assuming LTE for a input NLTE Tlusty model.
        (default None)
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

    # Synspec setup
    shutil.copy(atmod,'./'+os.path.basename(atmod))
    modelfile = os.path.basename(atmod)
    
    # Basic checks on the line list and model atmosphere
    #linelist, modelfile = checksynspec(linelist,modelfile,verbose=verbose)

    # Read model atmosphere
    atmostype, teff, logg, vmicro2, abu2, nd, atmos = models.read_model(modelfile,verbose=verbose)

    if vmicro is None: vmicro = vmicro2
    if abundances is None: abundances = abu2
    # We take a step of 1/3 of the Gaussian (thermal + micro) FWHM at the lowest T and for an atomic mass of 100
    space = np.mean(wrange) / clight * 2.355 / 3. * np.sqrt(0.1289**2 * np.min(atmos['t']) / 100. + vmicro** 2 / 2.) 
    
    # Check input parameters are valid
    imode = checkinput(wrange, vmicro, linelists)
  
    # Find out additional info for Tlusty models
    if atmostype == 'tlusty':
        madaffile, nonstdfile, nonstd, numpar, datadir, inlte, atommode, atominfo = read_tlusty_extras(modelfile)
        if (inlte == -1): nonstd['IBFAC'] = 1
        if (lte == True): inlte = 0
    else:
        nonstd = None
        inlte = 0
        atommode = None
        atominfo = None
 
    if verbose:
        print(modelfile,'is a',atmostype,' model')
        print('teff,logg,vmicro=',teff,logg,vmicro)
        print ('abundances=',abundances)
        print ('linelist=',linelist)
        print ('wrange=',wrange)

    logfile = 'syn.log'
    
    # Link data folder for used-provided model atoms
    dd = ''
    if atmostype == 'tlusty':
        # data dir
        hdd, dd = os.path.split(datadir)
        os.symlink(datadir,dd)
        if dd == 'data':
            for entry in isdf:
                assert (os.path.isfile(os.path.join(dd,entry))), 'Cannot find the data file:'+dd+'/'+entry     

    # Link the data folder (synspec data + default tlusty model atoms)
    modelatomdir = utils.datadir()    
    if dd != 'data':
        os.symlink(modelatomdir,'./data')          # data directory  

    write5(teff,logg,abundances,atom,inlte=inlte,atommode=atommode,atominfo=atominfo)   # abundance/opacity file  
    write8(teff,logg,nd,atmos,atmostype)                    # model atmosphere

    iprin = 0
    cutoff0 = 250.
    if logg > 3.0:
        if teff < 4000.: cutoff0 = 500.
        if teff < 3500.: cutoff0 = 1000.
        if teff < 3000.: cutoff0 = 1500. 
    write55(wrange,dw=space,imode=imode,iprin=iprin, inlte=inlte, hydprf=2,
            cutoff0=cutoff0,strength=strength,vmicro=vmicro,
            linelist=linelists,atmostype=atmostype)
    # Synspec control file
    writetas('tas',nd,linelists,nonstd=nonstd)               # non-std param. file
    create_links(linelists)                                  # auxiliary data

    synin = open('fort.5')
    synout = open(logfile,'w')

    p = subprocess.Popen(['synspec54'], stdin=synin, stdout = synout, stderr= synout, shell=True)
    p.wait()

    synout.flush()
    synout.close()
    synin.close()

    assert (os.path.isfile('fort.7')), 'Error: I cannot read the file *fort.7* in '+tmpdir+' -- looks like synspec has crashed, please look at '+logfile

    assert (os.path.isfile('fort.17')), 'Error: I cannot read the file *fort.17* in '+tmpdir+' -- looks like synspec has crashed, please look at '+logfile
    
    wave, flux = np.loadtxt('fort.7', unpack=True)
    if np.any(np.diff(wave) <= 0.0):
        wave, win = np.unique(wave,return_index=True)
        flux = flux[win] 
    wave2, flux2 = np.loadtxt('fort.17', unpack=True)  # low-res continuum
    if np.any(np.diff(wave2) <= 0.0):
        wave2,win = np.unique(wave2,return_index=True)
        flux2 = flux2[win]
    if dw is None and fwhm <= 0. and vrot <= 0.:
        cont = np.interp(wave, wave2, flux2)


    if fwhm > 0. or vrot > 0. or vmacro > 0.:
        print(vrot, fwhm, vmacro, space, steprot, stepfwhm)
        wave, flux = call_rotin(wave, flux, vrot, fwhm, vmacro,
                                space, 0.0, 0.0, clean=False,
                                reuseinputfiles=True, logfile=logfile)
        if dw is None:
            cont = np.interp(wave, wave2, flux2)

    if (dw is not None): 
        if (dw < 0.):
            ldw = np.abs(dw)/np.mean(wrange)
            nsamples = int((np.log(wrange[1]) - np.log(wrange[0]))/ldw) + 1
            wave3 = np.exp(np.arange(nsamples)*ldw + np.log(wrange[0]))
        else:
            nsamples = int((wrange[1] - wrange[0])/dw) + 1
            wave3 = np.arange(nsamples)*dw + wrange[0]
        cont = np.interp(wave3, wave2, flux2)
        flux = np.interp(wave3, wave, flux)
        wave = wave3


    #if save == True:
    #    out = ['MODEL   = '+modelfile+'\n']
    #    out.append('TEFF    = '+str(teff)+'\n')
    #    out.append('LOGG    = '+str(logg)+'\n')
    #    out.append('VMICRO  = '+str(vmicro)+'\n')
    #    out.append('WRANGE  = '+' '.join(map(str,wrange))+'\n')
    #    out.append('STRENGTH= '+str(strength)+'\n')
    #    out.append('LINELIST= '+' '.join(linelist)+'\n')
    #    out.append('ATOM    = '+atom+'\n')
    #    out.append('VROT    = '+str(vrot)+'\n')
    #    out.append('FWHM    = '+str(fwhm)+'\n')
    #    out.append('VMACRO    = '+str(vmacro)+'\n')
    #    out.append('STEPROT = '+str(steprot)+'\n')
    #    out.append('STEPFWHM= '+str(stepfwhm)+'\n')
    #    out.append('LTE     = '+str(lte)+'\n')
    #    out.append('ABU     = '+' '.join(map(str,abu))+'\n')
    #    header = ''.join(out)
    #    if synfile == None: 
    #        tmpstr = os.path.split(modelfile)[-1]
    #        synfile = tmpstr[:tmpstr.rfind('.')]+'.syn'
    #    np.savetxt(synfile,(wave,flux,cont),header=header)

    return flux, cont, wave


def write55(wrange,dw=1e-2,imode=0,iprin=0,inlte=0,hydprf=2,cutoff0=200.,
            strength=1e-4,vmicro=0.0,linelist=None, atmostype='kurucz'):

    # imode,idst,iprin
    # inmod,zero,ichang,ichemc
    # lyman,zero,zero,zero,zero
    # one,nlte,icontl,inlist,ifhe2
    # ihydpr,ihe1pr,ihe2pr
    # wstart,wend,cutoff,zero,strength,wdist 

    if (atmostype == 'tlusty' or atmostype == 'marcs'): inmod = 1 
    else: inmod = 0

    all_inlist = []
    for file in linelist:
        inlist = 10
        if file[-3:] == '.11' : inlist = 11
        all_inlist.append(inlist)
        assert (inlist - all_inlist[0] == 0), 'The line list files must be all either text or binary!'

    f = open('fort.55','w')
    f.write(" "+str(imode)+" "+zero+" "+str(iprin)+"\n")
    f.write(" "+str(inmod)+3*zero+"\n")
    f.write(5*zero+"\n")
    f.write(one+str(abs(inlte))+zero+str(inlist)+zero+"\n")
    f.write(str(hydprf)+2*zero+"\n")
    if imode <= -3:
        f.write( ' %f %f %f %i %e %f \n ' % (wrange[0],  -wrange[1], cutoff0, 0, strength, dw) )
    else:
        f.write( ' %f %f %f %i %e %f \n ' % (wrange[0],   wrange[1], cutoff0, 0, strength, dw) )
    ll = len(linelist)
    if ll < 2: f.write(2*zero)
    else: f.write(str(ll-1) + ' ' + ' '.join(map(str,np.arange(ll-1)+20)))
    f.write("\n")
    f.write( ' %f  \n' % (vmicro) )
    f.close()


def write5(teff,logg,abu, atom='ap18', ofile='fort.5', inlte=0, atommode=None, atominfo=None):
    
    symbol, mass, sol = atomic.elements()

    f = open(ofile,'w')
    f.write(' '+str(teff)+" "+str(logg).format('%7.4f')+"       ! TEFF, GRAV \n")
    if inlte == 0:
        f.write(" T  F               ! LTE, GRAY \n")
    else:
        f.write(" F  F               ! LTE, GRAY \n")

    f.write(" 'tas'              ! name of non-standard flags \n")
    f.write(" 50                 ! frequencies \n")

    natom = len(abu)
    f.write(" "+str(natom)+"        ! NATOMS \n")  

    assert (atom == 'hhm' or atom == 'ap18' or atom == 'yo19' or atom == 'test'), 'atom must be one of: hhm/ap18/yo19/test!'
    ex = np.ones(natom)
    if atom == 'hhm' : 
        zex = [1]  # atomic numbers of elements included explicitly (contributing cont. opacity)
    elif atom == 'yo19':
        zex = [1,11,12,19,20]
    elif atom == 'test': 
        zex = [1,26]
    elif atom == 'ap18': 
        zex = [1,2,6,7,8,11,12,13,14,20,26]

    for i in zex: ex[i-1] = 2

    # tlusty models provide atommode and atominfo that override the defaults for atom and ex
    if atommode is not None: ex = np.array(atommode)
    if atominfo is not None: atom = 'own'

    for i in range(natom):
        f.write(' %2d %e %i %s\n' %  (ex[i], abu[i], 0, '  ! ' +symbol[i]) )

    for i in range(3): f.write("* \n")
  
    if atom == 'hhm':  # highly simplified continuum opacities -- just H and H-
        f.write("* ../data_atom for ions  \n")
        f.write("   1   -1     1      0     0     1    ' H 1' 'data/hm.dat' \n" )
        f.write("   0    0     3      0 \n")
        f.write("   1    0     9      0     0     0    ' H 1' 'data/h1s.dat'  \n")
        f.write("   1    1     1      1     0     0    ' H 2' ' '  \n")
        f.write("   0    0     0     -1     0     0    '    ' ' '  \n")
        f.write("* \n")
        f.write("* end \n")
    elif atom == "yo19": # set for NLTE calculations for APOGEE (see Osorio+ 2019 A&A paper)
        f.write("* ../data_atom for ions  \n")
        f.write("  1    -1     1      0     0     1    ' H 0' 'data/hm.dat'  \n")
        f.write("  0     0     3      0   \n")
        f.write("  1     0     16     0     0     0    ' H 1' 'data/h1_16lev2.dat'  \n")
        f.write("  1     1     1      1     0     0    ' H 2' ' '  \n")
        f.write("  11    0     42     0     0     0    'Na 1' 'data/NaIkas.tl'  \n")
        f.write("  11    1     1      1     0     0    'Na 2' '' \n")
        f.write("  12    0     96     0     0     0    'Mg 1' 'data/Mg1kas_F_ccc.sy'  \n")
        f.write("  12    1     29     0     0     0    'Mg 2' 'data/Mg2kas_F_ccc.sy'  \n")
        f.write("  12    2     1      1     0     0    'Mg 3' ' '  \n")
        f.write("  19    0     31     0     0     0    'K  1' 'data/KIkas.tl'  \n")
        f.write("  19    1     1      1     0     0    'K  2' ''  \n")
        f.write("  20    0     66     0     0     0    'Ca 1' 'data/Ca1kas_F_zat.sy'  \n")
        f.write("  20    1     24     0     0     0    'Ca 2' 'data/Ca2kas_F_zat.sy'  \n")
        f.write("  20    2     1      1     0     0    'Ca 3' ' '  \n")
        f.write("   0    0     0     -1     0     0    '    ' ' '  \n")
        f.write("* \n")
        f.write("* end \n")
    elif atom == 'test': # generic set used in Allende Prieto+ (2018) A&A paper
        f.write("* ../data for ions  \n")
        f.write("   1   -1     1      0     0     1    ' H 1' 'data/hm.dat'  \n")
        f.write("   0    0     3      0 \n")
        f.write("   1    0     9      0     0     0    ' H 1' 'data/h1s.dat'  \n")
        f.write("   1    1     1      1     0     0    ' H 2' ' '  \n")
        f.write("  26    0     49     0     0     0    'Fe 1' 'data/tlusty_fe1_topmod.dat'  \n")
        f.write("  26    1     41     0     0     0    'Fe 2' 'data/tlusty_fe2_topmod.dat'  \n")
        f.write("  26    2     1      1     0     0    'Fe 3' ' '  \n")
        f.write("  0     0     0     -1     0     0    '    ' ' '  \n")
        f.write("* \n")
        f.write("* end \n")
    elif atom == 'ap18': # generic set used in Allende Prieto+ (2018) A&A paper
        f.write("* ../data for ions  \n")
        f.write("   1   -1     1      0     0     1    ' H 1' 'data/hm.dat'  \n")
        f.write("   0    0     3      0 \n")
        f.write("   1    0     9      0     0     0    ' H 1' 'data/h1s.dat'  \n")
        f.write("   1    1     1      1     0     0    ' H 2' ' '  \n")
        f.write("   2    0     14     0     0     0    'He 1' 'data/he1.dat'  \n")
        f.write("   2    1     14     0     0     0    'He 2' 'data/he2.dat '  \n")
        f.write("   2    2     1      1     0     0    'He 3' ' '  \n")
        f.write("   6    0     104    0     0     0    ' C 1' 'data/c1.t'  \n")
        f.write("   6    1     40     0     0     0    ' C 2' 'data/c2.t'  \n")
        f.write("   6    2     1      1     0     0    ' C 3' ' '  \n")
        f.write("   7    0     89     0     0     0    ' N 1' 'data/n1.t'  \n")
        f.write("   7    1     51     0     0     0    ' N 2' 'data/n2.t'  \n")
        f.write("   7    2     1      1     0     0    ' N 3' ' '  \n")
        f.write("   8    0     54     0      0    0    ' O 1' 'data/o1.t'  \n")
        f.write("   8    1     74     0      0    0    ' O 2' 'data/o2.t'  \n")
        f.write("   8    2     1      1      0    0    ' O 3' ' '  \n")
        f.write("  11    0     32     0     0     0    'Na 1' 'data/na1.t'  \n")
        f.write("  11    1     8      0     0     0    'Na 2' 'data/na2.t'  \n")
        f.write("  11    2     1      1     0     0    'Na 3' ' '  \n")
        f.write("  12    0     71     0     0     0    'Mg 1' 'data/mg1.t'  \n")
        f.write("  12    1     31     0     0     0    'Mg 2' 'data/mg2.t'  \n")
        f.write("  12    2     1      1     0     0    'Mg 3' ' '  \n")
        f.write("  13    0     33     0     0     0    'Al 1' 'data/al1.t'  \n")
        f.write("  13    1     81     0     0     0    'Al 2' 'data/al2.t'  \n")
        f.write("  13    2     1      1     0     0    'Al 3' ' '  \n")
        f.write("  14    0     57     0     0     0    'Si 1' 'data/si1.t'  \n")
        f.write("  14    1     46     0     0     0    'Si 2' 'data/si2.t'  \n")
        f.write("  14    2     1      1     0     0    'Si 3' ' '  \n")
        f.write("  20    0     79     0     0     0    'Ca 1' 'data/ca1.t'  \n")
        f.write("  20    1     32     0     0     0    'Ca 2' 'data/ca2.t'  \n")
        f.write("  20    2     1      1     0     0    'Ca 3' ' '  \n")
        f.write("  26    0     49     0     0     0    'Fe 1' 'data/tlusty_fe1_topmod.dat'  \n")
        f.write("  26    1     41     0     0     0    'Fe 2' 'data/tlusty_fe2_topmod.dat'  \n")
        f.write("  26    2     1      1     0     0    'Fe 3' ' '  \n")
        f.write("  0     0     0     -1     0     0    '    ' ' '  \n")
        f.write("* \n")
        f.write("* end \n")
    else:
        for line in atominfo: f.write(line)
    f.close()


def write8(teff, logg, nd, atmos, atmostype, ofile='fort.8'):
    """
    Writes the model atmosphere for synspec

    MARCS models can be passed in 'Tlusty' (default, after read with 
           read_marcs_models2) or 'Kurucz' format
    Phoenix and Kurucz models are passed to synspec formatted as 'Kurucz'

    """

    f = open(ofile,'w')

    if atmostype == 'tlusty':

        if ('n' in atmos.dtype.names):  # 4th column is number density n
            if ('pop' in atmos.dtype.names):   # explicit (usually NLTE) populations
                numpop = len(atmos['pop'][0]) 
                sformat = '  %f %e %e %e'
                i = 5
                for entry in atmos['pop'][0]: 
                    sformat = sformat + ' %e'
                    if i % 6 == 0: sformat = sformat + '  \n'
                    i = i + 1
                sformat = sformat + ' \n' 
                f.write(" "+str(nd)+" "+str(-(4+numpop))+"\n")
                for i in range(nd):
                    f.write(' %e ' % atmos['dm'][i])
                    if (i+1) % 5 == 0: f.write('\n')
                if (i+1) % 5 != 0: f.write('\n')
                for i in range(nd):
                    sdata = [atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['n'][i]]
                    for j in range(numpop):
                        sdata.append(atmos['pop'][i][j])
                    f.write( sformat % tuple(sdata) )                 
            elif ('dep' in atmos.dtype.names): # NLTE departure coefficients
                numpop = len(atmos['dep'][0]) 
                sformat = '  %f %e %e %e'
                i = 5
                for entry in atmos['dep'][0]: 
                    sformat = sformat + ' %e'
                    if i % 6 == 0: sformat = sformat + '  \n'
                    i = i + 1
                sformat = sformat + ' \n' 
                f.write(" "+str(nd)+" "+str(-(4+numpop))+"\n")
                for i in range(nd):
                    f.write(' %e ' % atmos['dm'][i])
                    if (i+1) % 5 == 0: f.write('\n')
                if (i+1) % 5 != 0: f.write('\n')
                for i in range(nd):
                    sdata = [atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['n'][i]]
                    for j in range(numpop):
                        sdata.append(atmos['dep'][i][j])
                    f.write( sformat % tuple(sdata) )         
            else:                              # LTE
                f.write(" "+str(nd)+" "+str(-4)+"\n")
                for i in range(nd):
                    f.write(' %e ' % atmos['dm'][i])
                    if (i+1) % 5 == 0: f.write('\n')
                if (i+1) % 5 != 0: f.write('\n')
                for i in range(nd):
                    f.write( '%f %e %e %e \n' % (atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['n'][i] ) )  
        else:
            pass

    else:

        if atmostype == 'marcs':
            f.write(" "+str(nd)+" "+str(-4)+"\n")
            for i in range(nd):
                f.write(' %e ' % atmos['dm'][i])
            f.write("\n")
            for i in range(nd):
                f.write( '%f %e %e %e \n' % (atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['rho'][i]/atmos['mmw'][i]/1.67333e-24 + atmos['ne'][i] ) )

        else:
            f.write( 'TEFF %7.0f  GRAVITY %7.5f  LTE \n' % (teff, logg) )
            for i in range(21): f.write('\n')
            f.write( 'READ DECK6%3i RHOX,T,P,XNE \n' % nd )
            for i in range(nd): 
                f.write( '%e %f %e %e \n' % (atmos['dm'][i], atmos['t'][i], atmos['p'][i], atmos['ne'][i]) )
      
    f.close()

    return()


def call_rotin(wave=None, flux=None, vrot=0.0, fwhm=0.0, vmacro=0.0, 
               space=1e-2, steprot=0.0, stepfwhm=0.0, 
               clean=True, reuseinputfiles=False, logfile='syn.log'):
    """
    Convolves a synthetic spectrum with a rotation and/or Gaussian kernel

    Interface to the fortran code rotin.

    Parameters
    ----------
    wave: numpy array of floats
      wavelengths (angstroms)
    flux: numpy array of floats
      flux 
    vrot: float
      projected rotational velocity (km/s)
      (default 0.)
    fwhm: float
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
    vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)
    space: float, optional
      characteristic wavelength scale for variations in the spectrum (angstroms)
      (default is 1e-2)
    steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
    stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
    clean: bool
      True by the default, set to False to avoid the removal of the rotin
      temporary files (default Tr<ue)
    reuseinputfiles: bool
      set to take the input data from the output synspec file (fort.7) rather than 
      from the input arrays (wave, flux)
    logfile: str
      name of the log file
      (default syn.log)

    Returns
    -------
    wave2: numpy array of floats
      wavelengths (angstroms)
    flux2: numpy array of floats
      flux 

    """
    if reuseinputfiles == False:
        f = open('fort.7','w')
        f2 = open('fort.17','w')
        maxflux = np.max(flux)
        for i in range(len(wave)):
            f.write( ' %f %f \n' % (wave[i], flux[i]) )
            f2.write( ' %f %f \n' % (wave[i], maxflux) )
        f.close()
        f2.close()

    f = open('fort.5','w')
    f.write( ' %s %s %s \n' % ("'fort.7'", "'fort.17'", "'fort.11'") )
    f.write( ' %f %f %f \n' % (vrot, space, steprot) )
    f.write( ' %f %f %f \n' % (fwhm, stepfwhm, vmacro) )
    print('stepfwhm=',stepfwhm)
    f.write( ' %f %f %i \n' % (np.min(wave), np.max(wave), 0) )
    f.close()

    synin = open('fort.5')
    synout = open(logfile,'a')
    p = subprocess.Popen([rotin], stdin=synin, stdout = synout, stderr = synout)
    p.wait()
    synout.flush()
    synout.close()
    synin.close()
  
    assert (os.path.isfile('fort.11')), 'Error: I cannot read the file *fort.11* in '+os.getcwd()+' -- looks like rotin has crashed, please look at syn.log'

    wave2, flux2 = np.loadtxt('fort.11', unpack=True)
    print(len(wave),len(wave2))
  
    if clean == True: cleanup_fort()

    return(wave2, flux2)


def checksynspec(linelist,modelfile,verbose=False):
    """
    Checking that executables and data are where it should be. Prepend
    default directories to linelist and model atmosphere files when necessary

    Parameters
    ----------
    linelist: array of str
      file names of the line lists to be used. The first string should correspond
      to the atomic line list and is mandatory. The remainder are optional and
      correspond to molecular line lists. All files should be in synspec format.
      (see documentation at http://nova.astro.umd.edu/Synspec43/synspec.html)

    """

    #dirs = [synpledir,modelatomdir,linelistdir,bindir]
    #for entry in dirs: assert (os.path.isdir(entry)), 'dir '+entry+' missing'

    linelist = checklinelistpath(linelist)

    files = [synspec,rotin]
    for entry in linelist: files.append(entry)
    for entry in files: assert (os.path.isfile(entry)), 'file '+entry+' missing'

    if not os.path.isfile(modelfile):
        mf = os.path.join(modeldir,modelfile)
        if os.path.isfile(mf): modelfile = mf

    if verbose:
        print(modeldir)
        print(modelfile)
    assert (os.path.isfile(modelfile)),'model atmosphere file '+modelfile+' missing'

    return(linelist,modelfile)

def checklinelistpath(linelist):
    """
    Checking whether the line lists in the array linelists are in the current folder
    or in linelistdir and returning the same array with an absolute path
    """

    i = 0 
    for entry in linelist: 
        if os.path.isfile(entry+".11"): 
            entry = entry+".11"      # give preference to the binary file
        elif os.path.isfile(entry):
            pass
        else:
            entry = os.path.join(linelistdir,entry)
            if os.path.isfile(entry+".11"): entry = entry+".11"

        assert(os.path.isfile(entry)), 'The line list '+entry+' is neither accessible in the working directory nor in the linelistdir ('+linelistdir+') folder'
    
        linelist[i] = entry
        if not os.path.isabs(entry): 
            ll = os.path.join (os.getcwd(), entry)
            if os.path.isfile(ll): linelist[i] = ll

        i = i + 1

    return(linelist)


def checkinput(wrange, vmicro, linelist):
    """
    Checking input parameters from user

    Parameters
    ----------
    wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
    vmicro: float, optional
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
    linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)

    Returns
    ------
    imode: int
      appropriate value for the variable imode, which specifies whether
      one will use many atomic lines (imode=0), just a few (imode=1),
      or none (H and HeII lines are an exception; imode=2)

    """

    # Determine imode
    #  imode = 0  is default, atoms and molecules, at least 2 line lists 
    #  synple sets IFMOL = 1 in 'tas' when an input molecular line list is used
    #  but does not set it when only an atomic line list is given
    #  imode = 2 for pure continuum
    #  imode = 1 for few-lines mode
    #  there are two more values which are possible but not considered in this routine
    #  imode = -3 for regular opacity tables (TLUSTY)
    #  imode = -4 for continuum-only (+ H and HeII lines) opacity tables

    assert (wrange[1] > wrange[0]),'the ending wavelength for the spectral range must be larger than the starting wavelength'
    if len(linelist) == 0: 
        imode = 2  # no atomic or molecular line list -> pure continuum and no molecules
    else:
        nlines, minlambda, maxlambda = getlinelistrange(linelist[0])

    # check
    if nlines > 10:
        assert (wrange[0] > minlambda and wrange[1] < maxlambda),'wrange exceeds the allow range ('+str(minlambda)+' to '+str(maxlambda)+')'
        imode = 0
    else:
        imode = 1

        assert (vmicro >= 0.0),'vmicro = '+str(vmicro)+' but cannot < 0.'
  
    return(imode)


def getlinelistrange(atomiclinelist):
    # Finds out min and max wavelengths for a line list

    if atomiclinelist[-3:] == '.11':
        file00 = atomiclinelist[:-3]+'.00'
        assert (os.path.isfile(file00)),'The file '+file00+' reporting linelist statistics is missing'
        f = open(file00,'r')
        lines = f.readlines()
        entries = lines[0].split()
        nlines = np.int64(entries[0])
        entries = lines[1].split()
        minlambda = float(entries[0])*10.
        entries = lines[2].split()
        maxlambda = float(entries[0])*10.
        f.close()
    else:
        f = open(atomiclinelist,'r')
        line = f.readline()
        entries = line.split()
        minlambda = float(entries[0])*10.
        fsize = os.path.getsize(atomiclinelist)
        f.seek(fsize-103)
        line = f.readline()
        f.close()
        entries = line.split()
        maxlambda = float(entries[0])*10.
        nlines = int(0.01 * fsize)

    return(nlines, minlambda,maxlambda)

def create_links(linelist):
    # Create soft links for line lists

    for i in range(len(linelist)):
        file = linelist[i]
        binaryfile = linelist[i][:-2]+'11'
        if os.path.isfile(binaryfile): file = binaryfile
        if i == 0: os.symlink(file,'fort.19')
        else: os.symlink(file,'fort.'+str(20-1+i))
    return()

def cleanup_fort():
    # Cleanup all fort* files

    files = os.listdir('.')
    for entry in files: 
        if os.path.islink(entry) and entry.startswith('fort'): os.unlink(entry)
        if os.path.isfile(entry) and entry.startswith('fort'): os.remove(entry)

    return()

def writetas(filename,nd,linelist,nonstd=None):
    #write non-std input parameters
    # input: filename -- str -- name of the non-std. param. file to print
    #        nd -- int -- number of layers in the model
    #        linelist -- list -- names of the linelist files (atomic first, then one 
    #				or more molecular ones
    #        nonstd -- dict -- additional entries to add to the tas file 
    #
  
    f = open(filename,'w')

    if nonstd is not None:
        for entry in nonstd.keys():
            f.write(entry + "=" + str(nonstd[entry])+"\n")

    f.write("ND= "+str(nd)+" \n")
    if len(linelist) > 1:  f.write("IFMOL= "+one+" \n")
    f.write("TMOLIM= 8000. \n")

    f.close()

    return()
