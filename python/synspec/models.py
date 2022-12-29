import os
import numpy as np

def read_model(modelfile,verbose=False):
    """
    Reads a model atmosphere into a structure
  
    Parameters
    ----------  
    modelfile : str
      file with a model atmosphere
      
    Returns
    -------
    atmostype :  str
      type of model atmosphere (kurucz/marcs/phoenix/tlusty)
    teff : float
      effective temperature (K)
    logg : float
      log10 of the surface gravity (cm s-2)
    vmicro : float
      microturbulence velocity (km/s)
    abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
    nd: int
      number of depths (layers) of the model
    atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density
    """

    # Check
    if not os.path.isfile(modelfile):
        mf = os.path.join(modeldir,modelfile)
        if os.path.isfile(mf): modelfile = mf

    atmostype = identify_atmostype(modelfile,verbose=verbose)

    if atmostype == 'kurucz':
        teff, logg, vmicro, abu, nd, atmos = read_kurucz_model(modelfile) 
    if atmostype == 'marcs':
        teff, logg, vmicro, abu, nd, atmos = read_marcs_model2(modelfile)
    if atmostype == 'phoenix':
        teff, logg, vmicro, abu, nd, atmos = read_phoenix_model(modelfile)
    if atmostype == 'tlusty':
        teff, logg, vmicro, abu, nd, atmos = read_tlusty_model(modelfile)

    return (atmostype,teff,logg,vmicro,abu,nd,atmos)


def identify_atmostype(modelfile,verbose=False):
    """
    Idenfies the type of model atmosphere in an input file

    Valid options are kurucz, marcs, tlusty (.7) or phoenix

    Parameters
    ----------
    modelfile: str
      file with a model atmosphere
    
    Returns
    -------
    atmostype: str
      can take the value 'kurucz', 'marcs', 'tlusty' or 'phoenix' 

    """

    if ('PHOENIX' in modelfile and 'fits' in modelfile):
        atmostype = 'phoenix'
    else: 
        if modelfile[-3:] == '.gz':
            f = gzip.open(modelfile,'rt')
        else:
            f = open(modelfile,'r')
        line = f.readline()
        if verbose:
            print('modelfile / line=',modelfile,line)
        #type(line)
        if ('TEFF' in line):
            atmostype = 'kurucz'
        else: 
            line = f.readline()
            if ('Teff' in line):
                atmostype = 'marcs'
            else:
                atmostype = 'tlusty'
        f.close()
   
    return(atmostype)


def read_kurucz_model(modelfile):
    """
    Reads a Kurucz model atmospheres
  
    Parameters
    ----------
    modelfile: str
      file name  
  
    Returns
    -------
    teff : float
      effective temperature (K)
    logg : float
      log10 of the surface gravity (cm s-2)
    vmicro : float
      microturbulence velocity (km/s)
    abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
    nd: int
      number of depths (layers) of the model
    atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density  
  
    """

    f = open(modelfile,'r')
    line = f.readline()
    entries = line.split()
    assert (entries[0] == 'TEFF' and entries[2] == 'GRAVITY'), 'Cannot find Teff and logg in the file header'
    teff = float(entries[1])
    logg = float(entries[3])

    while entries[0] != 'ABUNDANCE':  
        line = f.readline()
        entries = line.split()

    abu = []

    if entries[1] == 'SCALE': 
        scale = float(entries[2])
    

    while entries[0] == 'ABUNDANCE':
        i = 0
        for word in entries: 
            if (word == 'CHANGE'): w = i
            i = i + 1 
        for i in range(int((len(entries)-w-1)/2)):
            z = int(entries[w+1+2*i])
            if (z == 1): nhntot = float(entries[w+2+2*i])
            if (z < 3): abu.append(float(entries[w+2+2*i]) / nhntot) 
            else: abu.append(scale*10.**(float(entries[w+2+2*i])) / nhntot)

        line = f.readline()
        entries = line.split() 

    assert (entries[0] == 'READ'), 'I cannot find the header of the atmospheric table in the input Kurucz model'

    nd = int(entries[2]) - 1
    line = f.readline()
    entries = line.split()
    line = f.readline()
    entries = line.split()
    vmicro = float(entries[6])/1e5

    dm = [ float(entries[0]) ]
    t = [ float(entries[1]) ]
    p = [ float(entries[2]) ]
    ne = [ float(entries[3]) ] 

    for i in range(nd-1):
        line = f.readline()
        entries = line.split()
        dm.append( float(entries[0]))
        t.append(  float(entries[1]))
        p.append(  float(entries[2]))
        ne.append( float(entries[3]))

    atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                                'formats':('f', 'f', 'f','f')}) 
    atmos['dm'] = dm
    atmos['t'] = t
    atmos['p'] = p
    atmos['ne'] = ne

    return (teff,logg,vmicro,abu,nd,atmos)


def read_marcs_model(modelfile):
    """
    Reads a MARCS model atmospheres
  
    Parameters
    ----------
    modelfile: str
      file name. It can be a gzipped (.gz) file
  
    Returns
    -------
    teff : float
      effective temperature (K)
    logg : float
      log10 of the surface gravity (cm s-2)
    vmicro : float
      microturbulence velocity (km/s)
    abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
    nd: int
      number of depths (layers) of the model
    atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density  
  
    """  

    if modelfile[-3:] == '.gz':
        f = gzip.open(modelfile,'rt')
    else:
        f = open(modelfile,'r')
    line = f.readline()
    line = f.readline()
    entries = line.split()
    assert (entries[1] == 'Teff'), 'Cannot find Teff in the file header'
    teff = float(entries[0])
    line = f.readline()
    line = f.readline()
    entries = line.split()
    assert (entries[1] == 'Surface' and entries[2] == 'gravity'), 'Cannot find logg in the file header'
    logg = np.log10(float(entries[0]))
    line = f.readline()
    entries = line.split()
    assert (entries[1] == 'Microturbulence'), 'Cannot find vmicro in the file header'
    vmicro = float(entries[0])

    while entries[0] != 'Logarithmic':  
        line = f.readline()
        entries = line.split()

    abu = []
    line = f.readline()
    entries = line.split()

    i = 0
    while entries[1] != 'Number':
        for word in entries: 
            abu.append( 10.**(float(word)-12.0) )
            i = i + 1 
        line = f.readline()
        entries = line.split() 

    if i < 99: 
        for j in range(99-i):
            abu.append(1e-111)
            i = i + 1

    nd = int(entries[0])
    line = f.readline()
    entries = line.split()

    assert (entries[0] == 'Model'), 'I cannot find the header of the atmospheric table in the input MARCS model'
    
    line = f.readline()
    line = f.readline()
    entries = line.split()

    t = [ float(entries[4]) ]
    p = [ float(entries[6]) ]
    ne = [ float(entries[5]) / bolk / float(entries[4]) ] 

    for i in range(nd-1):
        line = f.readline()
        entries = line.split()

        t.append(  float(entries[4]))
        p.append(  float(entries[6]))
        ne.append( float(entries[5]) / bolk / float(entries[4]))

    line = f.readline()
    line = f.readline()
    entries = line.split()

    dm = [ float(entries[-1]) ]

    for i in range(nd-1):
        line = f.readline()
        entries = line.split()

        dm.append(  float(entries[-1]))

    atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                                'formats':('f', 'f', 'f','f')}) 
    atmos['dm'] = dm
    atmos['t'] = t
    atmos['p'] = p
    atmos['ne'] = ne

    return (teff,logg,vmicro,abu,nd,atmos)


def read_marcs_model2(modelfile): 
    """
    Reads a MARCS model atmospheres. 
    While read_marcs_model returns T, Pg and Ne in the structure 'atmos'
    read_marcs_model2 returns T, rho, mmw, and Ne.
  
    Parameters
    ----------
    modelfile: str
      file name. It can be a gzipped (.gz) file
  
    Returns
    -------
    teff : float
      effective temperature (K)
    logg : float
      log10 of the surface gravity (cm s-2)
    vmicro : float
      microturbulence velocity (km/s)
    abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
    nd: int
      number of depths (layers) of the model
    atmos: numpy structured array
      array with the run with depth of column mass, temperature, density, 
      mean molecular weight and electron number density  
  
    """  

    if modelfile[-3:] == '.gz':
        f = gzip.open(modelfile,'rt')
    else:
        f = open(modelfile,'r')
    line = f.readline()
    line = f.readline()
    entries = line.split()
    assert (entries[1] == 'Teff'), 'Cannot find Teff in the file header'
    teff = float(entries[0])
    line = f.readline()
    line = f.readline()
    entries = line.split()
    assert (entries[1] == 'Surface' and entries[2] == 'gravity'), 'Cannot find logg in the file header'
    logg = np.log10(float(entries[0]))
    line = f.readline()
    entries = line.split()
    assert (entries[1] == 'Microturbulence'), 'Cannot find vmicro in the file header'
    vmicro = float(entries[0])

    while entries[0] != 'Logarithmic':  
        line = f.readline()
        entries = line.split()

    abu = []
    line = f.readline()
    entries = line.split()

    i = 0
    while entries[1] != 'Number':
        for word in entries: 
            abu.append( 10.**(float(word)-12.0) )
            i = i + 1 
        line = f.readline()
        entries = line.split() 

    if i < 99: 
        for j in range(99-i):
            abu.append(1e-111)
            i = i + 1

    nd = int(entries[0])
    line = f.readline()
    entries = line.split()

    assert (entries[0] == 'Model'), 'I cannot find the header of the atmospheric table in the input MARCS model'

    line = f.readline()
    line = f.readline()
    entries = line.split()

    t = [ float(entries[4]) ]
    p = [ float(entries[6]) ]
    ne = [ float(entries[5]) / bolk / float(entries[4]) ] 

    for i in range(nd-1):
        line = f.readline()
        entries = line.split()

        t.append(  float(entries[4]))
        p.append(  float(entries[6]))
        ne.append( float(entries[5]) / bolk / float(entries[4]))

    line = f.readline()
    line = f.readline()
    entries = line.split()

    rho = [ float(entries[3]) ]
    dm = [ float(entries[-1]) ]
    mmw = [ float(entries[4]) ]

    for i in range(nd-1):
        line = f.readline()
        entries = line.split()

        rho.append( float(entries[3]))
        dm.append(  float(entries[-1]))
        mmw.append(  float(entries[4]))

    atmos = np.zeros(nd, dtype={'names':('dm', 't', 'rho','mmw','ne'),
                                'formats':('f', 'f', 'f','f','f')}) 
    atmos['dm'] = dm
    atmos['t'] = t
    atmos['rho'] = rho
    atmos['mmw'] = mmw
    atmos['ne'] = ne

    return (teff,logg,vmicro,abu,nd,atmos)


def read_tlusty_model(modelfile,startdir=None): 
    """
    Reads a Tlusty model atmosphere. 

    Parameters
    ----------
    modelfile: str
      file name (.7, .8, or .22). It will look for the complementary .5 file to read
      the abundances and the micro (when specified in the non-std. parameter file)

    startdir: str
      directory where the calculations are initiated. The code will look at that
      location to find the tlusty model atom directory and the non-std. parameter
      file when a relative path is provided
      (default is None, indicating it is the current working directory)
  
    Returns
    -------
    teff : float
      effective temperature (K)
    logg : float
      log10 of the surface gravity (cm s-2)
    vmicro : float
      microturbulence velocity (km/s), by default 0.0 unless set with the parameter
      VTB in the non-std. parameter file specified in the .5 file
    abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
    nd: int
      number of depths (layers) of the model
    atmos: numpy structured array
      array with the run with depth of column mass, temperature, density
      (other variables that may be included, e.g. populations for NLTE models, 
      are ignored). 

    """  

    assert ((modelfile[-2:] == ".8") | (modelfile[-2:] == ".7") | (modelfile[-3:] == ".22")), 'Tlusty models should end in .7, .8, or .22'
    if modelfile[-2] == ".":
        madaffile = modelfile[:-1]+"5"
    else:
        madaffile = modelfile[:-2]+"5"    
    assert (os.path.isfile(madaffile)),'Tlusty model atmosphere file '+modelfile+' should come with an associated .5 file'

    if startdir is None: startdir = os.getcwd()

    #we start reading the .5
    f = open(madaffile,'r')
    line = f.readline()
    entries = line.split()
    teff = float(entries[0])
    logg = float(entries[1])
    line = f.readline()
    line = f.readline()
    entries = line.split()
    nonstdfile = entries[0][1:-1]

    nonstdfile0 = nonstdfile
    if nonstdfile != '':
        if not os.path.isabs(nonstdfile): 
            mf = os.path.join(startdir,nonstdfile)
            if os.path.isfile(mf): 
                nonstdfile = mf
            else:
                mf = os.path.join(modeldir,nonstdfile)
                nonstdfile = mf

        assert (os.path.exists(nonstdfile)), 'The non-std parameter file indicated in the tlusty model, '+nonstdfile0+', is not present' 

    nonstd = {}
    if nonstdfile != '':
        assert (os.path.isfile(nonstdfile)),'Tlusty model atmosphere file '+modelfile+' invokes non-std parameter file, '+nonstdfile+' which is not present'

        ns = open(nonstdfile,'r')
        nonstdarr = ns.readlines()
        ns.close()
        for entry in nonstdarr:
            entries = entry.replace('\n','').split(',')
            for piece in entries:
                sides = piece.split('=')
                nonstd[sides[0].replace(' ','')]= sides[1].replace(' ','')

        print('Tlusty nonstd params=',nonstd)

    # the micro might be encoded as VTB in the nonstdfile!!
    # this is a temporary patch, but need to parse that file
    vmicro = 0.0
    if 'VTB' in nonstd: vmicro = float(nonstd['VTB'])

    line = f.readline()
    line = f.readline()
    entries = line.split()
    natoms = int(entries[0])
  
    abu = []
    for i in range(natoms):
        line = f.readline()
        entries = line.split()
        abu.append( float(entries[1]) )

    if i < 98: 
        for j in range(98-i):
            abu.append(1e-111)
            i = i + 1

    f.close()

    # now the .8
    f = open(modelfile,'r')
    line = f.readline()
    entries = line.split()
    nd = int(entries[0])
    numpar = int(entries[1])
    if (numpar < 0): 
        numpop = abs(numpar) - 4 
    else:
        numpop = numpar - 3

    assert (len(entries) == 2), 'There are more than two numbers in the first line of the model atmosphere'

    dm = read_multiline_fltarray(f,nd)
    atm = read_multiline_fltarray(f,nd*abs(numpar))
    f.close()

    atm = np.reshape(atm, (nd,abs(numpar)) )

    if (numpar < 0):  # 4th column is number density n
        if (numpop > 0): # explicit (usually NLTE) populations
            if modelfile[-2] == ".":  # NLTE populations or departure coefficients
                tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('n','f'), ('pop', 'f', (numpop))])
            else: 
                tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('n','f'), ('dep', 'f', (numpop))])
        else:
            tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('n','f')])  
    else:
        if (numpop > 0):
            if modelfile[-2] == ".": # NLTE populations or departure coefficients
                tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('pop', 'f', (numpop))])
            else:
                tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('dep', 'f', (numpop))])
        else:
            tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f') ])

    atmos = np.zeros(nd, dtype=tp)

    atmos['dm'] = dm
    atmos['t'] = atm [:,0]
    atmos['ne'] = atm [:,1]
    atmos['rho'] = atm [:,2]
    if (numpar < 0): atmos['n'] = atm [:,3]
    if (numpop > 0): 
        if modelfile[-2] == ".":
            atmos['pop'] = atm [:,4:]
        else:
            atmos['dep'] = atm [:,4:]

    return (teff,logg,vmicro,abu,nd,atmos)


def read_tlusty_extras(modelfile,startdir=None):
    """
    Identifies and reads the non-std parameter file and its content, finds out the 
    number of parameters in the model, whether the file contains populations or departure
    coefficients, and the name of the data directory for Tlusty 
    model atmospheres. 
    
    Parameters
    ----------
    modelfile: str
      file name (.8, .7 or .22). It will look for the complementary .5 file to read
      the abundances and other information

    startdir: str
      directory where the calculations are initiated. The code will look at that
      location to find the tlusty model atom directory and the non-std. parameter
      file when a relative path is provided
      (default is None, indicating it is the current working directory)
  
    Returns
    -------
    madaffile: str
       model atom data and abundance file (.5 Tlusty file)
    nonstdfile: str
       non-std parameter file 
    nonstd: dict
       content of the non-std parameter file
    numpar: int
       number of parameters (can be negative when the model includes number density)
    datadir: str
       name of the model atom directory
    inlte: int
       0 when the populations are to be computed internally by synspec (LTE)
       1 the Tlusty model contains populations
      -1 the Tlusty model contains departure coefficients
    atommode: list
       mode for each of the atoms included. The code indicates
       0= not considered
       1= implicit (no cont. opacity)
       2= explicit  (see synspec man.)
       4= semi-explicit (see synspec man.)
       5= quasi-explicit  (see synspec. man)
    atominfo: list
       all the lines in the file that provide info on the model atoms used
  
    """  

    assert ((modelfile[-2:] == ".8") | (modelfile[-2:] == ".7") | (modelfile[-3:] == ".22")), 'Tlusty models should end in .7, .8, or .22'
    if modelfile[-2] == ".":
        madaffile = modelfile[:-1]+"5"
    else:
        madaffile = modelfile[:-2]+"5"    
    assert (os.path.isfile(madaffile)),'Tlusty model atmosphere file '+modelfile+' should come with an associated .5 file'

    if startdir is None: startdir = os.getcwd()

    # we start reading the .5
    f = open(madaffile,'r')
    line = f.readline()
    line = f.readline()
    line = f.readline()
    entries = line.split()
    nonstdfile = entries[0][1:-1]

    nonstdfile0 = nonstdfile  
    if nonstdfile != '':
        if not os.path.isabs(nonstdfile): 
            mf = os.path.join(startdir,nonstdfile)
            if os.path.isfile(mf): 
                nonstdfile = mf
            else:
                mf = os.path.join(modeldir,nonstdfile)
                nonstdfile = mf

        assert (os.path.exists(nonstdfile)), 'The non-std parameter file indicated in the tlusty model, '+nonstdfile0+', is not present' 

    nonstd = {}
    if nonstdfile != '':
        assert (os.path.isfile(nonstdfile)),'Tlusty model atmosphere file '+modelfile+' invokes non-std parameter file, '+nonstdfile+' which is not present'

        ns = open(nonstdfile,'r')
        nonstdarr = ns.readlines()
        ns.close()
        for entry in nonstdarr:
            entries = entry.replace('\n','').split(',')
            for piece in entries:
                sides = piece.split('=')
                nonstd[sides[0].replace(' ','')]= sides[1].replace(' ','')


    line = f.readline()
    line = f.readline()
    entries = line.split()
    natoms = int(entries[0])
  
    atommode = []
    for i in range(natoms):
        line = f.readline()
        entries = line.split()
        atommode.append(int(entries[0]))

    atominfo = []
    #keep reading until you find 'dat' to identify data directory 
    line = f.readline()
    while True: 
        atominfo.append(line)
        if '.dat' in line: break
        line = f.readline()

    entries = line.split()
    cadena = entries[-1][1:-1]
    datadir, file = os.path.split(cadena)


    datadir0 = datadir
    if datadir != '':
        if not os.path.isabs(datadir): 
            mf = os.path.join(startdir,datadir)
            if os.path.exists(mf): 
                datadir = mf
            else:
                mf = os.path.join(synpledir,datadir)
                datadir = mf

    # continue reading the rest of the file into atominfo
    line = f.readline()
    while True:
        if line == '': break
        atominfo.append(line)
        line = f.readline()

        assert (os.path.exists(datadir)), 'The datadir indicated in the tlusty model, '+datadir0+', is not present' 
        
    f.close()

    # now the .8
    f = open(modelfile,'r')
    line = f.readline()
    entries = line.split()
    nd = int(entries[0])
    numpar = int(entries[1])
    if abs(numpar) > 4: 
        inlte = 1 
    else: 
        inlte = 0

    if (modelfile[-3:] == ".22"): inlte = -1

    f.close()

    return (madaffile, nonstdfile, nonstd, numpar, datadir, inlte, atommode, atominfo)


def read_phoenix_model(modelfile):
    """
    Reads a FITS Phoenix model atmospheres
  
    Parameters
    ----------
    modelfile: str
      file name  
  
    Returns
    -------
    teff : float
      effective temperature (K)
    logg : float
      log10 of the surface gravity (cm s-2)
    vmicro : float
      microturbulence velocity (km/s)
    abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
    nd: int
      number of depths (layers) of the model
    atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density  
  
    """  

    from astropy.io import fits

    h = fits.open(modelfile)[0].header
    f = fits.open(modelfile)[1].data

    nd = len(f['temp'])

    teff = float(h['PHXTEFF'])
    logg = float(h['PHXLOGG'])
    vmicro = float(h['PHXXI_L'])

    m_h = float(h['PHXM_H'])
    alpha = float(h['PHXALPHA'])
  
    symbol, mass,sol = elements(reference='husser') 
    abu = sol 
    z_metals = np.arange(97,dtype=int) + 3
    z_alphas = np.array([8,10,12,14,16,20,22],dtype=int)
    for i in range(len(z_metals)): abu[z_metals[i] - 1] = abu[z_metals[i] - 1] + m_h
    for i in range(len(z_alphas)): abu[z_alphas[i] - 1] = abu[z_alphas[i] - 1] + alpha

    atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                                'formats':('f', 'f', 'f','f')}) 

    atmos['dm'] = f['pgas'] / 10.**logg
    atmos['t'] = f['temp']
    atmos['p'] = f['pgas']
    atmos['ne'] = f['pe']/ bolk / f['temp']

    return (teff,logg,vmicro,abu,nd,atmos)


def read_phoenix_text_model(modelfile):
    """
    Reads a plain-text Phoenix model atmospheres
  
    Parameters
    ----------
    modelfile: str
      file name  
  
    Returns
    -------
    teff : float
      effective temperature (K)
    logg : float
      log10 of the surface gravity (cm s-2)
    vmicro : float
      microturbulence velocity (km/s)
    abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
    nd: int
      number of depths (layers) of the model
    atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density  
  
    """  

    f = open(modelfile,'r')
    line = f.readline()
    while line[0:4] != " no.":
        line = f.readline()
    entries = line.split()
    nd = int(entries[5])
    print('nd=',nd)
    while line[0:14] != " model:   teff":
        line = f.readline()
    entries = line.split()
    teff = float(entries[3])
    print('teff=',teff)
    line = f.readline()
    line = f.readline()
    entries = line.split()
    assert (entries[0] == 'log(g):' and entries[2] == '[cm/s**2]'), 'Cannot find logg in the file header'
    logg = float(entries[1])
    print('logg=',logg)
    line = f.readline()
    while line[0:22] !=  "  Element abundances :":  
        line = f.readline()

    symbol,mass,sol = elements()
    
    sy = []
    ab = []

    while line[0:29] !=  "  Element abundances relative":  
        line = f.readline()
        #print(line)
        if line[0:9] == ' element:':
            entries = line.split()
            for word in entries[1:]: sy.append(word)
        if line[0:11] == ' abundance:':
            entries = line.split()
            for word in entries[1:]: ab.append(word)

    assert (len(sy) == len(ab)), 'different elements in arrays sy (elemental symbols) and ab (abundances)'

    abu = np.ones(99)*1e-99
    i = 0
    for item in sy:
        try:
            index = symbol.index(item)
            abu[index] =  10.**(float(ab[i])-12.) 
        except ValueError:
            print("the symbol ",item," is not recognized as a valid element")
        i = i + 1

    print('abu=',abu)

    while line[0:72] !=  "   l        tstd temperature        pgas          pe     density      mu":  
        line = f.readline()

    line = f.readline()
    entries = line.split()

    t = [ float(entries[2].replace('D','E')) ]
    p = [ float(entries[3].replace('D','E')) ]
    ne = [ float(entries[4].replace('D','E')) / bolk / float(entries[2].replace('D','E')) ] 
    dm = [ float(entries[3].replace('D','E')) / 10.**logg ] #assuming hydrostatic equil. and negliglible radiation and turb. pressure

    for i in range(nd-1):
        line = f.readline()
        entries = line.split()

        t.append(  float(entries[2].replace('D','E')))
        p.append(  float(entries[3].replace('D','E')))
        ne.append( float(entries[4].replace('D','E')) / bolk / float(entries[2]))
        dm.append ( float(entries[3].replace('D','E')) / 10.**logg )

    vmicro = 0.0
    while (line[0:6] != " greli"):
        line = f.readline()
        if line == '':
            print('Cannot find a value for vmicro (vturb) in the model atmosphere file ',modelfile)
            break
  
    if line != '':
        entries = line.split()
        vmicro = float(entries[5])

    atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                                'formats':('f', 'f', 'f','f')}) 
    atmos['dm'] = dm
    atmos['t'] = t
    atmos['p'] = p
    atmos['ne'] = ne

    return (teff,logg,vmicro,abu,nd,atmos)

