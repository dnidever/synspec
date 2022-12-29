import os
import gzip
import numpy as np
from astropy.table import Table

# The data directory
def datadir():
    """ Return the data/ directory."""
    fil = os.path.abspath(__file__)
    codedir = os.path.dirname(fil)
    datadir = codedir+'/data/'
    return datadir

# The test directory
def testdir():
    """ Return the test/ directory."""
    fil = os.path.abspath(__file__)
    codedir = os.path.dirname(fil)
    testdir = codedir+'/test/'
    return testdir

def readlines(fil=None,comment=None,raw=False,nreadline=None,noblank=False):
    """
    Read in all lines of a file.
    
    Parameters
    ----------
    file : str
         The name of the file to load.
    comment : str
         Comment line character to ignore (e.g., "#").
    raw : bool, optional, default is false
         Do not trim \n off the ends of the lines.
    nreadline : int, optional
         Read only this number of lines.  Default is to read all lines.
    noblank : boolean, optional
         Remove blank lines or lines with only whitespace.  Default is False.

    Returns
    -------
    lines : list
          The list of lines from the file

    Example
    -------

    .. code-block:: python

       lines = readlines("file.txt")

    """
    if fil is None: raise ValueError("File not input")
    # Read gzipped file
    if fil.endswith('gz'):
        fp = gzip.open(fil)
        contents = fp.read() # contents now has the uncompressed bytes of foo.gz
        fp.close()
        lines = contents.decode('utf-8') # u_str is now a unicode string
        lines = lines.split('\n')
    # Read normal ASCII file
    else:
        # Read all lines of normal ASCII file
        if nreadline is None:
            with open(fil,'r') as f:
                lines = f.readlines()
        # Only read a 
        else:
            with open(fil,'r') as f:
                lines = []
                for i in range(nreadline):
                    lines.append( f.readline() )
    # Remove blank lines
    if noblank:
        lines = [l for l in lines if l.strip()!='']
    # Strip newline off
    if raw is False: lines = [l.rstrip('\n') for l in lines]
    # Check for comment string:
    if comment is not None:
        lines = [l for l in lines if l[0]!=comment]                
    return lines

def writelines(filename=None,lines=None,overwrite=True,raw=False):
    """
    Write a list of lines to a file.
    
    Parameters
    ----------
    filename : str
        The filename to write the lines to.
    lines : list
         The list of lines to write to a file.
    overwrite : bool, optional, default is True
        If the output file already exists, then overwrite it.
    raw : bool, optional, default is False
        Do not modify the lines. Write out as is.

    Returns
    -------
    Nothing is returned.  The lines are written to `fil`.

    Example
    -------

    .. code-block:: python

       writelines("file.txt",lines)

    """
    # Not enough inputs
    if lines is None: raise ValueError("No lines input")
    if filename is None: raise ValueError("No file name input")
    # Check if the file exists already
    if os.path.exists(filename):
        if overwrite is True:
            os.remove(filename)
        else:
            print(filename+" already exists and overwrite=False")
            return
    # Modify the input as needed
    if raw is False:
        # List, make sure it ends with \n
        if type(lines) is list:
            for i,l in enumerate(lines):
                if l.endswith('\n') is False:
                    lines[i] += '\n'
            # Make sure final element does not end in \n
            #n = size(lines)
            #if n>1:
            #    if lines[-1].endswith('\n'):
            #        lines[-1] = lines[-1][0:-1]
            #else:
            #    if lines[0].endswith('\n'):
            #        lines = lines[0][0:-1]
    # Convert string to list
    if (type(lines) is str) | (type(lines) is np.str_): lines=list(lines)
    # Convert numpy array and numbers to list of strings
    if type(lines) is not list:
        if hasattr(lines,'__iter__'):
            lines = [str(l)+'\n' for l in lines]
            # Make sure final element does not end in \n        
            #if lines[-1].endswith('\n'): lines[-1] = lines[-1][0:-1]        
        else:
            lines = str(lines)
    # Write the file
    f = open(filename,'w')
    f.writelines(lines)
    f.close()
      
def read_synthfile(synthfile):
    """ Read the raw synthetic spectrum file."""

    lines = readlines(synthfile)
    nlines = len(lines)
    
    #ALL abundances NOT listed below differ from solar by   0.10 dex
    #element Li:  abundance =  0.28
    #MODEL:           Teff = 4150           log g = 2.5           vt= 2.00 M/H= 0.10 
    #   6695.000   6718.758      0.020      1.000
    # 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    # 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000    

    # Find which line MODEL is one
    for i in range(nlines):
        if lines[i].find('MODEL')>-1:
            break
    linestart = i
    
    # Get the fluxes
    fline = ' '.join(lines[linestart+2:])
    flux = np.array(fline.split()).astype(float)
    flux = 1-flux
    n = len(flux)

    # Get the wavelengths
    wline = lines[linestart+1]
    wstart,wend,wstep = wline.split()[0:3]
    wave = np.arange(n)*float(wstep)+float(wstart)

    return wave,flux

def read_sumfile(sumfile,solar=False,verbose=False):
    """
    This loads the MOOG "sumout" file 
    
    Parameters
    ----------
    sumfile    The MOOG "sumout" filename 
    /silent    No output to the screen. 
    /verbose   Output all of the MOOG "sumout" information. 
    
    Returns
    -------
    tab        IDL structure of the summary file.  For each line there 
                  are: WAVE, SPECIES, EP, LOGGF, EW, LOGRW, ABUND, 
                  DELAVG, ELEMENT. 
    elements  A structure for each unique element: ELEMENT, ATOMNUM, 
                    NLINES, ABUND, M_H, M_FE. 
    atmos      IDL structure with the atmosphere parameters (Teff,logg,microtrub,metal) 
                  used to run MOOG. 
    solar     The MOOG solar abundance structure. 

    Example
     -------

    tab,elements,atmos = moog_loadsum('sumout.txt')
    
    By D.Nidever  Sep 2009 
    """ 
     
    # Check the sumfile file 
    if os.path.exists(sumfile) == False:
        raise FileNotFoundError(sumfile)
    
    # Read the file
    lines = dln.readlines(sumfile)
    nlines = len(lines)
    
    # Loop through the lines 
    species = ''
    tab = []
    for i in range(nlines):      
        iline = lines[i] 
         
        # Beginning of new species 
        if iline[0:17] == 'Abundance Results': 
            species = iline[29:29+10]
            species = species.replace(' ','') 
        arr = iline.split()
        isnum = dln.isnumber(arr[9])
        #isnum = valid_num(arr[0]) 
         
        # Header line 
        if arr[0].strip() == 'Teff':
            eq1 = iline.find('=')
            teff = float(iline[eq1+1:eq1+7])
            eq2 = iline[eq1+1:].find('=')
            logg = float(iline[eq2+1:eq2+7])
            eq3 = iline[eq2+1:].find('=')
            microturb = float(iline[eq3+1:eq3+7])
            eq4 = iline[eq3+1:].find('=')
            metal = float(iline[eq3+1:eq3+7])
            #equalind1 = strpos(iline,'=',0) 
            #teff = float(strmid(iline,equalind1+1,6)) 
            #equalind2 = strpos(iline,'=',equalind1+1) 
            #logg = float(strmid(iline,equalind2+1,6)) 
            #equalind3 = strpos(iline,'=',equalind2+1) 
            #microturb = float(strmid(iline,equalind3+1,6)) 
            #equalind4 = strpos(iline,'=',equalind3+1) 
            #metal = float(strmid(iline,equalind4+1,6)) 
            atmos = {'teff':teff,'logg':logg,'microturb':microturb,'metal':metal} 
         
        if isnum == 1: 
            ncol = len(arr) 
            #fieldnames = ['WAVE','EP','LOGGF','EW','LOGRW','ABUND','DELAVG'] 
            #str = arr2str(reform(arr,ncol,1),fieldnames=fieldnames,/silent) 
            arr2 = np.array(arr).astype(float)
            tab1 = {'wave':arr2[0],'species':species,'ep':arr2[1],'loggf':arr2[2],'ew':arr2[3],
                    'logrw':arr2[4],'abund':arr2[5],'delavg':arr2[6]}
            tab.append(tab1)
 
    ntab = len(tab)
    tab = Table(tab)    
    if verbose:
        print('N lines = {:d}'.format(ntab))
     
    # Verbose output
    if verbose:
        print('MOOG OUTPUT')
        print('')
        for f in lines: print(f)
        print('')
     
    # Atmospheric parameters
    if verbose:
        print('Atmospheric Parameters:')
        print('  Teff={0:.1f} logg={1:.2f} vt={2:.2f} [M/H]={3:.3f}'.format(atmos['teff'],atmos['logg'],
                                                                            atmos['microturb'],atmos['metal']))
    
    # Get element names
    two = [sp[0:2] for sp in tab['species']]
    first = [sp[0:1] for sp in tab['species']]
    second = [sp[1:2] for sp in tab['species']]
    #two = strmid(tab['species'],0,2) 
    #first = strmid(tab['species'],0,1) 
    #second = strmid(tab['species'],1,1) 
    element = two 
    bd, = np.where(np.array(second)=='I')
    nbd = len(bd)
    if nbd > 0: 
        element[bd] = first[bd] 

    uelement = np.unique(element)
    nel = len(uelement)
     
    # Add element names
    #add_tag,str,'ELEMENT','',str     
    tab['element'] = '  '
    tab['element'] = element 

    if verbose:
        print('Nelements = {:f}'.format(nel))
     
    # MOOG abundance and element information 
    # 
    # MOOG gives the abundances in log eps(X) = log( N(X)/N(H) ) + 12.0 
    #  It also uses the Anders and Grevesse (1989) solar abundances 
    #  and a meteoritic Fe abundances of 7.52. 
    # 
    # metal deficiencies: [M/H] = log(N(X)/N(H)) - log(N(X)/N(H))_solar 
    #  or [X/H] = log eps(x) - log eps(x)_solar 
    #  e.g. [Fe/H] = log eps(Fe) - 7.52 
     
    # From MOOG src/Batom/f file 
    elnames = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
               'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
               'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
               'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
               'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
               'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
               'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
               'Lu','Hf','Ta','Wl','Re','Os','Ir','Pt','Au','Hg',
               'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
               'Pa','U ','Np','Pu','Am'] 
    #c 
    #c  xabu = the set of current solar (when available) or meteorite 
    #c  abundances, scaled to log(h) = 12.00 .  The data are from Anders 
    #c  and Grevesse (1989, Geochim.Cosmichim.Acta, v53, p197) and the solar 
    #c  values are adopted except for a) uncertain solar data, or b) Li, Be, 
    #c  and B, for which the meteoritic values are adopted. 
    #c  I was told to use the new Fe value of 7.52 as adopted in Sneden 
    #c  et al. 1992 AJ 102 2001. 
    #      data xabu/ 
    # 
    elabund = [12.00,10.99, 3.31, 1.42, 2.88, 8.56, 8.05, 8.93, 4.56, 8.09,
               6.33, 7.58, 6.47, 7.55, 5.45, 7.21, 5.5 , 6.56, 5.12, 6.36,
               3.10, 4.99, 4.00, 5.67, 5.39, 7.52, 4.92, 6.25, 4.21, 4.60,
               2.88, 3.41, 2.37, 3.35, 2.63, 3.23, 2.60, 2.90, 2.24, 2.60,
               1.42, 1.92, 0.00, 1.84, 1.12, 1.69, 1.24, 1.86, 0.82, 2.0,
               1.04, 2.24, 1.51, 2.23, 1.12, 2.13, 1.22, 1.55, 0.71, 1.50,
               0.00, 1.00, 0.51, 1.12, 0.33, 1.1 , 0.50, 0.93, 0.13, 1.08,
               0.12, 0.88, 0.13, 0.68, 0.27, 1.45, 1.35, 1.8 , 0.83, 1.09,
               0.82, 1.85, 0.71, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.12,
               0.00, 0.00, 0.00, 0.00, 0.00] 
     
    #abund_solar = replicate({atomnum:0L,element:'',abund:0.0},len(elnames))
    dts = [('atomnum',int),('element',str,2),('abund',float)]
    abund_solar = np.zeros(len(elnames),dtype=np.dtype(dts))
    abund_solar['atomnum'] = np.arange(len(elnames))+1
    abund_solar['element'] = elnames 
    abund_solar['abund'] = elabund 
     
     
    # Add ATOMNUM to structure
    tab['atomnum'] = -1
    for i in range(ntab): 
        ind, = np.where(abund_solar['element']==tab['element'][i])
        nind = len(ind)
        if nind > 0: 
            tab['atomnum'][i] = abund_solar['atomnum'][ind[0]]
 
    # Getting Fe abundance 
    ind, = np.where(tab['element'] == 'Fe')
    nind = len(ind)
    if nind > 0: 
        abundfe = np.median(tab['abund'][ind]) 
    else: 
        print('NO Fe abundance - using solar')
        abundfe = 7.52 
    Fe_H = abundfe - 7.52  # [Fe/H] 
     
    # Start the "element" structure
    dt = [('atomnum',int),('element',str,2),('nlines',int),('abund',float),('M_H',float),('M_Fe',float)]
    elemtab = np.zeros(nel,dtype=np.dtype(dt))
    elemtab = Table(elemtab)
    
    # Loop through the elements 
    for i in range(nel): 
        ind, = np.where(tab['element'] == uelement[i])
        nind = len(ind)            
        tab1 = tab[ind]
        if nind > 1: 
            medabund = np.median(tab1['abund']) 
        else: 
            medabund = tab1['abund'][0]
         
        # Getting solar abundance 
        indsol, = np.where(abund_solar['element'] == uelement[i])
        atomnum = abund_solar['atomnum'][indsol[0]]
        abund_sol = abund_solar['abund'][indsol[0]]
         
        # Metal deficiency 
        #  metal deficiencies: [M/H] = log(N(X)/N(H)) - log(N(X)/N(H))_solar 
        #    or [X/H] = log eps(X) - log eps(X)_solar 
        #    e.g. [Fe/H] = log eps(Fe) - 7.52 
        M_H = medabund - abund_sol 
         
        # Metal deficiency versus Fe 
        # [X/H] = [X/Fe] + [Fe/H] and therefore 
        # [X/Fe] = [X/H] - [Fe/H] 
        M_Fe = M_H - Fe_H 
         
        # Put it in the structure 
        elemtab['atomnum'][i] = atomnum 
        elemtab['element'][i] = uelement[i] 
        elemtab['nlines'][i] = nind 
        elemtab['abund'][i] = medabund 
        elemtab['M_H'][i] = M_H 
        elemtab['M_Fe'][i] = M_Fe 
         
        #form='(A-4,I-4,I-4,3F8.2)' 
        #print,format=form,uelement[i],atomnum,nind,medabund,M_H,M_Fe 
        #form='(A-4,A-14,A-15)' 
        #print,format=form,uelement[i],'Abund = '+strtrim(string(medabund,format='(F9.2)'),2),  ;     'Nlines = '+strtrim(nind,2) 
 
     
    # Sort the element array by atomnum
    elemtab.sort('atomnum')
    #si = np.argsort(elemtab.atomnum) 
    #elemtab = elemtab[si] 
     
     
    # Print out the Element information
    if verbose:
        print('-------------------------------------')
        print('Elem A.Num Nline Abund  [M/H] [M/Fe]')
        print('-------------------------------------') 
        for i in range(nel): 
            etab = elemtab[i] 
            formt = '{0:3fs}{1:5d}{2:5d}{3:8.2f}{4:7.2f}{5:7.2f}'
            print(form.format(etab['element'],etab['atomnunm'],etab['nlines'],etab['abund'],etab['M_H'],etab['M_Fe']))
        print('-------------------------------------' )
     
     
    # Delete the MOOG sumout file 
    #if keyword_set(delete) then begin 
    #  if not keyword_set(silent) then print,'Deleting MOOG "sumout" file' 
    #  FILE_DELETE,sumfile,/allow 
    #endif 
 
