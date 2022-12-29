import numpy as np

def periodic(n):
    """ Routine to get element name / atomic number conversion """
    elem = np.array(['H','He','Li','Be','B','C','N','O','F','Ne',
                     'Na','Mg','Al','Si','P','S','Cl','Ar',
                     'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
                     'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
                     'Cs','Ba','La','Ce','Pr','Nd'])
    if isinstance(n,str):
        j = np.where(elem == n)[0]
        return j+1
    else:
        if n == 0:
            return ''    
        else:
            return elem[n-1]

def solar(el=None):
    """ Return solar abundances """
    sunabund_2007  =  np.array([
       12.00, 10.93,  1.05,  1.38,  2.70,  8.39,  7.78,  8.66,  4.56,   #  1 -  9
        7.84,  6.17,  7.53,  6.37,  7.51,  5.36,  7.14,  5.50,  6.18,   # 10 - 18
        5.08,  6.31,  3.17,  4.90,  4.00,  5.64,  5.39,  7.45,  4.92,   # 19 - 27
        6.23,  4.21,  4.60,  2.88,  3.58,  2.29,  3.33,  2.56,  3.25,   # 28 - 36
        2.60,  2.92,  2.21,  2.58,  1.42,  1.92, -99.0,  1.84,  1.12,   # 37 - 45
        1.66,  0.94,  1.77,  1.60,  2.00,  1.00,  2.19,  1.51,  2.24,   # 46 - 54
        1.07,  2.17,  1.13,  1.70,  0.58,  1.45, -99.0,  1.00,  0.52,   # 55 - 63
        1.11,  0.28,  1.14,  0.51,  0.93,  0.00,  1.08,  0.06,  0.88,   # 64 - 72
       -0.17,  1.11,  0.23,  1.25,  1.38,  1.64,  1.01,  1.13,  0.90,   # 73 - 81
        2.00,  0.65, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.06,   # 82 - 90
       -99.0, -0.52, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0 ]) # 91 - 99
    if el is None : return sunabund_2007
    n = periodic(el)
    return sunabund_2007[n-1]

def rydberg(n1,n2):
    """ Rydberg formula to give H wavelengths for transitions between levels n1 and n2. """
    r = 1.0973731568e7
    me = 9.109382e-31
    mprot = 1.672621e-27
    rm = r/(1+me/mprot)

    l = rm*(1./n1**2-1./n2**2)
    w = 1./l*1.e10
    return w

def hlines(plot=None,yloc=0.,n1=4,n2=range(11,22)):
    """ Return approximate (Rydberg) location of H lines, defaulting to lines in APOGEE spectra. """
    h = []
    for n in n2:
        h.append(rydberg(n1,n))
    h = np.array(h)
    if plot is not None :
        plots.plotp(plot,h,h*0.+yloc)
    return h

def elements(reference=None):
    """
    Reads the solar elemental abundances
  
    Parameters
    ----------
    reference: string, optional
        set to 'husser; for the abundances adopted for Phoenix models by Huser et al. (2013),
        set to 'basti' for abundances adopted for the BaSTI stellar models (Hidalgo et al. 2018),
        set to 'ags2005' for Asplund et al. (2005) are used -- consistent with
        the MARCS (Gustafsson et al. 2008) models and and Kurucz (Meszaros et al. 2012)
        Kurucz model atmospheres.
        (default, reference=None, will trigger 'ags2005')
        
    Returns
    -------
    symbol: numpy array of str
        element symbols
    mass: numpy array of floats
        atomic masses (elements Z=1-99)
    sol: numpy array of floats
        solar abundances N/N(H)
  
    """

    symbol = [
        'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne', 
        'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca', 
        'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn', 
        'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ,'Zr', 
        'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', 
        'Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd', 
        'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', 
        'Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg', 
        'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', 
        'Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es' ]
    
    mass = [ 1.00794, 4.00260, 6.941, 9.01218, 10.811, 12.0107, 14.00674, 15.9994,
             18.99840, 20.1797, 22.98977, 24.3050, 26.98154, 28.0855, 30.97376, 
             32.066, 35.4527, 39.948, 39.0983, 40.078, 44.95591, 47.867, 50.9415, 
             51.9961, 54.93805, 55.845, 58.93320, 58.6934, 63.546, 65.39, 69.723, 
             72.61, 74.92160, 78.96, 79.904, 83.80, 85.4678, 87.62, 88.90585, 
             91.224, 92.90638, 95.94, 98., 101.07, 102.90550, 106.42, 107.8682, 
             112.411, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.29, 
             132.90545, 137.327, 138.9055, 140.116, 140.90765, 144.24, 145, 150.36, 
             151.964, 157.25, 158.92534, 162.50, 164.93032, 167.26, 168.93421, 
             173.04, 174.967, 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 
             195.078, 196.96655, 200.59, 204.3833, 207.2, 208.98038, 209., 210., 
             222., 223., 226., 227., 232.0381, 231.03588, 238.0289, 237., 244., 
             243., 247., 247., 251., 252. ]
    
    if reference == 'husser':
        #a combination of meteoritic/photospheric abundances from Asplund et al. 2009
        #chosen for the Husser et al. (2013) Phoenix model atmospheres
        sol = [  12.00, 10.93,  3.26,  1.38,  2.79,  8.43,  7.83,  8.69,  4.56,  7.93, 
                 6.24,  7.60,  6.45,  7.51,  5.41,  7.12,  5.50,  6.40,  5.08,  6.34, 
                 3.15,  4.95,  3.93,  5.64,  5.43,  7.50,  4.99,  6.22,  4.19,  4.56, 
                 3.04,  3.65,  2.30,  3.34,  2.54,  3.25,  2.36,  2.87,  2.21,  2.58, 
                 1.46,  1.88, -9.99,  1.75,  1.06,  1.65,  1.20,  1.71,  0.76,  2.04, 
                 1.01,  2.18,  1.55,  2.24,  1.08,  2.18,  1.10,  1.58,  0.72,  1.42, 
                 -9.99,  0.96,  0.52,  1.07,  0.30,  1.10,  0.48,  0.92,  0.10,  0.92, 
                 0.10,  0.85, -0.12,  0.65,  0.26,  1.40,  1.38,  1.62,  0.80,  1.17,
                 0.77,  2.04,  0.65, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99,  0.06,   
                 -9.99, -0.54, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99 ]

    elif reference == 'basti':
        #recommended solar abundances from Lodders (2011) https://ui.adsabs.harvard.edu/abs/2010ASSP...16..379L/abstract
        #except for C, N, O, P, S, K, and Fe for which values are from Caffau et al. (2011)
        #https://ui.adsabs.harvard.edu/abs/2011SoPh..268..255C/abstract
        sol = [ 12.00, 10.925, 3.28, 1.32, 2.81,  8.50,  7.86,  8.76, 4.44, 8.05, 
                6.29,  7.54,  6.46,  7.53,  5.46,  7.16,  5.25,  6.50,  5.11, 6.31, 
                3.07,  4.93,  3.99,  5.65,  5.50,  7.52,  4.90,  6.22,  4.27, 4.65, 
                3.10,  3.59,  2.32,  3.36,  2.56,  3.28,  2.38,  2.90,  2.20, 2.57, 
                1.42,  1.94, -9.99,  1.78,  1.10,  1.67,  1.22,  1.73,  0.78, 2.09, 
                1.03,  2.20,  1.57,  2.27,  1.10,  2.18,  1.19,  1.60,  0.77, 1.47, 
                -9.99,  0.96,  0.53,  1.09,  0.34,  1.14,  0.49,  0.95,  0.14, 0.94, 
                0.11,  0.73, -0.14,  0.67,  0.28,  1.37,  1.36,  1.64,  0.82, 1.19, 
                0.79,  2.06,  0.67, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99, 0.08, 
                -9.99, -0.52, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99]

    elif (reference is None or reference == 'ags2005'):
        #Asplund, Grevesse and Sauval (2005), basically the same as 
        #Grevesse N., Asplund M., Sauval A.J. 2007, Space Science Review 130, 205
        sol = [  0.911, 10.93,  1.05,  1.38,  2.70,  8.39,  7.78,  8.66,  4.56,  7.84, 
                 6.17,  7.53,  6.37,  7.51,  5.36,  7.14,  5.50,  6.18,  5.08,  6.31, 
                 3.05,  4.90,  4.00,  5.64,  5.39,  7.45,  4.92,  6.23,  4.21,  4.60, 
                 2.88,  3.58,  2.29,  3.33,  2.56,  3.28,  2.60,  2.92,  2.21,  2.59, 
                 1.42,  1.92, -9.99,  1.84,  1.12,  1.69,  0.94,  1.77,  1.60,  2.00, 
                 1.00,  2.19,  1.51,  2.27,  1.07,  2.17,  1.13,  1.58,  0.71,  1.45, 
                 -9.99,  1.01,  0.52,  1.12,  0.28,  1.14,  0.51,  0.93,  0.00,  1.08, 
                 0.06,  0.88, -0.17,  1.11,  0.23,  1.45,  1.38,  1.64,  1.01,  1.13,
                 0.90,  2.00,  0.65, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99,  0.06,   
                 -9.99, -0.52, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99 ]

    else:
        print('NOT a valid reference for the solar composition (ags2005, husser, basti)')
        nans = np.empty((99,))
        nans[:] = np.nan
        return (symbol, mass, nans )
              
    sol[0] = 1.0
    for i in range(len(sol)-1): sol[i+1] = 10.**(sol[i+1]-12.0)

    return (symbol,mass,sol)


  
