import os
import sys
import fileinput
import platform
import subprocess
from setuptools.command.install import install
from setuptools import setup, find_packages

from platform import system as current_platform
import shutil
from glob import glob

def get_bin_path():
    # Get environment bin/ directory
    bindir = None
    # check if the --user option was set
    userdir = None
    for u in install.user_options:
        if u[0]=='user':
            # "install in user site-package '/Users/nidever/.local/lib/python3.7/site-packages'"
            uline = u[2]
            userdir = uline[uline.find('site-package')+13:]
            userdir = userdir.replace("'","")
    if userdir is not None:
        # /Users/nidever/.local/bin
        bindir = os.path.dirname(os.path.dirname(os.path.dirname(userdir)))+'/bin'
        if os.path.exists(bindir)==False:
            bindir = False
    # Try virtual environment using sys.prefix
    if bindir is None:
        venv = get_virtualenv_path()
        if venv is not None:
            bindir = venv+'/bin'
            if os.path.exists(bindir)==False:
                bindir = None
    # Get bin/ directory from python executable
    if bindir is None:
        out = subprocess.run(['which','python'],shell=True)
        if type(out) is bytes:
            out = out.decode()
        bindir = os.path.dirname(out)
        if os.path.exists(bindir)==False:
            bindir = None

    if bindir is None:
        raise Exception('No bin/ directory found')

    return bindir


def compilesynspec():
    """ Compile the Synspec Fortran code."""
    src_path = './src/'

    # Identify the platform
    platform = current_platform()
    # Check for platform first
    #if platform not in ('Darwin', 'Linux'):
    #    sys.stderr.write("Platform '%s' not recognised!\n" % platform)
    #    sys.exit()
        
    # Install the software
    print('Compiling the Fortran code')
    ret = subprocess.run(['make'], cwd=src_path, shell=True)

    # Get the path for the binaries
    bindir = get_bin_path()
        
    # Copy fortran binaries to bin/ directory
    for f in ['synspec54','rotin','list2bin']:
        if os.path.exists('bin/'+f):
            if os.path.exists(bindir+f): os.remove(bindir+f)
            print('Copying bin/'+f+' -> '+bindir+'/'+f)
            shutil.copyfile('bin/'+f,bindir+'/'+f)
            # Make executable
            os.chmod(bindir+'/'+f,0o755)
        else:
            print('bin/'+f+' NOT FOUND')

    # Download the linelists and convert to binary
    print('Downloading linelists')
    ret = subprocess.run(['make'], cwd='./python/synspec/linelists/', shell=True)
    

# We need to build Synspec
if 'install' in sys.argv or 'develop' in sys.argv or 'bdist_wheel' in sys.argv:
    compilesynspec()
    
setup(name='synspec',
      version='1.0.3',
      description='Synspec spectral synthesis code and python wrapper',
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/dnidever/synspec',
      requires=['numpy','astropy(>=4.0)','scipy','matplotlib','dlnpyutils'],
      include_package_data=True,
      packages = ['synspec'],
      package_dir={"": "python"}   
)
