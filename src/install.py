import os, sys, commands, glob, getopt, shutil

# Start installation
print "Installing mechafm (Mechanical AFM)"
print "- Python implementation build from Hapala et al., PRB 90:085421 (2014)"

# Determine python version and where we are
pythonexec  = sys.executable
mechafminstdir = os.getcwd()

# Programlist
proglist = ['mechafm']
binnames = ['mechafm-python']

def error(errorstr):
    if errorstr!='':
        print 'Error! '+errorstr+'\n'
    print 'Usage: /path/to/python install.py'
    sys.exit(1)

args = sys.argv[1:]
numpy, scipy, matplotlib = True, True, True

# Test if numpy exists
if numpy:
    try:
        from numpy import *
    except:
        error("This program requires numpy, which could not be imported. It probably is not installed (it's not part of the standard Python distribution). See the Python site (http://www.python.org) for information on downloading source or binaries.")

# Test if numpy exists
if scipy:
    try:
        from scipy import *
    except:
        error("This program requires scipy, which could not be imported. It probably is not installed (it's not part of the standard Python distribution). See the Python site (http://www.python.org) for information on downloading source or binaries.")

# Test if numpy exists
if matplotlib:
    try:
        from matplotlib import *
    except:
        error("This program requires matplotlib, which could not be imported. It probably is not installed (it's not part of the standard Python distribution). See the Python site (http://www.python.org) for information on downloading source or binaries.")

# Check if python scripts exist
print "Installing scripts..."
for script in range(len(proglist)):
    if not glob.glob("%s.py" % proglist[script]):
        error("No such %s.py script!" % proglist[script])

# Create bin directory
bindir = "bin"
if not glob.glob(bindir):
    os.mkdir(bindir)

# Create bash scripts
for script in range(len(proglist)):
    binprog = os.path.join(bindir,binnames[script])
    f = open(("%s" % binprog),"w")
    f.write("#!/bin/bash\n")
    realprog = os.path.join(mechafminstdir,("%s.py" % proglist[script]))
    f.write("%s %s $@" % (pythonexec, realprog))
    f.close()
    os.chmod(binprog, 0755)

# Finish installation
print "NOTE: Add %s/bin to your path!" % mechafminstdir
print "Installation completed..."



