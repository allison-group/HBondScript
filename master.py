import os
import configparser
import cleanData
import calcOcc
import plotHbonds

config = configparser.ConfigParser()
config.sections()
config.read('config.cfg')

# Parameters
vmd = config['DEFAULT']['vmd']
view = config['DEFAULT']['view']
structureFile = config['SIMULATION']['structureFile']
trajectoryFile = config['SIMULATION']['trajectoryFile']
nFrames = config['SIMULATION']['nFrames']
HBondDist = config['HBONDS']['HBondDist']
HBondAngle = config['HBONDS']['HBondAngle']
bond_files = config['PATH']['bond_files']

##################################
if not os.path.exists(bond_files):
    os.mkdir(bond_files)

### CLEAN UP FILES FROM PREVIOUS RUN TO AVOID DOUBLE-APPENDING ########
# If there's a bonds.dat file from a previous iteration, delete it so we don't continue to append to it
if os.path.exists(bond_files + '/bonds.dat'):
    os.remove(bond_files + '/bonds.dat')
# If bonds1.parsed (created by cleanData.py) exists, delete it
if os.path.exists(bond_files + '/bonds1.parsed'):
    os.remove(bond_files + '/bonds1.parsed')
if os.path.exists(bond_files + '/bonds.edges'):
    os.remove(bond_files + '/bonds.edges')
if os.path.exists(bond_files + '/bonds.info'):
    os.remove(bond_files + '/bonds.info')
if os.path.exists(bond_files + '/bonds.list'):
    os.remove(bond_files + '/bonds.list')
######################################################################

# Use VMD to run calcHbonds.tcl to:
# - load the structure file and trajectory, count the number of frames and write that to frames.txt
# - check that the total number of atoms accounted for in index.dat matches the total number of atoms in the molecule
# - assign the groups of atoms outlined in index.dat to chains
# - write out the first frame of the trajectory as a PDB file (which adds chains if we're working with GROMACS output)
# - loop through pairs of chains, calculate the hydrogen bonds, and append to bonds.dat
os.system(vmd + ' -dispdev text -eofexit -args < calcHbonds.tcl ' + structureFile + ' ' + trajectoryFile + ' ' + nFrames + ' ' + HBondDist + ' ' + HBondAngle)

cleanData.cleanData()
calcOcc.calcOcc()

# If we are visualising the hbonds, do that
if os.path.exists(bond_files + '/bonds.edges') and os.path.exists(bond_files + '/bonds.list') and os.path.exists(bond_files + '/bonds.info'):
    if view == 'True':
        os.system(vmd + ' -e viewHbonds.tcl')
    else:
        plotHbonds.plotHbonds()
