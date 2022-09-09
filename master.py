#!/opt/homebrew/bin/python3

import os,sys,glob
import subprocess 

#### SETUP ########################
# Define where to find VMD. Note: you also have to define this in calcOcc.py
vmd = '/Applications/VMD\ 1.9.4a57-arm64-Rev12.app/Contents/vmd/vmd_MACOSXARM64'
# Choose whether to view the hbonds in VMD or not
# Note: You can either view or make the graph. If you choose to view, you can make the graph later using plotHbonds.py
#       If you choose to make the graph, you can view the hbonds later using viewHbonds.tcl
view = False
# Define psf/gro file name
structureFile = 'PI3K/npt_protein_only.gro'
# Define trajectory file name
trajectoryFile = 'PI3K/onlyprot_10000frames.xtc'
##################################

### CLEAN UP FILES FROM PREVIOUS RUN TO AVOID DOUBLE-APPENDING ########
# If there's a bonds.dat file from a previous iteration, delete it so we don't continue to append to it
if os.path.isfile('bonds.dat'):
	os.system('rm bonds.dat')
# If bonds0.parsed (created by cleanData.py) exists, delete it
if os.path.isfile('bonds0.parsed'):
	os.system('rm bonds0.parsed')
# If bonds1.parsed (created by cleanData.py) exists, delete it
if os.path.isfile('bonds1.parsed'):
	os.system('rm bonds1.parsed')
######################################################################

# Use VMD to run calcHbonds.tcl to:
# - load the structure file and trajectory, count the number of frames and write that to frames.txt
# - check that the total number of atoms accounted for in index.dat matches the total number of atoms in the molecule
# - assign the groups of atoms outlined in index.dat to chains
# - write out the first frame of the trajectory as a PDB file (which adds chains if we're working with GROMACS output)
# - loop through pairs of chains, calculate the hydrogen bonds, and append to bonds.dat
os.system(vmd + ' -dispdev text -eofexit -args < calcHbonds.tcl ' + structureFile + ' ' + trajectoryFile)

# Run cleanData.py to:
# - Read in the hbond data in bonds.dat
# - Select the donor and acceptor atoms, whatever h represents, and the chain ID (?)
# - Fill any gaps with -99999 and write to bonds0.parsed
# - Remove the gaps that were just filled and write to bonds1.parsed
# JRA: No idea what the '2' is used for (no arguments are read in). py.dat may contain note that no hbonds are present?
# JRA: CAN THIS (the 2 and output) BE BINNED? slave1M.py runs fine without the 2. py.dat is only used for error-reporting that could just go to stdout?
#os.system('./slave1M.py  2> py.dat')
os.system('./cleanData.py')

# JRA: CAN THIS BE BINNED?
# Write any compaints about no hbonds being found to stdout (why not do this in the first place?)
#line = ""
#with open('py.dat','r') as f:
#	for line in f:
#		line = line.rsplit()[0]
#if line == "TypeError:":
#	raise TypeError('No H-bonds found in system')

# Run calcOcc.py to:
# - run structure_stats.tcl, which:
#   - reads in temp.pdb (first frame in PDB format)
#   - checks for residue names > 3 characters
#   - checks that residue numbers always increment by 1
#   - writes a header line with the total number if bad residues (and bad increments) to temp_CA.pdb (rest is same as temp.pdb)
#   - writes output to vmd.log (including any residue names to fix)
# - writes the complaint written above to stdout but with less info
# - reads in bonds1.parsed, sorts this data into a list of hbonds that in turn contains a list of the frames in which each hbond occurs
# - processes this list to calculate the occupancy of each hbond
# - writes these hbonds and their occupancies to file
# - this is where I altered the output to write bonds.info with residue and atom names instead of bonds.edges
os.system('./calcOcc.py')

# If we are visualising the hbonds, do that
if (view == True):
  if os.path.isfile('bonds.edges') and os.path.isfile('bonds.list') and os.path.isfile('bonds.info'):
	  os.system(vmd + ' -e viewHbonds.tcl')
elif os.path.isfile('bonds.edges') and os.path.isfile('bonds.list') and os.path.isfile('bonds.info'):
  os.system('./plotHbonds.py')
