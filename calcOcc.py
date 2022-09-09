"""
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
# - this is where I altered the output to write bond_files/bonds.info with residue and atom names instead of bonds.edges
"""
from __future__ import division
import os
from tqdm import tqdm 
import numpy as np
import configparser

config = configparser.ConfigParser()
config.sections()
config.read('config.cfg')
vmd = config['DEFAULT']['vmd']

bond_files = config['PATH']['bond_files']

# Define a class that contains information about the hydrogen bonds themselves
# JRA: added rn (residue name) and an (atom name)
class edges:
  def __init__(self,occ,frames,list_H,d,a,h,r1,rn1,an1,c1,r2,rn2,an2,c2):
    self.occ = occ
    self.frames = frames
    self.list_H = list_H
    self.d = d
    self.a = a
    self.h = h
    self.r1 = r1
    self.rn1 = rn1
    self.an1 = an1
    self.c1 = c1
    self.r2 = r2
    self.rn2 = rn2
    self.an2 = an2
    self.c2 = c2

def calcOcc():
    # JRA: changed so VMD is called according to the definition above
    os.system(vmd + ' -dispdev text -eofexit < structure_stats.tcl > vmd.log')

    # Read in temp_CA.pdb, extract the second item in the first line (after that it bails).
    # This is the value of res_const_chk, written to temp_CA.pdb by structure_stats.tcl
    # a non-zero value indicates the number of residues that have names > 3 characters (or instances of increments between residue numbers > 1)
    # JRA: deleting this chunk of code doesn't seem to affect the overall running of the program
    # and four-letter residue names don't seem to stop VMD analysing the hydrogen bonds (it still produces output...but may miss hbonds for any four-letter residues?)
    #c = 0
    #with open('temp_CA.pdb') as f:
    #  for line in f:
    #    c += 1
    #    line = line.split('\n')[0].split(' ')[1]
    #    if c == 1:
    #			break
    #if int(line) > 0:
    #	raise ValueError('The structure is missing a residue -- or has improper residue names see vmd.log in working directory for more info.')

    # Read in the atom numbers (idx), residue numbers (res) and chain ID (chx)
    vals_idx = np.genfromtxt('temp.pdb',usecols=[1],dtype=None,skip_header=1,skip_footer=1,delimiter=[6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,10,2,2],encoding=None)
    vals_res = np.genfromtxt('temp.pdb',usecols=[8],dtype=None,skip_header=1,skip_footer=1,delimiter=[6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,10,2,2],encoding=None)
    vals_chx = np.genfromtxt('temp.pdb',usecols=[7],dtype=None,skip_header=1,skip_footer=1,delimiter=[6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,10,2,2],encoding=None)
    # JRA: Added code to read in the residue names (rname) and atom names (aname)
    vals_rname = np.genfromtxt('temp.pdb',usecols=[5],dtype=None,skip_header=1,skip_footer=1,delimiter=[6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,10,2,2],encoding=None)
    vals_aname = np.genfromtxt('temp.pdb',usecols=[3],dtype=None,skip_header=1,skip_footer=1,delimiter=[6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,10,2,2],encoding=None)
    all_e = []

    # Import hbond info
    data = np.genfromtxt(bond_files + '/bonds1.parsed',usecols=[0,1,2])
    # import frames
    frames = np.genfromtxt(bond_files + '/bonds1.parsed',usecols=3)
    # unique-ofy frames
    fr = np.unique(frames)
    # Iterate over frames
    # JRA: deleted as not used
    #c = 0
    #r1r2 = []
    for mem in tqdm(fr):
    # hbonds in current frame
      hbonds = data[np.where(mem == frames)[0]]
      for bonds_h in hbonds:
        d = bonds_h[0]
        a = bonds_h[1]
        h = bonds_h[2]
        # Checking for entries that imply gaps in the original hbond data (but these should have been removed in
        # going from bonds0.parsed to bonds1.parsed)
        if d != -99999 and a != -99999 and h != -99999:
          r1 , r2 =  vals_res[np.where(vals_idx == d+1)[0][0]] , vals_res[np.where(vals_idx == a+1)[0][0]]
          c1 , c2 =  vals_chx[np.where(vals_idx == d+1)[0][0]] , vals_chx[np.where(vals_idx == a+1)[0][0]]
          # JRA: Added two lines to get residue names (rn1, rn2 for donor and acceptor) and atom names (an1,
          # an2 for donor and acceptor)
          rn1 , rn2 = vals_rname[np.where(vals_idx == d+1)[0][0]], vals_rname[np.where(vals_idx == a+1)[0][0]]
          an1 , an2 = vals_aname[np.where(vals_idx == d+1)[0][0]], vals_aname[np.where(vals_idx == a+1)[0][0]]

          # JRA: Ashar commented out the line below
                #def __init__(self,occ,frames,list_H,d,a,h,r1,c1,r2,c2):
                # check current values of all_e
          flag = 0
          # if there are entries in all_e
          if len(all_e) > 1:
            # check if the donor and acceptor atoms and chain match those of the hbond currently being processed,
            # i.e., has this hbond already been encountered in a previous frame
            for it in range(0,len(all_e)):
              if d == all_e[it].d and a == all_e[it].a and h == all_e[it].h:
                            # increment list_H
                            #nidx_ = np.where(r1r2 == [r1,r2])[0]
                # if they do, add the frame number (?) to the list_H column - ah, this is a list of the frames in which this hbond occurs
                all_e[it].list_H.append(mem)
                            #print nidx_,it
                flag = 1
          # if there are no entries in all_e (and the above section hasn't been run)
          if flag == 0:
            # JRA: Ashar comment out the line below
                    #r1r2.append([r1,r2])
            # Add the current hbond to the list of hbonds (edges)
            # JRA: altered to include rn and an
            # JRA: note that the first two entries, occupancy and (total number of frames) are set to zero - these are filled later
            # JRA: list_H is set to the current frame being processed - if this hbond is encountered again, the section above is run and that frame number appended
            # JRA: hence list_H is a list of the frames in which this Hbond occurs
            #all_e.append(edges(0,0,[mem],d,a,h,r1,c1,r2,c2))
            all_e.append(edges(0,0,[mem],d,a,h,r1,rn1,an1,c1,r2,rn2,an2,c2))

    # Add frames, occ, and make log
    # This file lists key features of each hbond: occupancy, donor atom, acceptor atom, h(?), donor {residue number, [residue name], [atom name], chain}, acceptor {residue number, [residue name], [atom name], chain}; where the items in square brackets were added by me to make bond_files/bonds.info. It's used to create the y-axis tick labels.
    f = open(bond_files + '/bonds.edges','a')
    # This file lists the frames in which each hbond occurs. It's used for placing the dots onto the eventual hbond plot.
    g = open(bond_files + '/bonds.list','a')
    # JRA: Writing a new file with the additional information we need for better y-axis labels
    h = open(bond_files + '/bonds.info','a')
    # Loop through entries in all_e (list of all edges/hbonds)
    for it in range(0,len(all_e)):
      # set frames entry to number of unique frames
        all_e[it].frames = len(fr)
      # calculate % occupancy (list_H is number of frames in which this hbond occurs), fill this entry
        all_e[it].occ = round(len(all_e[it].list_H)*100/all_e[it].frames,2)
        # Writing out minimal edge data
        f.write("{:<6} {:<7} {:<7} {:<7} {:<7} {:<3} {:<7} {:<3}\n".format(\
            all_e[it].occ,\
            all_e[it].d,all_e[it].a,all_e[it].h,\
            all_e[it].r1,all_e[it].c1,all_e[it].r2,all_e[it].c2))
        for frx in all_e[it].list_H:
            g.write("{:<9} ".format(frx))
        g.write('\n')
      # JRA: Writing out the additional information
        h.write("{:<6} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7} {:<3} {:<7} {:<7} {:<7} {:<3}\n".format(\
            all_e[it].occ,\
            all_e[it].d,all_e[it].a,all_e[it].h,\
            all_e[it].r1,all_e[it].rn1,all_e[it].an1,all_e[it].c1,\
        all_e[it].r2,all_e[it].rn2,all_e[it].an2,all_e[it].c2))
    f.close()
    g.close()
    h.close()

    print("Bond Statistics:")
    data = np.genfromtxt(bond_files + '/bonds.edges',usecols=0)
    print("0--25:",((0 <= data) & (data < 25)).sum(),"|| 25--50:",((25 <= data) & (data < 50)).sum(),"|| 50--75:",\
        ((50 <= data) & (data < 75)).sum(),"|| 75--100:",((75 <= data) & (data <= 100)).sum()," || Total = ",len(data))

