#!/opt/homebrew/bin/python3

import os,sys,glob
from tqdm import tqdm 
import numpy as np
from operator import itemgetter

# Define a class for holding hydrogen bond information
class bonds:
	def __init__(self,d,a,h,f):
		self.d = d # donor
		self.a = a # acceptor
		self.h = h # hydrogen atom in the middle
		self.f = str(f) # chain ID

c = 0
all_ = []
# Read in bonds.dat, extract donor, acceptor and hydrogen atom numbers
# See below for example of formatting of bonds.dat
with open('bonds.dat') as f:
	for line in f:
		d, a, h = "","",""
		try:
      # split by newline, take first entry, then split by '*', take second entry, then split by '} {', take first/second/third entry,
      # then for d and h, split again by brackets and take second or first entry, respectively
			d = line.split('\n')[0].split('*')[2].split('} {')[0].split('{')[1]
			a = line.split('\n')[0].split('*')[2].split('} {')[1]
			h = line.split('\n')[0].split('*')[2].split('} {')[2].split('}')[0]
		except:
      # and if the above doesn't work, that's likely because there was only one hydrogen bond, so split by space instead of brackets
			d = line.split('\n')[0].split('*')[2].split(' ')[0]
			a = line.split('\n')[0].split('*')[2].split(' ')[1]
			h = line.split('\n')[0].split('*')[2].split(' ')[2]
    # frame
		f = line.split('\n')[0].split('*')[1].split('frame_')[1] 
    # define this set of hbond descriptors and add to all_
		a = bonds(d,a,h,f)
		all_.append(a)

# JRA: example of VMD's hbond output (in bonds.dat):
#*Chain:A-Chain:B:frame_4*{3882 3732 3419} {20797 20858 20859} {3883 3735 3420}
#*Chain:A-Chain:B:frame_5*{3787 3871} {20796 20764} {3788 3872}
#*Chain:A-Chain:B:frame_6*3882 20796 3883
# in each bracket (no brackets if only one hbond) are listed the donor, acceptor and hydrogen atoms involved in the hbond(s) formed at that frame

# Open bonds0.parsed for appending
b2 = open('bonds0.parsed','a')
# Loop through the set of hbonds in all_
for mem in all_:
  # Extract the items (why not write to file directly above?)
	d = mem.d.split(' ')
	a = mem.a.split(' ')
	h = mem.h.split(' ')
	f = mem.f.split(' ')
  # Write to file with tab separation; note only writing the first item in f (frame number)
	for i in range(0,len(d)):
		b2.write(d[i] + '\t' + a[i] + '\t' + h[i] + '\t' + f[0] + '\n')
b2.close()

# Read in the data we just wrote to file, fill any gaps with -99999
b3_f = np.genfromtxt('bonds0.parsed',delimiter='\t',filling_values=-99999)
# No idea what this achieves, but b3 isn't used again anyway
#b3 = np.matrix(sorted(b3_f, key=lambda x: x[3]))
# Loop through data just read from file, append each item to idx_
idx_ = []
for mem in b3_f:
	idx_.append(mem[0])
# Write all lines that don't contain gaps (that were filled with -99999) to bonds1.parsed
if max(idx_) != -99999:
	np.savetxt('bonds1.parsed',b3_f, fmt='%-4d',delimiter='\t')
else:
	raise TypeError('No H-bonds found in system')

