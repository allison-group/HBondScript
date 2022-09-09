#!/opt/homebrew/bin/python3

import os,sys,glob
import numpy as np

# get list of donor residue numbers for each hydrogen bond
data_ = np.genfromtxt('bonds.edges',dtype=None,usecols=4,encoding=None)

# back up the unsorted  bonds.edges and bonds.list
os.system('mv bonds.edges bonds.edges.o')
os.system('mv bonds.list bonds.list.o')
# JRA: also back up unsorted bonds.info
os.system('mv bonds.info bonds.info.o')

# sort the list of donor residue numbers - outputs an array of indices that would correspond to the input array being sorted
# and run through that list of indices
for mem in np.argsort(data_):
  c_ = 0
  d_ = 0
  # JRA: added so I can sort my data file
  e_ = 0
  # search through the backed-up bonds.edges file to find the line that corresponds to the current sorted index and write it to file (?)
  with open('bonds.edges.o') as f:
    for line in f:
      if c_ == mem:
        with open('bonds.edges','a') as g:
          g.write(line)
      c_ += 1
  # do the same for bonds.list
  with open('bonds.list.o') as f:
    for line in f:
      if d_ == mem:
        with open('bonds.list','a') as g:
          g.write(line)
      d_ += 1
  # JRA: do the same for bonds.info
  with open('bonds.info.o') as f:
    for line in f:
      if e_ == mem:
        with open('bonds.info','a') as g:
          g.write(line)
      e_ += 1

# JRA: Why not sort the whole array with np.argsort() and write it to file - you can tell it which column to sort against and then use take_along_axis to apply the resulting index_array
# JRA: Why does he treat bonds.edges and bonds.list separately?
