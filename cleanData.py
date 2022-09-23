"""
Run cleanData.py to:
- Read in the hbond data in bonds.dat
- Select the donor and acceptor atoms, whatever h represents, and the chain ID (?)
- Fill any gaps with -99999 and write to bonds0.parsed
- Remove the gaps that were just filled and write to bonds1.parsed
JRA: No idea what the '2' is used for (no arguments are read in). py.dat may contain note that no hbonds are present?
JRA: CAN THIS (the 2 and output) BE BINNED? slave1M.py runs fine without the 2. py.dat is only used for error-reporting that could just go to stdout?
"""
import numpy as np
import configparser

config = configparser.ConfigParser()
config.sections()
config.read('config.cfg')
bond_files = config['PATH']['bond_files']

# Define a class for holding hydrogen bond information
class bonds:
    def __init__(self, d, a, h, f):
        self.d = d  # donor
        self.a = a  # acceptor
        self.h = h  # hydrogen atom in the middle
        self.f = str(f)  # chain ID


def cleanData():
    all_ = []
    # Read in bonds.dat, extract donor, acceptor and hydrogen atom numbers
    # # See below for example of formatting of bonds.dat
    # JRA: example of VMD's hbond output (in bonds.dat):
    # *Chain:A-Chain:B:frame_4*{3882 3732 3419} {20797 20858 20859} {3883 3735 3420}
    # *Chain:A-Chain:B:frame_5*{3787 3871} {20796 20764} {3788 3872}
    # *Chain:A-Chain:B:frame_6*3882 20796 3883
    # in each bracket (no brackets if only one hbond) are listed the donor, acceptor and hydrogen atoms involved in the hbond(s) formed at that frame
    with open(bond_files + '/bonds.dat') as file:
        for line in file:
            d, a, h = "", "", ""
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
            file = line.split('\n')[0].split('*')[1].split('frame_')[1]
            a = bonds(d, a, h, file)
            all_.append(a)

    b2 = open(bond_files + '/bonds1.parsed', 'a')
    for mem in all_:
        d = mem.d.split(' ')
        a = mem.a.split(' ')
        h = mem.h.split(' ')
        f = mem.f.split(' ')
        for i in range(0, len(d)):
            if d[i] == "":
                d[i] = "-99999"
            if a[i] == "":
                a[i] = "-99999"
            if h[i] == "":
                h[i] = "-99999"
            b2.write(d[i] + '\t' + a[i] + '\t' + h[i] + '\t' + f[0]+'\n')
    b2.close()
