from __future__ import division
import sys, os
import numpy as np
import matplotlib.pyplot as pl
import configparser

config = configparser.ConfigParser()
config.sections()
config.read('config.cfg')
bond_files = config['PATH']['bond_files']
frames = int(config['SIMULATION']['nFrames'])


def plotHbonds():
    # JRA: Added reading of bond_files/bonds.info, removed reading of bonds.edges with dtype=None
    # This is a new file I added to the output with more information to allow the y-axis tick labels to be changed
    data = np.genfromtxt(bond_files + '/bonds.info', dtype=None, encoding='utf-8')
    # Read in occupancy data
    data_occ = np.genfromtxt(bond_files + '/bonds.info', usecols=0, encoding=None)
    # Read in list of donor atom numbers
    data_don = np.genfromtxt(bond_files + '/bonds.info', usecols=4, encoding=None)
    # Read in list of acceptor atom numbers
    data_acc = np.genfromtxt(bond_files + '/bonds.info', usecols=8, encoding=None)

    # A list for storing the y-axis tick labels
    tick_ = []

    # JRA: easily change the x-axis label and tick sizes and the y-axis tick size
    xlabelFontsize = 12
    xtickFontsize = 12
    ytickFontsize = 12
    # JRA: set the desired figure size (for consistency between plots)
    xsize = 12
    ysize = 6

    # Select plotting protocol
    method = input("\nThreshold Occ or Residue Num. Enter 'Occ' or 'Num':")

    # If 'occupancy' is chosen, request and read in minimum and maximum occupancy to plot
    if method == "Occ":
        thresh = input("\nThreshold Occ by % (minimum and maximum occupancy to plot, e.g., 40 50): ").split(' ')
        # If only one value provided, use it as the minimum occupancy and set the maximum to 100%
        if len(thresh) == 1:
            print(
                "Only one value provided. Using it as the minimum occupancy and setting the maximum occupancy to 100%")
            t_min = float(thresh[0])
            t_max = 100.0
        # If two values are provided, use the smallest (regardless of order) as the minimum and the largest as the maximum
        elif len(thresh) == 2:
            t_min = min(float(thresh[0]), float(thresh[1]))
            t_max = max(float(thresh[0]), float(thresh[1]))
        # If more than two values are provided, use only the first two as outlined above
        else:
            print("More than two values provided. Taking only the first two.")
            t_min = min(float(thresh[0]), float(thresh[1]))
            t_max = max(float(thresh[0]), float(thresh[1]))

        # JRA: Set the figure size
        pl.figure(figsize=(xsize, ysize))

        # Each row of bond_files/bonds.list contains a list of the frames in which the hbond referenced by the row number exists
        with open(bond_files + '/bonds.list', 'r') as f:
            # cnt counts lines in f; c counts instances of the occupancy meeting the criteria
            cnt, c = 0, 0
            # loop through the lines of bond_files/bonds.list
            for line in f:
                # map the floats in line to integers representing the frames for which this hbond occurs
                # note that in python3, map no longer returns a list so this has to be done explicitly as later functions expect a list
                line = list(map(int, map(float, line.rsplit())))
                # Check occupancy data read in from bonds.edges at the start against occupancy criteria to see whether we plot this hbond
                if data_occ[cnt] >= t_min and data_occ[cnt] <= t_max:
                    # Initialise y as an array of '0.5s', with the number of entries equal to the number of items in line (number of frames in which this hbond occurs)
                    y = np.zeros(len(line)) + 0.5
                    # plot the array of 0.5's created above (y), incremented by the index of the hbond we're looking at (+c), against the list of frames in which this hbond occurs ('line')
                    pl.plot(line, y + c, '.')
                    # JRA: Ashar's y-axis tick labels
                    # tick_.append(str(data[cnt][4]) + ':' + str(data[cnt][5]) + ':' + str(int(data[cnt][1])) + '-' +\
                    #	 str(data[cnt][6]) + ':' + str(data[cnt][7]) + ':' + str(int(data[cnt][2])))
                    # JRA: Potential tick label items
                    DonResName = str(data[cnt][5])
                    DonResNum = str(int(data[cnt][4]))
                    DonAtomName = str(data[cnt][6])
                    DonAtomNum = str(int(data[cnt][1]))
                    DonChain = str(data[cnt][7])
                    AccResName = str(data[cnt][9])
                    AccResNum = str(int(data[cnt][8]))
                    AccAtomName = str(data[cnt][10])
                    AccAtomNum = str(int(data[cnt][2]))
                    AccChain = str(data[cnt][11])
                    # JRA: Choose which items you want to use and their separators
                    tick_.append(
                        DonResName + ' ' + DonResNum + '-' + DonAtomName + ' : ' + AccResName + ' ' + AccResNum + '-' + AccAtomName)
                    c += 1
                cnt += 1
        # Place ticks at a spacing of 1, starting from 0.5, for as many ticks as there are entries in tick_, make the labels the entries in tick_, and set the font size to 8
        pl.yticks(np.arange(0.5, len(tick_), 1), tick_, fontsize=ytickFontsize)
        # JRA: Specify the xtick font size
        pl.xticks(fontsize=xtickFontsize)
        # JRA: Ashar set the x-axis limits to -100 and 100 more than the number of frames (why?!), and the y-axis limits to 0 and the number of entries in tick_
        # pl.axis([-100,frames+100,0,len(tick_)])
        # JRA: I set the x-axis limits to 0 and the number of frames, and the y-axis limits to 0 and the number of entries in tick_
        pl.axis([0, frames, 0, len(tick_)])
        # JRA: Add an x-axis label
        pl.xlabel('Time (ns)', fontsize=xlabelFontsize)
        # JRA: Ensure the tick and axis labels fall within the figure area
        pl.tight_layout()
        # Show the plot
        pl.show()

    # If 'number' is chosen
    elif method == 'Num':
        print("Atoms from the following residues were found interacting. Select one of these to plot its details.")
        print("Don:", np.unique(data_don))
        print("Acc:", np.unique(data_acc))
        # Ask for and read in chosen residue number
        resid = float(input("\nResid: (provide one number from above list)"[0]))
        # Get the indices of the times when this residue occurs in the donor and acceptor list
        # NOTE: the final [0] does nothing as far as I can tell
        val_at = np.where(data_don == resid)[0]
        val_at2 = np.where(data_acc == resid)[0]
        # Ashar initialised both counters here. But cnt needs to be re-initialised for each item in val_at or val_at2
        # cnt,c = 0,0
        c = 0

        # JRA: Set the figure size
        pl.figure(figsize=(xsize, ysize))

        # Loop through the indices for when this residue is a donor
        for mem in val_at:
            # reset cnt here instead
            # cnt counts lines in bond_files/bonds.list, and we need to check each line in bond_files/bonds.list for each item in val_at
            cnt = 0
            with open(bond_files + '/bonds.list', 'r') as f:
                for line in f:
                    line = list(map(int, map(float, line.rsplit())))
                    # If the current value of cnt is equal to one of the donor indices, plot the hbond timeseries
                    if cnt == mem:
                        y = np.zeros(len(line)) + 0.5
                        pl.plot(line, y + c, '.')
                        # JRA: Ashar's y-axis tick labels
                        # tick_.append(str(data[cnt][4]) + ':' + str(data[cnt][5]) + ':' + str(int(data[cnt][1])) + '-' +\
                        #	 str(data[cnt][6]) + ':' + str(data[cnt][7]) + ':' + str(int(data[cnt][2])))
                        # JRA: Potential tick label items
                        DonResName = str(data[cnt][5])
                        DonResNum = str(int(data[cnt][4]))
                        DonAtomName = str(data[cnt][6])
                        DonAtomNum = str(int(data[cnt][1]))
                        DonChain = str(data[cnt][7])
                        AccResName = str(data[cnt][9])
                        AccResNum = str(int(data[cnt][8]))
                        AccAtomName = str(data[cnt][10])
                        AccAtomNum = str(int(data[cnt][2]))
                        AccChain = str(data[cnt][11])
                        # JRA: Choose which items you want to use and their separators
                        tick_.append(
                            DonResName + ' ' + DonResNum + '-' + DonAtomName + ' : ' + AccResName + ' ' + AccResNum + '-' + AccAtomName)
                        # And increment c - this increments the position on the y-axis
                        c += 1
                    # Increment cnt (essentially counts through lines in bond_files/bonds.list)
                    cnt += 1

        # Reset cnt (but not c, as we want to continue moving up the y-axis)
        for mem in val_at2:
            # reset cnt here instead
            cnt = 0
            with open(bond_files + '/bonds.list', 'r') as f:
                for line in f:
                    line = list(map(int, map(float, line.rsplit())))
                    if cnt == mem:
                        y = np.zeros(len(line)) + 0.5
                        pl.plot(line, y + c, '.')
                        # JRA: Ashar's y-axis tick labels
                        # tick_.append(str(data[cnt][4]) + ':' + str(data[cnt][5]) + ':' + str(int(data[cnt][1])) + '-' +\
                        #	 str(data[cnt][6]) + ':' + str(data[cnt][7]) + ':' + str(int(data[cnt][2])))
                        # JRA: Potential tick label items
                        DonResName = str(data[cnt][5])
                        DonResNum = str(int(data[cnt][4]))
                        DonAtomName = str(data[cnt][6])
                        DonAtomNum = str(int(data[cnt][1]))
                        DonChain = str(data[cnt][7])
                        AccResName = str(data[cnt][9])
                        AccResNum = str(int(data[cnt][8]))
                        AccAtomName = str(data[cnt][10])
                        AccAtomNum = str(int(data[cnt][2]))
                        AccChain = str(data[cnt][11])
                        # JRA: Choose which items you want to use and their separators
                        tick_.append(
                            DonResName + ' ' + DonResNum + '-' + DonAtomName + ' : ' + AccResName + ' ' + AccResNum + '-' + AccAtomName)
                        c += 1
                    cnt += 1
        # Place ticks at a spacing of 1, starting from 0.5, for as many ticks as there are entries in tick_,
        # make the labels the entries in tick_, and set the font size
        pl.yticks(np.arange(0.5, len(tick_), 1), tick_, fontsize=ytickFontsize)
        # JRA: Specify the xtick font size
        pl.xticks(fontsize=xtickFontsize)
        # JRA: Ashar set the x-axis limits to -100 and 100 more than the number of frames (why?!), and the y-axis limits to 0 and the number of entries in tick_
        # pl.axis([-100,frames+100,0,len(tick_)])
        # JRA: I set the x-axis limits to 0 and the number of frames, and the y-axis limits to 0 and the number of entries in tick_
        pl.axis([0, frames, 0, len(tick_)])
        # JRA: Add an x-axis label
        pl.xlabel('Time (ns)', fontsize=xlabelFontsize)
        # JRA: Ensure the tick and axis labels fall within the figure area
        pl.tight_layout()
        # Show the plot
        pl.show()
