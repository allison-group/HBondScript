# Initiate the timer
proc progress_init {tot} {
   set ::progress_start     [clock seconds]
   set ::progress_last      0
   set ::progress_last_time 0
   set ::progress_tot       $tot
}

# Monitor progress
# We update if there's a 5% difference or a 5 second difference
proc progress_tick {cur} {
   set now [clock seconds]
   set tot $::progress_tot

   if {$cur > $tot} {
       set cur $tot
   }
   if {($cur >= $tot && $::progress_last < $cur) ||
       ($cur - $::progress_last) >= (0.05 * $tot) ||
       ($now - $::progress_last_time) >= 5} {
       set ::progress_last $cur
       set ::progress_last_time $now
       set percentage [expr round($cur*100/$tot)]
       set ticks [expr $percentage/2]
       if {$cur == 0} {
           set eta   ETA:[format %7s Unknown]
       } elseif {$cur >= $tot} {
           set eta   TOT:[format %7d [expr int($now - $::progress_start)]]s
       } else {
           set eta   ETA:[format %7d [expr int(($tot - $cur) * ($now - $::progress_start)/$cur)]]s
       }
       set lticks [expr 50 - $ticks]
       set str "[format %3d $percentage]%|[string repeat = $ticks]"
       append str "[string repeat . $lticks]|[format %8d $cur]|$eta\r"
       puts -nonewline stdout $str
       if {$cur >= $tot} {
           puts ""
       }
       flush stdout
   }
}

########################

# Define the structure file (.psf/.gro) and trajectory (.dcd/.xtc)
# Structure file is the first argument, trajectory is the second argument
set str [lindex $argv 0]
set traj [lindex $argv 1]
set gx_ 0
set ch_ 0
# Is it GROMACS? And if so, load .gro file.
if {[lindex [split $str .] 1] eq "gro"} {
	mol load gro $str
	incr gx_
# Or is it CHARMM/NAMD? And if so, load .psf file.
} elseif {[lindex [split $str .] 1] eq "psf"} {
	mol load psf $str
	incr ch_
}
# Load trajectory
mol addfile $traj waitfor all
# Count number of frames and write to file
set nf [molinfo top get numframes]
# exec echo $nf > frames.txt

# If GROMACS...
if {$gx_ eq 1} {
  # Define a list of chain IDs (.gro files don't contain chain IDs)
	set list_chx_ {A B C D E F G H I J K L M N O P Q R S T U V W X Y Z a b c d e f g h i j k l m n o p q r s t u v w 0 1 2 3 4 5 6 7 8 9}
	############################
	# Read in the index file that lists the first and last atom of the two "chains" between which hydrogen bonds are to be calculated
  # NOTE: according to the "sanity check" below, these two chains must cover all atoms in the system
	set fp [open index.dat r]
	set count 0
	set lnx_ 0
  while { [gets $fp data] >= 0 } {
		incr lnx_
    # JRA: changed from tab to space seapration as tabs aren't consistent between systems
		set st_ [lindex [split $data "\ "] 0]
		set en_ [lindex [split $data "\ "] 1]
		set diff [expr $en_ - $st_ + 1]
		set count [expr $count + $diff]
	}
	close $fp
	############################
	# Sanity check: does the sum of the differences between the first and last atoms defining the chains (count) equal the total number of atoms in the system?
  # JRA: [Why should it? Is it not possible to only calculate hydrogen bonds for a subset of the system? And why not just take the final atom number supplied above rather than calculate all these differences and sums?]
	set sel_all [atomselect top all]
	set all_idx [llength [$sel_all get index]]
	if {$all_idx ne $count} {
		error 'Some atom indices are unaccounted for. Make sure your index.dat rows account for all atom indices. $all_idx atoms in system. $count atoms in index'
	}
	############################
  # Read in the index file again (?!), this time convert the atom numbers to vmd-style atom indices (-1),
  # define each group of atoms (from the index file) as a chain
	set fp [open index.dat r]
	set lnx_ 0
	set diff_c 0
	while { [gets $fp data] >= 0 } {
    # JRA: same here
		set st_ [expr [lindex [split $data "\ "] 0] -1]
		set en_ [expr [lindex [split $data "\ "] 1] -1]
		set sel [atomselect top "index $st_ to $en_"]
		$sel set chain [lindex $list_chx_ $lnx_]
		$sel delete
		incr lnx_
	}
	close $fp
	############################
}
# Done with if GROMACS

#  Write the first frame of the trajectory to file in PDB format - for PSF+DCD this just extracts a single frame, for GRO+XTC this has added chain IDs
set sel [atomselect top "all and not waters and not ions" frame 0]
$sel writepdb temp.pdb
$sel delete

# Set up a list of chains
set chx [lsort -unique [[atomselect top "all and not waters and not ions"] get chain]]
# Set progress counter
progress_init [expr [llength $chx]*$nf]

# Loop through pairs of chains
set n_c 0
foreach chain0 $chx {
	foreach chain1 $chx {
		if {$chain0 ne $chain1} {
			set f [open bond_files/bonds.dat a]
            # Loop through frames
			for {set i 0} {$i < $nf} {incr i} {
				progress_tick $n_c
				incr n_c
				set ch0 [atomselect top "chain $chain0" frame $i]
				set ch1 [atomselect top "chain $chain1" frame $i]
                # Write hydrogen bonds between chain0 and chain1 to the file bonds.dat. Criteria are 3 A and 20 deg.
				puts $f *Chain:$chain0-Chain:$chain1:frame_$i*[measure hbonds 3 20 $ch0 $ch1]
				$ch0 delete
				$ch1 delete
			}
			close $f
		}
	}
}
