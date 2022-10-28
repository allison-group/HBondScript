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

# Define the structure file (first argument, .psf/.gro) and trajectory file (second argument, .dcd/.xtc)
set str [lindex $argv 0]
set traj [lindex $argv 1]
# Read number of frames in trajectory (defined in config.cfg)
set nf [lindex $argv 2]
# Define criteria for formation of a hyddrogen bond
set hbd [lindex $argv 3]
set hba [lindex $argv 4]

# GROMACS/CHARMM
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

# If GROMACS...then there is no chain entry in the PDB, so we need to set them up
if {$gx_ eq 1} {
  # Define a list of chain IDs (.gro files don't contain chain IDs)
	set list_chx_ {A B C D E F G H I J K L M N O P Q R S T U V W X Y Z a b c d e f g h i j k l m n o p q r s t u v w 0 1 2 3 4 5 6 7 8 9}
  # JRA: this could be set up in config.cfg, e.g. if we provide the atom numbers there (rather than in
  # index.dat), they could be straight away linked to a chain ID and a more meaningful name (e.g., PI3K)
  #
  # Read in the index file, convert the atom numbers to vmd-style atom indices (-1),
  # define each group of atoms (from the index file) as a chain
	set fp [open index.dat r]
	set lnx_ 0
	while { [gets $fp data] >= 0 } {
    # JRA: space, not tab, separated for transferability 
		set st_ [expr [lindex [split $data "\ "] 0] -1]
		set en_ [expr [lindex [split $data "\ "] 1] -1]
		set sel [atomselect top "index $st_ to $en_"]
		$sel set chain [lindex $list_chx_ $lnx_]
		$sel delete
		incr lnx_
	}
	close $fp
}
# Done with if GROMACS

# JRA: this list could be made above rather than here
# Set up a list of chains
set chx [lsort -unique [[atomselect top "all and not waters and not ions"] get chain]]
# Set progress counter
progress_init [expr [llength $chx]*$nf]

# Loop through pairs of chains
set n_c 0
foreach chain0 $chx {
	foreach chain1 $chx {
		if {$chain0 ne $chain1} {
      # JRA: here we are appending to bonds.dat; safer to write instead?
			set f [open bond_files/bonds.dat a]
      # Loop through frames
			for {set i 0} {$i < $nf} {incr i} {
				progress_tick $n_c
				incr n_c
				set ch0 [atomselect top "chain $chain0" frame $i]
				set ch1 [atomselect top "chain $chain1" frame $i]
        # JRA: criteria now set in config.cfg
        # Write hydrogen bonds between each pair of chains to bonds.dat
				puts $f *Chain:$chain0-Chain:$chain1:frame_$i*[measure hbonds $hbd $hba $ch0 $ch1]
				$ch0 delete
				$ch1 delete
			}
			close $f
		}
	}
}
