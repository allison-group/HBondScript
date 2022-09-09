# Load up the PDB file of the first frame (written by new_hbond_M.tcl)
mol load pdb temp.pdb
# Select all atoms
set sel [atomselect top all]

# Set up a list of residue numbers
set residue_num [lsort -unique -dict [$sel get resid]]
# Total number of residues
set num_o_res [llength $residue_num]
# Check the residue names to make sure they're all standard 3-letter names (why?)
set res_const_chk 0
# Make a list of residue names and loop through them
set rex [lsort -unique [$sel get resname]]
foreach i $rex {
  # If the length of a residue name is greater than 3, complain and increment res_const_chk
	if {[llength [split $i ""]] > 3} {
		incr res_const_chk
    # NOTE: the above only checks if a residue name is > 3 characters. Is 4 the most that is possible, or just an assumption?
		error "There exist in the system residue names with 4 characters. Adjust these to 3 characters. E.g. fix : $i to [ join [lrange [split $i ""] 0 2] ""]"
	}
}

# Loop through the total number of residues
for {set i 0} {$i < $num_o_res} {incr i} {
	if {$i > 0} {
    # Check that the residue numbers always increment by 1; if not, increment res_const_chk
		if {[expr [lindex $residue_num $i] - [lindex $residue_num [expr $i-1]]] != 1} {
			incr res_const_chk
		}
	}
}
# JRA: No error message is thrown by the residue number incrementation check above; does it matter at all??

# Select all atoms (again...)
set all [atomselect top all]
# Write to PDB file (again...)
$all writepdb temp_CA.pdb
$all delete

# Cleaning stuff
# JRA: Adding a line to the top of temp_CA.pdb that reads 'HEADER' plus the value of res_const_chk - why??
exec echo HEADER    $res_const_chk >> foo.txt
exec cat foo.txt temp_CA.pdb > foo1.txt
exec mv foo1.txt temp_CA.pdb
exec rm foo.txt
