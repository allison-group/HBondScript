# Delete any already-loaded molecules
mol delete all
# Load the PDB file of the first frame
mol load pdb temp.pdb

# read in the hbond info
set h_bonds [open bond_files/bonds.edges r]
# draw the hydrogen bonds as red cylinders between the donor and acceptor atoms
draw color red
while { [gets $h_bonds hb_idx] >= 0 } {
	set acc [lindex $hb_idx 2]
	set hyd [lindex $hb_idx 3]
	set p1 [atomselect top "index $acc"]
	set p2 [atomselect top "index $hyd"]
	draw cylinder [lindex [$p1 get {x y z}] 0] [lindex [$p2 get {x y z}] 0] radius [expr [lindex $hb_idx 0]/100] filled yes
}
close $h_bonds
puts "Visualizing hBonds"

