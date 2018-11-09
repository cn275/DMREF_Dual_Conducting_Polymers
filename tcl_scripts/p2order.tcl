# P2 order
# Mayank Misra
# mm3897@columbia.edu
# 3 March 2017

# Requires Linear Algebra package from http://www.hume.com/la/la.html
################################################################################
# Purpose: Calculate P2 bond order for a given ordered list of index of
# of the selection (clist) and cut-off distance (cutl) for the calcualtion
################################################################################


proc p2order {name clist cutl} {

	#loop makes bond vector
		for {set i 2} {$i<$n} {incr i} {
		set j [expr $i-2]
		set a [lindex $clist $j]
		set b [lindex $clist $i]
		set c1 [atomselect top "index $a"]
		set c2 [atomselect top "index $b"]
		set x1 [$c1 get {x}]
		set y1 [$c1 get {y}]
		set z1 [$c1 get {z}]
		set x2 [$c2 get {x}]
		set y2 [$c2 get {y}]
		set z2 [$c2 get {z}]
		$c1 delete
		$c2 delete
		set x($b) [expr $x2 - $x1]
		set y($b) [expr $y2 - $y1]
		set z($b) [expr $z2 - $z1]
		set l [expr sqrt(pow($x($b),2)+pow($y($b),2)+pow($z($b),2))]
		set x($b) [expr $x($b)/$l]
		set y($b) [expr $y($b)/$l]
		set z($b) [expr $z($b)/$l]
	}

	#----------define default value to initial----------
	set a [lindex $clist 0]
	set b [lindex $clist 2]
	set x($a) $x($b)
	set y($a) $y($b)
	set z($a) $z($b)

	set a [lindex $clist 1]
	set b [lindex $clist 3]
	set x($a) $x($b)
	set y($a) $y($b)
	set z($a) $z($b)
	#--------------------

	#calculate bond orientation
	set p2 []
	foreach m $clist {
		set p 0
		set neigh [atomselect top "(name $name pbwithin $cutl of index $m) and not index $m"]
		set nlist [$neigh get index]
		set k [$neigh num]
		$neigh delete

		foreach mm $nlist {
			set ax [expr $x($m)*$x($mm)]
			set ay [expr $y($m)*$y($mm)]
			set az [expr $z($m)*$z($mm)]
			set ad [expr $ax+$ay+$az]
			set pt [expr ((3*pow($ad,2))-1)/2]
			set p [expr $p + $pt]
		}
		lappend p2 [expr $p/$k]
	}
	return [expr ([join $p2 +])/[llength $p2]]
}
