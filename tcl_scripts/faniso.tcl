# Fractional Anisotropy
# Mayank Misra
# mm3897@columbia.edu
# 2 January 2018

# Requires Linear Algebra package from http://www.hume.com/la/la.html
################################################################################
# Purpose: Find the Fractional Anisotropy (FA) of a selection (fragment)

# Example:
# This will loop through each fragment in the molecule ($mol) and find the
# average FA for the system
#
# set f_list [lsort -unique [[atomselect $mol all] get fragment]]
#
# fa = 0.0
# foreach frag in f_list {
# 	set sel [atomselect $mol "fragment $frag"]
# 	set fa [expr $fa + [faniso $sel]]
# 	$sel delete
# }
# set fa [expr $fa / [llength $flist]]
# puts "$fa"
################################################################################


package require La
namespace import -force La::*

proc faniso {mol} {
	set com [measure center $mol]
	set x [ $mol get x ]
	set y [ $mol get y ]
	set z [ $mol get z ]
	set m [ $mol get mass ]
	set Ixx 0
	set Ixy 0
	set Ixz 0
	set Iyy 0
	set Iyz 0
	set Izz 0
	foreach xx $x yy $y zz $z mm $m {
		set mm [expr abs($mm)]
		set xx [expr $xx - [lindex $com 0]]
		set yy [expr $yy - [lindex $com 1]]
		set zz [expr $zz - [lindex $com 2]]
		set rr [expr $xx + $yy + $zz]
		set Ixx [expr $Ixx + $mm*($yy*$yy+$zz*$zz)]
		set Ixy [expr $Ixy - $mm*($xx*$yy)]
		set Ixz [expr $Ixz - $mm*($xx*$zz)]
		set Iyy [expr $Iyy + $mm*($xx*$xx+$zz*$zz)]
		set Iyz [expr $Iyz - $mm*($yy*$zz)]
		set Izz [expr $Izz + $mm*($xx*$xx+$yy*$yy)]
	}
	set I [list 2 3 3 $Ixx $Ixy $Ixz $Ixy $Iyy $Iyz $Ixz $Iyz $Izz]
	mevsvd_br I c
	set l1 [lindex $c 3]
	set l2 [lindex $c 4]
	set l3 [lindex $c 5]
	set l [expr ($l1 + $l2 + $l3)/3.0]
	set p1 [expr ($l1-$l)*($l1-$l)]
	set p2 [expr ($l2-$l)*($l2-$l)]
	set p3 [expr ($l3-$l)*($l3-$l)]
	set p4 [expr sqrt($p1+$p2+$p3)*sqrt(1.5)]
	set p5 [expr sqrt($l1*$l1+$l2*$l2+$l3*$l3)]
	return [expr $p4/$p5]
}
