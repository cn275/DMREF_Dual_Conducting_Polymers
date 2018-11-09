# Chain Vector Direction
# Mayank Misra
# mm3897@columbia.edu
# 25 December 2017

# Requires Linear Algebra package from http://www.hume.com/la/la.html
################################################################################
# Purpose: Find the unit vector towards which a selection (fragment) is
# pointing
#
# Usage:
# set vector [axisvec [atomselect $mol "fragment $f"]]
################################################################################

package require La
namespace import -force La::*

proc axisvec {mol} {
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
	return [[lindex $I 5] [lindex $I 8] [lindex $I 11]]
}
