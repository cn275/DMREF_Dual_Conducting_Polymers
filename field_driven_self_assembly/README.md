The files in this folder are for use with the LAMMPS software package. Modified from the fix_addforce codes (17Nov16 release). This fix is intended to divide the simulation box into a grid (spacing is specified by user) to generate a grid of cells. The user can import a linear table which specifies a force vector to be applied to atoms when they are inside a given unit cell. 
This is an extension of previous work where an array of "ghost" particles was used with gaussian pair potentials to apply forces to atoms in a molecular dynamics simulation.
These forces can be turned on/off to allow for the second order self assembly of molecules in the shape of the provided field. See https://doi.org/10.1021/acs.jctc.8b00419 for full paper

This formulation removes the need for ghost particles and instead a table can be used, much like tabular pair potentials.

an example implementation:

fix     fix_name group addforcefield 0 0 0 file pot_file.txt pot scale_factor scale_value

where:
fix_name=the fix name
group=the group which the fix is applied to
addforcefield=this fixes specific name
0, 0, 0 = placeholders (leave as is)
file = specifying the file option
pot_file.txt=name of forcefield file
pot = name of potential inside the forcefield file to use (same convention as tabular potentials)
scale_factor = specify the scale factor option
scale_value = the factor to scale the field by

to create a forcefield file examine the example potential provided. The header of each potential has a name followed by:

N n NX nx NY ny NZ nz
n=number of total grid cells
ni=number of cells along i axis

The forcefield is required to be a linear table for ease of parsing, and computation.

The position of cell nc can be calculated as follows:
x=nc%nx*xspacing
y=int((nc%(nx*ny))/nx)*yspacing
z=int(nc/(nx*ny))*zspacing

the spacings are simply the simulation box dimensions divided by th number of grid points along their respective dimensions

The columns in the potential file are, from left to right: grid point, placeholder, Energy, Force, fx, fy, fz

With this formulation the computations are approximately 4 times faster and are as adjustable as the original ghost particle formulation.



