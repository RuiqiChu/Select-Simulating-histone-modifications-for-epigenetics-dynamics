#Types: 
#1:unmarked
#2:TF ON(sphere)
#3:Reader
#4:Writer
#5:ATAC OFF
#6:K27ac
#7:ATAC ON
#8:TF OFF

# This simulation uses 2 spheres, 100 RW pairs, with morse potential as morse 12.0 1.75 0.0 1.3 between readers and ATAC atoms which can be changed
# in 'in.setup' file
# To run:

# First 

lmp < in.Equilib_polymer_spheres_dumbells

# This 'unpacks' the initial polymer conditions
# Then

./run_simulation

# This will give a file named dump.1.DNA which can be viewed in VMD, a text file named chosen_polymers which gives the details of swichting
# atoms: the reader, bead and the distance between them
# Then

python process_dump_file

# This will give one text named atom_positions.txt ad another text named atom_types.txt which records the corresponding information from the dump file
# Then

python type_heatmap

# This plots and saves two plots: one for the heatmap which shows the type changes, and another one shows the fraction of ON state during 
# the whole simulation for each atom
# Then

python bonded_fraction

# This calculates the distance between two atoms one by one to find out whether two atoms are bonded
# it gives two outputs: bonded_atoms.txt which has the timestep and bonded two atoms as well as their positions,
# sphere_percentage_data.txt which gives the percentage of bonded atoms in each timestep
# which can be plotted using

python plot_fraction

# Several parameters are hardcoded and should be mannually changed when changing parameters:
# Number of atoms in run_simulation.py line 77-79
# Number of atoms in process_dump_file line 142-144
# Number of atoms in bonded_fraction.py line 29 and 31
# Name and title of plots in type_heatmap.py and plot_fraction.py
