#!/usr/bin/env python3

# Program to read in a dump file frame by frame
# and do some calculations.

# Example here is radius of gyration
# Here Rg of all atoms in the file are calculated

# Assumes the dump file has the following format
# id type x y z ix iy iz
# were x  etc. are atom coordinates, and are numbers between -L/2 and L/2
#      ix etc. are image flags for taking periodic boundaries into account

# Assumes that atom ids start at 1 and are continuous between 1...N

#***************************************************************************
# NOTE : Be careful NOT to use this code for the case of a dump line with  *
#        xs ys zs as in that case the atom coordinates have been scaled    *
#        to lie between 0 and 1                                            *
#***************************************************************************


import numpy as np
import operator


class Atom:
    """ A Class for storing atom information """

    def __init__(self):
        """ Initialise the class """
        self.id = 0                                              # id of the atom
        self.type = 0                                            # type of the atom
        self.L = np.array([0.0,0.0,0.0],dtype=np.float64)        # size of simulation box
        self.x = np.array([0.0,0.0,0.0],dtype=np.float64)        # position of the atom
        self.image = np.array([0,0,0],dtype=np.int32)            # image flags for atoms
        self.x_unwrap = np.array([0.0,0.0,0.0],dtype=np.float64) # position of the atom - unwrapped coords
        self.unwrap_flag = False

        
    def sep(self,B):
        """ Find separation of this Atom and Atom B """
        return np.sqrt( (self.x[0]-B.x[0])**2 + 
                        (self.x[1]-B.x[1])**2 + (self.x[2]-B.x[2])**2 )

    def minus(self,B):
        """ Subtract B.x vector from self.x vector """
        AminusB = np.array([0.0,0.0,0.0],dtype=np.float64)
        for i in range(3):
            AminusB[i] = self.x[i] - B.x[i]
        return AminusB

    def xdot(self,B):
        """ Find dot product of position x of this Atom and Atom B """
        AdotB = np.array([0.0,0.0,0.0],dtype=np.float64)
        # TO DO : find AdotB
        return AdotB

    def unwrap(self):
        """ Unwraps the coordinates for periodic box to generate x_unwrap array """
        if not self.unwrap_flag:   # first check it has not already been done
            for j in range(3):
                self.x_unwrap[j] = self.x[j] + self.image[j]*self.L[j] # unwrap
            unwrap_flag = True



def readframe(infile,N):
    """ Read a single frame of N atoms from a dump file 
        Expects coordinates to be in range -L/2 to L/2
        DOES NOT Unwrap corrdinates for periodic box """

    atoms = [Atom() for i in range(N)]
    L = np.array([0.0,0.0,0.0],dtype=np.float64)

    # read in the 9 header lines of the dump file
    # get box size
    for i in range(9):
        line = infile.readline()
        if i==1:  ## second line of frame is timestep
            timestep = np.int32(line)
        if i==5 or i==6 or i==7:   # 6th-8th lines are box size in x,y,z dimensions
            # get the box size
            line = line.split()
            L[i-5]=np.float64(line[1]) - np.float64(line[0]);

    # now read the atoms, putting them at the correct index (index=id-1)
    for i in range(N):
        line = infile.readline()
        line = line.split()
        index = int(line[0])-1  # LAMMPS atom ids start from 1, python index starts from 0
        atoms[index].id = int(line[0])
        atoms[index].type = int(line[1])
        atoms[index].L = L
        for j in range(3):
            atoms[index].x[j] = np.float64(line[j+2])
        for j in range(3):
            atoms[index].image[j] = np.int32(line[j+5])

    return atoms,timestep


def lines_in_file(filename):
    """ Get the number of lines in the file """

    with open(filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1


def radius_of_gyration(atoms,Npoly):
    """ Calculate the radius of gytation -- Rg^2 = (1/N) sum ( r_k - r_mean)^2 
    remember to unwrap periodic boundaries """

    # requires unwrapped coords, so first calculate these if not already
    for i in range(len(atoms)):
        atoms[i].unwrap()
    
    # get mean position
    r_mean = np.zeros(3,dtype=np.float64)
    for i in range(Npoly):    ## assumes that atoms at index 0...Npoly-1 inclusive are in the polymer
        r_mean += atoms[i].x_unwrap
    r_mean = r_mean/Npoly
    
    # get Rg2
    Rg2 = 0.0
    for i in range(Npoly):
        Rg2 += np.sum( np.square( atoms[i].x_unwrap - r_mean ) )
    Rg2 = Rg2/Npoly
    
    return np.sqrt( Rg2 )




############################################################################
#### Start of the main program

dumpfilename = 'dump.1.DNA' #'dump.DNA'   # this is hardcoded here, could instead read as command line argument
Natoms = 2042   # this is hardcoded here, could instead read as command line argument or read from dump file
                # ALERT: this one is for total atoms
Npoly = 2000    # this code assumes that the first Npoly atoms are in the polymer



Nlines = lines_in_file(dumpfilename)  # get length of file
Nframes = int(Nlines / (Natoms+9))         # there are 9 header lines in each frame

#outfile_Rg = 'r_g_1.dat'  # this is hardcoded here, could instead read as command line argument
inf = open(dumpfilename, 'r')  
# open the intput file
# =============================================================================
# 
# 
# # open the output files and print a header
# ouf_rg = open(outfile_Rg, 'w')  
# ouf_rg.write("# frame number, radius of gyration\n")
# 
# # go through the file frame by frame
# for frame in range(Nframes):
#     # read the frame, unwrapping periodic coordinates
#     atoms, timestep = readframe(inf,Natoms)
#     
#     # unwarp period boundary coordinates -- needed for radius of gyration
#     for i in range(len(atoms)):
#         atoms[i].unwrap()
# 
#     # calculate radius of gyration
#     Rg = radius_of_gyration(atoms,Npoly)
#     
#     # output some results
#     ouf_rg.write( "%i %.5f\n"%(timestep,Rg) )
# 
# print('Done, bro')
# =============================================================================

# close the files

# parameters for output file
outfile_types = 'atom_types.dat'
outfile_positions = 'atom_positions.txt'  # Output file for recording atom positions

# open the output files and print headers
ouf_types = open(outfile_types, 'w')
ouf_types.write("# frame number, atom id, atom type\n")

ouf_positions = open(outfile_positions, 'w')
ouf_positions.write("# frame number, atom id, atom type, x, y, z\n")

# go through the file frame by frame
for frame in range(Nframes):
    # read the frame, unwrapping periodic coordinates
    atoms, timestep = readframe(inf, Natoms)

    # output atom types
    for atom in atoms:
        ouf_types.write(f"{timestep} {atom.id} {atom.type}\n")

    # output atom positions
    for atom in atoms:
        ouf_positions.write(f"{timestep} {atom.id} {atom.type} {atom.x[0]} {atom.x[1]} {atom.x[2]}\n")

# close the files
inf.close()
ouf_types.close()
ouf_positions.close()

print('Done, bro')




# Finished!
