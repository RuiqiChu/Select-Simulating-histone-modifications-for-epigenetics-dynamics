#!/bin/env python3

### Driver program for LAMMPS with switching proteins
### Runs in serial, but still uses mpi4py so
### parallelisation later should be easy

### MPI is 'message passing interface' for parallel programming
### mpi4py is the MPI python package

### In MPI the same code is run on each processor.
### For lammps, each processor must pass the same commands to lammps.
### Each processor gets assigned an integer, stored in variable 'me'
### Any 'random number' stuff, we only want to do one one processor
### (because each processor gets a different seed, so would generate different
### random numbers) and then share the results with the other processors
### So we do 'if me==0' for that code

### If running only on a single processors, me=0 always
### Currently I am assuming single processor, so have commented out
### the 'share the results with the other processors' code

### We will assume
#         Bead ids 1 to Npoly are polymer
#         Bead ids Npoly+1 to Npoly+Nprot are proteins
#         Proteins switch randomly back and forward between two types with rate k

###############################################################################
###############################################################################

import sys, getopt
import random
import numpy
import array
from mpi4py import MPI
from lammps import lammps
import ctypes

### -------- Functions here



### ---------------------------------------------------------------------------
### --------- Main program ----------------------------------------------------
### ---------------------------------------------------------------------------

# set up MPI (it was initialised when the library was imported) 
comm = MPI.COMM_WORLD 
me = comm.Get_rank() 
nprocs = comm.Get_size()


### ----- Parse command line --------------------------------

# defaults
runnumber = 1  # an id number for input and output files
quiet = False  # flag to switch off screen output

try:
    opts, args = getopt.getopt(sys.argv[1:],"qn:",)
except getopt.GetoptError:
    print('run_simulation.py [-n <runnumber>] [-q]')
    print('      -n <runnumber>  specifies a run number, default 1.')
    print('      -q              runs quietly without output to screen')
    sys.exit(1)
for opt, arg in opts:
    if (opt=="-n"):
        runnumber = int(arg)
    elif (opt=="-q"):
        quiet = True

        
### ----- Set parameters ---------------------------------------------------

seed = 1651654  # probably want to change this for different runnumber's
                # - could have an array of seeds and select the runnumber-th one

Npoly=2000
Nsphereprot=2
Ndumbprot=20
Ndumbbeads=2*Ndumbprot
Nprot=Nsphereprot+Ndumbbeads


proteinONtype=2
proteinOFFtype=8

polymerunmarkedtype=1
polymerktype=6 #k27ac
polymerATACOFFtype=5
polymerATACONtype=7



totalruntime = 10000000              # units are timesteps
update_interval_steps = 15000        # this is how often switching is attempted
switch_interval_steps = 200000       # 1/k where k is switch rate in inverse steps

#totalruntime = 10000              # units are timesteps
#update_interval_steps = 15        # this is how often switching is attempted#
#switch_interval_steps = 200       # 1/k where k is switch rate in inverse steps

# need switch_interval_steps >> update_interval_steps to get a rate of ~k


### ----- Set variables based on parameters --------------------------------

Nupdates = int(totalruntime/update_interval_steps)
switch_rate = 1.0/switch_interval_steps                            
switch_prob = update_interval_steps*switch_rate 
oneminusswitch_prob = 1.0-switch_prob


### ----- Set up arrays for run --------------------------------------------

protState = numpy.zeros((Nprot, 1), dtype=int)  # an array to store the current state of proteins
polyState = numpy.zeros((Npoly, 1), dtype=int)  # an array to store the current state of polymer beads
# lammps requires the arrays to be of c-type
temp_proteinIDs = range(Npoly+1,Npoly+Nprot+1)  # maps the index (all sphere proteins) to lammps atom id
sphereprotIDs = range(Npoly+1,Npoly+Nsphereprot+1)
temp_polymerIDs = range(1,Npoly+1)

proteinIDs = (Nprot*ctypes.c_int)()
for i in range(Nprot):
    proteinIDs[i]=temp_proteinIDs[i]
del temp_proteinIDs

polymerIDs = (Npoly*ctypes.c_int)()
for i in range(Npoly):
    polymerIDs[i]=temp_polymerIDs[i]
del temp_polymerIDs


### ----- Set up random numbers --------------------------------------------

Myrandoms = numpy.random.RandomState(seed)


### ----- Write some messages to scree and/or log file ---------------------

if me==0 and not quiet:
    sys.stdout.write("#********************************************************************************\n")
    sys.stdout.write("#  Running simulation %i on %i processors\n"%(runnumber,nprocs))
    sys.stdout.write("#\n")
    sys.stdout.write("#  bridges switch on average every %f timesteps\n"%switch_interval_steps)
    sys.stdout.write("#\n")
    sys.stdout.write("# Probability: (try to keep this small)\n")
    sys.stdout.write("#  switch_prob = %f\n"%switch_prob)
    sys.stdout.write("#********************************************************************************\n\n")



### ----- Set up LAMMPS ----------------------------------------------------

# initialize LAMMPS 
lmp = lammps(cmdargs=["-screen","screen.%i.log"%runnumber,"-log","log.%i.log"%runnumber])

# pass runnumber and seed to lammps
lmp.command("variable runnumber equal %i"%runnumber)
lmp.command("variable seed equal %i"%seed)

# load lammps script for setup
lmp.file("in.setup")

# get number of atoms
natoms = lmp.get_natoms()

# get box size (assume square)
lx = lmp.extract_global("boxxhi",1) - lmp.extract_global("boxxlo",1)
hl = 0.5*lx

### ----- Do some checks ---------------------------------------------------

if  natoms != Nprot+Npoly :
    sys.exit("Error, atom numbers do not match ics\n")
if nprocs>1:
    sys.stdout.write("Error, python script not set up for parallel\n")
    sys.exit(1)


### ----- Set the initial types of the proteins ----------------------------

# switch half on, half off to start
mid=numpy.int(len(sphereprotIDs)/2)
first=[ str(id) for id in sphereprotIDs[0:mid] ]
first=" ".join(first)   # make a big string with a list of protein ids
second=[ str(id) for id in sphereprotIDs[mid:] ]
second=" ".join(second)

protState[0:mid]=1  # ON
protState[mid:Nsphereprot]=0   # OFF
RWID=proteinIDs[Nsphereprot:]
readerIDs=numpy.zeros(int(Ndumbprot))
writerIDs=numpy.zeros(int(Ndumbprot))


a=0
b=0
# Create a NumPy array to store the even integers from RWID
even_indices = numpy.zeros(int(Ndumbprot))
odd_indices = numpy.zeros(int(Ndumbprot))

# Get the ID for readers and writers based on odd or even, be aware that this only works when Npoly+Nsphere is even
# (which is normally the case) 

for i in (RWID):
    if i%2 ==0:
        even_indices[a]=int(i)-Npoly-1
        a += 1
    else:
        odd_indices[b]=int(i)-Npoly-1
        b += 1
readerIDs=even_indices.astype(int)
writerIDs=odd_indices.astype(int)




protState[readerIDs]=3   # RW pair which won't switch
protState[writerIDs]=4 
#print(protState)
#exit()


# send commands to lammps
lmp.command("group toON id "+first)
lmp.command("set group toON type %i"%proteinONtype)
lmp.command("group toON delete")

lmp.command("group toOFF id "+second)
lmp.command("set group toOFF type %i"%proteinOFFtype)
lmp.command("group toOFF delete")


### ----- Set the initial types of the polymers ----------------------------

# get one ATACOFF bead after every 99 unmarked beads
    
indices = [i for i in range(len(polymerIDs)) if (i+1) % 100 == 0 and i != 0]

# make an array of beads IDs then track all non-ATAC beads

polyID_array=numpy.zeros(Npoly)
for a in range(Npoly):
    polyID_array[a]=a
# get an array of all index apart from indices choosen for the ATACOFF beads
unmark = [element for element in polyID_array if element not in indices]
unmark = list(map(int, unmark))
# update the polymerstate, 1 for ATACOFF and 0 for unmarked
for i in range(len(polymerIDs)):
    polyState[i] = 1 if i in indices else 0
third=[ str(polymerIDs[id]) for id in indices ]
third=" ".join(third)   # make a big string with a list of polymer beads ids
fourth=[ str(polymerIDs[id]) for id in unmark]
fourth=" ".join(fourth)

#print(polyState)

# send commands to lammps
lmp.command("group toATACOFF id "+third)
lmp.command("set group toATACOFF type %i"%polymerATACOFFtype)
lmp.command("group toATACOFF delete")

lmp.command("group unmarked id "+fourth)
lmp.command("set group unmarked type %i"%polymerunmarkedtype)
lmp.command("group unmarked delete")



def calculate_distance(array1, array2):
    """ Calculate Euclidean distance between two points in 3D """
    if abs(array1[0]-array2[0])>hl:
        dx=2*hl-abs(array1[0]-array2[0])
    else:
        dx=array1[0]-array2[0]
    if abs(array1[1]-array2[1])>hl:
        dy=2*hl-abs(array1[1]-array2[1])
    else:
        dy=array1[1]-array2[1]
    if abs(array1[2]-array2[2])>hl:
        dz=2*hl-abs(array1[2]-array2[2])
    else:
        dz=array1[2]-array2[2]
    return numpy.sqrt(dx**2+dy**2+dz**2)

def decide(candidates,output):
    with open('chosen_polymers.txt', 'a') as file:
        for i in (numpy.array(polymerIDs)[candidates]):         # here i is the ID, which is 1 more than the index
        #print(i)
            polymer = polyPositions_reshape[i-1]
            # iterate over polymer beads
            for j in range(len(writerpositions)):
                writer = writerpositions[j]
                # calculate distance
                distance = calculate_distance(polymer, writer)
                #line = 'ID: {}, writer:{},distance: {}\n'.format(i, j+1,distance)
                #file.write(line)
                # Record changing for beads that are within a distance of 5 from a writer
                if distance <5:
                    # with a probability of 0.5
                    if random.random() < 0.5:
                        # record the polymer id
                        output=numpy.append(output,int(i-1))
                        line = 'ID: {},writer:{}, distance: {}\n'.format(i,2004+2*j, distance)
                        file.write(line)
                    break
                
    return output
### ----- Do the run --------------------------------------------------------

if me==0 and not quiet:
    sys.stdout.write("Starting the main run\n")
    sys.stdout.flush()

    
# run 
# Clear the testing file
#with open ('writeme.txt', 'w') as file: 
#    pass
#with open('distancetest.txt','w') as file1:
#    pass
with open('chosen_polymers.txt','w') as file:
    pass
for step in range(Nupdates):
    if (me == 0 and step%1==0) and not quiet:   # output some info now and again
        sys.stdout.write("running %i   ( %.2f %% done )     \r"%(step*update_interval_steps,100.0*(step)/float(Nupdates)))
        sys.stdout.flush()
    # count the number of lammps commands issued in this round
    flag_didAcommand = 0
    # get atom positions of proteins
    protPositions = lmp.gather_atoms_subset("x",1,3,Nprot,proteinIDs)  # 1 means floats; 3 is the number of floats per atom
    polyPositions = lmp.gather_atoms_subset("x",1,3,Npoly,polymerIDs)
    protPositions_np = numpy.array(protPositions)
    polyPositions_np = numpy.array(polyPositions)

    

    # open the file in append mode and write the timestep and atoms positions, these are for testing and debugging propose and
    # can be commented out in actual run
    with open ('chosen_polymers.txt', 'a') as file:  
         file.write(f"\nTimestep: {step}\n")
        

    polyPositions_reshape = polyPositions_np.reshape(-1,3)
    protPositions_reshape = protPositions_np.reshape(-1,3)
    writerpositions = protPositions_reshape[writerIDs]

    # choose which proteins to switch
    nowOn = numpy.where(protState==1)[0]
    nowOff = numpy.where(protState==0)[0]   
    # choose which polymers to switch, number here stands for bead type, not the actual all atom types in README file
    nowATACOFF = numpy.where(polyState==1)[0]
    nowunmark = numpy.where(polyState==0)[0]    
    nowATACON = numpy.where(polyState==2)[0]
    nowk27 = numpy.where(polyState==3)[0]

    # empty array to store the ID of beads that needs to be changed
    choosen_ATAC = numpy.array([])
    choosen_unmark = numpy.array([])


    choosen_ATAC=decide(nowATACOFF,choosen_ATAC)
    #print('chosen ATAC:',choosen_ATAC)
    choosen_unmark=decide(nowunmark,choosen_unmark)

    choosen_ATAC=choosen_ATAC.astype(int)
    choosen_unmark=choosen_unmark.astype(int)

    
    if me==0:
        #print("Choosing on->off")
        on2off = nowOn[ numpy.where(Myrandoms.choice([0, 1], size=len(nowOn), p=[oneminusswitch_prob, switch_prob])==1) ]
        #print("Choosing off->on")
        off2on = nowOff[ numpy.where(Myrandoms.choice([0, 1], size=len(nowOff), p=[oneminusswitch_prob, switch_prob])==1) ]
        # Myrandoms.choice() generates an array with 1 or 0 randomly with the given probabilities
        # numpy.where( ...==1 ) selects out the ones where a 1 was chosen
        # so on2off is a list of on proteins we are going to switch off etc.
        
        # get the id list as a string (probably there is a more efficient way)
        # and update the states array
        on2offgroup = ""
        for i in on2off:     
            on2offgroup = on2offgroup + " " + str(proteinIDs[i])
            protState[i] = 0
            
        off2ongroup = ""
        for i in off2on:
            off2ongroup = off2ongroup + " " + str(proteinIDs[i])
            protState[i] = 1

    # if running in parallel, uncomment to share this info between processors
    #comm.Bcast(protState, root=0) # fixed length numpy array so can use fast MPI Bcast version
    #on2offgroup = comm.bcast(on2offgroup, root=0) # to send a string need, to use slow python pickle version of bcast
    #off2ongroup = comm.bcast(off2ongroup, root=0)

    # Send the commands to lammps to do the switching
    
    # ON->OFF switch
    if on2offgroup != "":
        #print("on2off %s"%on2offgroup)
        lmp.command("group swON2OFF id "+on2offgroup)
        lmp.command("set group swON2OFF type %i"%proteinOFFtype)
        lmp.command("group swON2OFF delete")
        flag_didAcommand += 1
        
    # OFF->ON switch
    if off2ongroup != "":
        #print("off2on %s"%off2ongroup)
        lmp.command("group swOFF2ON id "+off2ongroup)
        lmp.command("set group swOFF2ON type %i"%proteinONtype)
        lmp.command("group swOFF2ON delete")
        flag_didAcommand += 1

    # For testing, change this if to a true statement
    if 1==2 :
        actualStates = lmp.gather_atoms_subset("type",0,1,Nprot,proteinIDs) # 0 means ints; 1 is the number of ints per atom
        actualIDs = lmp.gather_atoms_subset("id",0,1,Nprot,proteinIDs)
        for i in range( len(actualStates) ):            
            if protState[i]==0:
                atype = proteinOFFtype
            else:
                atype = proteinONtype
            if atype != actualStates[i]:
                sys.stdout.write("Error, lammps and python types don't match %i: %i,%i : %i %i\n"%(i,actualIDs[i],proteinIDs[i],atype,actualStates[i]))
                sys.stdout.flush()
                sys.exit(1)
        


    #print(polymerIDs[:])
    
    
    if me==0:
        ATACoff2on = choosen_ATAC
        unmark2k = choosen_unmark
        ATACon2off = nowATACON[ numpy.where(Myrandoms.choice([0, 1], size=len(nowATACON), p=[oneminusswitch_prob, switch_prob])==1) ]
        k2unmark = nowk27[ numpy.where(Myrandoms.choice([0, 1], size=len(nowk27), p=[oneminusswitch_prob, switch_prob])==1) ]

        ATACoff2ongroup = ""
        for i in ATACoff2on:     
            ATACoff2ongroup = ATACoff2ongroup + " " + str(polymerIDs[i])
            polyState[i] = 2
            #print(ATACoff2ongroup)   
        unmark2kgroup = ""
        for i in unmark2k:
            unmark2kgroup = unmark2kgroup + " " + str(polymerIDs[i])
            polyState[i] = 3
            #print(unmark2kgroup)
            

        ATACon2offgroup = ""
        for i in ATACon2off:     
            ATACon2offgroup = ATACon2offgroup + " " + str(polymerIDs[i])
            polyState[i] = 1
            #print(ATACon2offgroup)
        k2unmarkgroup = ""
        for i in k2unmark:
            k2unmarkgroup = k2unmarkgroup + " " + str(polymerIDs[i])
            polyState[i] = 0
    #print(polyState)
    # if running in parallel, uncomment to share this info between processors
    #comm.Bcast(protState, root=0) # fixed length numpy array so can use fast MPI Bcast version
    #on2offgroup = comm.bcast(on2offgroup, root=0) # to send a string need, to use slow python pickle version of bcast
    #off2ongroup = comm.bcast(off2ongroup, root=0)

    # Send the commands to lammps to do the switching
    
    # ATACOFF->ATACON switch
    if ATACoff2ongroup != "":
        #print("ATACoff2on %s"%ATACoff2ongroup)
        lmp.command("group swATACOFF2ON id "+ATACoff2ongroup)
        lmp.command("set group swATACOFF2ON type %i"%polymerATACONtype)
        lmp.command("group swATACOFF2ON delete")
        flag_didAcommand += 1
        
    # unmark->k27ac switch
    if unmark2kgroup != "":
        #print('aaa')
        #print("unmark2k %s"%unmark2kgroup)
        lmp.command("group swunmark2k id "+unmark2kgroup)
        lmp.command("set group swunmark2k type %i"%polymerktype)
        lmp.command("group swunmark2k delete")
        flag_didAcommand += 1

    # ATACON->ATACOFF switch
    if ATACon2offgroup != "":
        #print("ATACon2off %s"%ATACon2offgroup)
        lmp.command("group swATACON2OFF id "+ATACon2offgroup)
        lmp.command("set group swATACON2OFF type %i"%polymerATACOFFtype)
        lmp.command("group swATACON2OFF delete")
        flag_didAcommand += 1
        
    # k27ac-> unmark switch
    if k2unmarkgroup != "":
        #print("k2unmark %s"%k2unmarkgroup)
        lmp.command("group swk2unmark id "+k2unmarkgroup)
        lmp.command("set group swk2unmark type %i"%polymerunmarkedtype)
        lmp.command("group swk2unmark delete")
        flag_didAcommand += 1
  
        
    # do the next lammps run   
    if flag_didAcommand == 0:    
        # if nothing has changed, can skip 'pre' setup before run          
        lmp.command("run %i pre no post no"%update_interval_steps)     
    else: # in both cases do set 'post no' to skip calculating timing stats     
        lmp.command("run %i post no"%update_interval_steps)


# Write a restart data file
lmp.command("write_data DNA.${runnumber}.restart nocoeff")

# We're done!
sys.stdout.write("Done")
