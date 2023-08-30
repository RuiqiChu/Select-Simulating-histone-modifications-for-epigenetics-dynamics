import pandas as pd
import numpy as np
import math

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
    return np.sqrt(dx**2+dy**2+dz**2)

def analyze_bonding(data, array2):
    # Convert the 'atom_id' column to integers
    data["atom_id"] = data["atom_id"].astype(int)
    unique_frames = data["frame"].unique()
    
    with open("bonded_atoms.txt", "w") as f:
        for i, frame in enumerate(unique_frames):
            frame_data = data[data["frame"] == frame]
            total_readers = len(frame_data[frame_data["atom_id"].between(2002, 2042, inclusive="both") & (frame_data["atom_id"] % 2 == 1)])
            bonded_readers = 0
            reader_proteins = frame_data[frame_data["atom_id"].between(2002, 2042, inclusive="both") & (frame_data["atom_id"] % 2 == 1)]
            polymer_beads = frame_data[frame_data["atom_id"] <= 2000]

            for _, reader in reader_proteins.iterrows():
                is_bonded = False
                for _, polymer in polymer_beads.iterrows():
                    distance = calculate_distance([reader[3], reader[4], reader[5]],
                                                  [polymer[3], polymer[4], polymer[5]])
                    if distance < 1.2:
                        is_bonded = True
                        # Record the bonded atom ID and timestep
                        f.write(f"Atom ID: {reader['atom_id']}, position:{[reader[3], reader[4], reader[5]]},Polymer ID: {polymer['atom_id']}, position:{[polymer[3], polymer[4], polymer[5]]},Timestep: {frame}\n")
                        break
                if is_bonded:
                    bonded_readers += 1
            reader_percentage_now = bonded_readers / total_readers
            array2[i] = reader_percentage_now

    return array2



# Read data from the file and skip the first line (header)
file_path = "atom_positions.txt"  # File is in the same directory
#print('data read')
column_names = ["frame", "atom_id", "atom_type", "x", "y", "z"]
atoms_data = pd.read_csv(file_path, sep=" ", header=None, names=column_names, skiprows=1)
# Calculate the fraction of bonded sphere proteins for each frame
hl = 50
# Plot the fraction of bonded sphere proteins versus time
frames = atoms_data["frame"].unique()
reader_percentage = np.zeros(len(frames))
reader_percentage = analyze_bonding(atoms_data,reader_percentage)

# Save the frames and sphere_percentage to a text file
with open("sphere_percentage_data.txt", "w") as data_file:
    for frame, percentage2 in zip(frames, reader_percentage):
        data_file.write(f"{frame} {percentage2}\n")
print('All done, bro')