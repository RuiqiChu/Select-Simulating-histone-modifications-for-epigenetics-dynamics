import numpy as np
import matplotlib.pyplot as plt

# Load the data from the text file
data = np.loadtxt("sphere_percentage_data.txt")

# Extract frames and sphere_percentage
frames = data[:, 0]/20000
#sphere_percentage = data[:, 1]
reader_percentage = data[:, 1]

average_reader_percentage = np.mean(reader_percentage)

# Plot the data
plt.plot(frames, reader_percentage, label="Reader")
#plt.plot(frames, sphere_percentage, label="Sphere")
plt.xlabel("Frame Number")
plt.ylabel("Fraction of Bonded Proteins")
plt.title("Fraction of Bonded Proteins vs. Time(20 pairs)")
# Add legend for protein types
first_legend = plt.legend(loc="upper right", title="Protein Type")

ax = plt.gca().add_artist(first_legend)
# Add another legend for the average value
plt.legend([f"Average Sphere Percentage: {average_reader_percentage:.3f}"], loc="upper left")
plt.savefig('0s20d12morse8!2-5lj_test')
# Set y-axis limits to range from 0 to 1
plt.ylim(0, 1)

# Show the plot
plt.show()
