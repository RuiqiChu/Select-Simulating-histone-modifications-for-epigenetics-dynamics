#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:06:47 2023

@author: s1969574
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

# # load the data into a pandas DataFrame
# df = pd.read_csv('atom_types.dat', delim_whitespace=True, comment='#', names=['timestep', 'id', 'type'])

# # pivot the data to create a matrix of types (values) at each timestep (rows) and atom id (columns)
# df_pivot = df.pivot(index='timestep', columns='id', values='type')

# # create a discrete colormap with specific colors for types 1 to 8
# colors = ['red', 'blue', 'green', 'yellow', 'purple', 'cyan', 'magenta', 'black']
# cmap = ListedColormap(colors) 

# # create the heatmap
# fig, ax = plt.subplots(figsize=(10, 8))
# sns.heatmap(df_pivot, cmap=cmap, cbar=False)

# # create a custom legend
# type_names = ['unmarked', 'TF ON', 'Reader', 'Writer', 'ATAC OFF', 'K27ac', 'ATAC ON', 'TF OFF']
# patches = [mpatches.Patch(color=color, label=type_name) for color, type_name in zip(colors, type_names)]
# plt.title('Type change heatmap(400s50d11-5morse8!2-5lj)')
# plt.legend(handles=patches, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)

# plt.savefig("400s50d11-5morse8!2-5lj.png")
# # show the plot
# plt.show()
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

# Load the data into a pandas DataFrame
df = pd.read_csv('atom_types.dat', delim_whitespace=True, comment='#', names=['timestep', 'id', 'type'])

# Pivot the data to create a matrix of types (values) at each timestep (rows) and atom id (columns)
df_pivot = df.pivot(index='timestep', columns='id', values='type')

# Define the types for ON and OFF states
on_types = [6, 7]
off_types = [1, 5]

# Calculate the fraction of ON state for each atom
on_fraction_per_atom = df_pivot.apply(lambda x: (x.isin(on_types).sum() / len(x)), axis=0)

# Create a discrete colormap with specific colors for types 1 to 8
colors = ['red', 'blue', 'green', 'yellow', 'purple', 'cyan', 'magenta', 'black']
cmap = ListedColormap(colors) 

# Create the heatmap
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(df_pivot, cmap=cmap, cbar=False)

# Create a custom legend
type_names = ['unmarked', 'TF ON', 'Reader', 'Writer', 'ATAC OFF', 'K27ac', 'ATAC ON', 'TF OFF']
patches = [mpatches.Patch(color=color, label=type_name) for color, type_name in zip(colors, type_names)]
plt.title('Type change heatmap(20 pairs)')
plt.legend(handles=patches, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
plt.savefig("0s20d12morse8!2-5lj")
# Plot the ON state fraction
plt.figure()
plt.plot(on_fraction_per_atom.index, on_fraction_per_atom.values, linestyle='-', color='blue')
plt.xlabel('Atom ID')
plt.ylabel('Fraction of ON State')
plt.title('Fraction of ON State for Each Atom(20 pairs)')
plt.xticks(range(1, 2001,100))
#plt.xticks(range(100,400,10))
plt.ylim(0, 1)
#plt.xlim(0,2100)
plt.savefig("ON_state_fraction(0s20d12morse8!2-5lj_test).png")

# Show the plots
plt.show()
