import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("./hists/muons/final_coord.txt")

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]


fig, axs = plt.subplots(figsize=(12,12))
plt.hist2d(x=x, y=z, bins=[25,100], range=[[0, 100], [0, 104]], cmin=1e-24)
plt.xlabel("X")
plt.ylabel("Z")
plt.xlim([0,100])
plt.ylim([0,104])
fig.savefig("histXZ.pdf")

fig, axs = plt.subplots(figsize=(12,12))
plt.hist2d(x=y, y=z, bins=[25,100], range=[[0, 100], [0, 104]], cmin=1e-24)
plt.xlabel("Y")
plt.ylabel("Z")
plt.xlim([0,100])
plt.ylim([0,104])
fig.savefig("histYZ.pdf")
