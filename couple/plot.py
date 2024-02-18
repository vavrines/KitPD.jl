import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy import interpolate
# read data
na = 'Damage'
file1 = f"{na}.txt"
a = np.loadtxt(file1)
Coordx = a[:, 0]
Coordy = a[:, 1]
Dmg = a[:, 2]
minX = -0.05
maxX = 0.05
minY = -0.05
maxY = 0.05
X = np.linspace(minX, maxX, 100)
Y = np.linspace(minY, maxY, 100)
X1, Y1 = np.meshgrid(X, Y)
D = interpolate.griddata((Coordx, Coordy), Dmg, (X1, Y1), method='linear')
fig, ax = plt.subplots(figsize=(12, 12))
levels = np.arange(0.0, 1, 0.01)
cset1 = ax.contourf(X1, Y1, D, levels, cmap=cm.jet)
ax.set_title("Damage", size=20)
ax.set_xlim(minX, maxX)
ax.set_ylim(minY, maxY)
ax.set_xlabel("X(m)", size=15)
ax.set_ylabel("Y(m)", size=15)
cbar = fig.colorbar(cset1)
cbar.set_label('Damage', size=18)
fig.savefig(f"{na}.jpg", bbox_inches='tight', dpi=300, pad_inches=0.0)
plt.clf()
