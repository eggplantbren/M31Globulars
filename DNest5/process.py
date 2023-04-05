import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

inner = pd.read_csv("M31.1.csv")

x = inner["radius.kpc"]*np.cos((inner["theta"]+90)*np.pi/180.0)
y = inner["radius.kpc"]*np.sin((inner["theta"]+90)*np.pi/180.0)

x = np.array(x)
y = np.array(y)
v = np.array(inner["V_M31"])
sig_v = np.array(inner["V_sig"])
z = np.array(inner["Fe_H"])

# Now load outer data
outer = np.loadtxt("/home/brewer/Projects/AndromedaMixture/Data/data.txt")
x = np.hstack([x, outer[:,0]])
y = np.hstack([y, outer[:,1]])
v = np.hstack([v, outer[:,2]])
sig_v = np.hstack([sig_v, outer[:,3]])
z = np.hstack([z, outer[:,4]+1000])

data = np.vstack([x, y, v, sig_v, z]).T
np.savetxt("merged.txt", data)

away = v > 0.0
towards = v <= 0.0

# Square root trick for plot
phi = np.arctan2(y, x)
rr  = np.sqrt(x**2 + y**2)
rr = rr**1.0
xx = rr*np.cos(phi)
yy = rr*np.sin(phi)
plt.scatter(xx[away], yy[away], marker="o", color="r", alpha=0.5,
                s=0.2*np.abs(v[away]))
plt.scatter(xx[towards], yy[towards], marker="o", color="b", alpha=0.5,
                s=0.2*np.abs(v[towards]))
plt.axis("equal")
plt.show()

