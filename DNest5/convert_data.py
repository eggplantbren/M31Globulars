import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load original dataset (CSV)
original = pd.read_csv("M31.1.csv")

# Remove ones with metallicity >= -0.4
data = original.loc[original["Fe_H"] < -0.4, :]

# Extract the six columns
cols = [data["x"], data["y"], data["Fe_H"], data["error"],
        data["V_M31"], data["V_sig"]]
data = np.vstack(cols).T

np.savetxt("data.txt", data)
plt.plot(data[:,0], data[:,1], "o")
plt.show()
