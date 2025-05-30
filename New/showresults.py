import dnest4.classic as dn4
import corner
import matplotlib.pyplot as plt
import pandas as pd
import os


dn4.postprocess()

os.system("R CMD BATCH convert.R")


posterior_samples = pd.read_csv("postdf.csv")
corner.corner(posterior_samples,
              plot_density=False, plot_contours=False, fontsize=14,
              hist_kwargs={"color":"blue", "alpha":0.3, "histtype":"stepfilled",
                           "edgecolor":"black","lw":"3"})

plt.show()

