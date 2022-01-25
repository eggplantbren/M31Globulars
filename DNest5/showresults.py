import matplotlib.pyplot as plt
import sys
import yaml

f = open("config.yaml")
dnest5_dir = yaml.load(f, Loader=yaml.SafeLoader)["dnest5_dir"]
f.close()

if __name__ == "__main__":
    sys.path.insert(0, dnest5_dir)
    import showresults as sr
    sr.standard_results(sys.argv)

    import numpy as np
    import pandas as pd
    samples = pd.read_csv("output/posterior.csv")
    data = np.loadtxt("data.txt")

    colors = ["red" if v < 0.0 else "blue" for v in data[:,3]]
    plt.scatter(data[:,0], data[:,1], s=30*np.abs(data[:,2]), marker="o",
                alpha=0.3, color=colors)
    plt.xlim([-5, 5])
    plt.ylim([-5, 5])
    plt.axis("square")

    phi1 = np.array(samples["phi1"])
    phi2 = np.array(samples["phi2"])
    phi1 = phi1[samples["M_crit"] < -1.5]
    phi2 = phi2[samples["M_crit"] < -1.5]

    xs = np.linspace(-5, 5, 11)
    for i in range(min(1000, len(phi1))):
        ys = np.tan(phi1[i])*xs
        plt.plot(xs, ys, "m-", linewidth=1, alpha=0.01)
        ys = np.tan(phi2[i])*xs
        plt.plot(xs, ys, "c-", linewidth=1, alpha=0.01)


    plt.show()

