import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from proPlot import *

ApplyStyle()

tPlot = [20, 25, 30]

tr = transforms.Affine2D().translate(-750,-750)
tr = tr.rotate_deg(180)

plt.figure()

# upper line
for i in range(3):
    figValues = plt.imread(f"comparaison_with_without_init_volume/movie_with/out{tPlot[i]}.png")

    plt.subplot(2,3,i+1)
    plt.imshow(figValues, transform=tr + plt.gca().transData, interpolation='none')
    plt.xlim(55, 695)
    plt.ylim(55, 695)
    plt.xticks([])
    plt.yticks([])
    
    if i==0: plt.ylabel("With initial volume")

for i in range(3):
    figValues = plt.imread(f"comparaison_with_without_init_volume/movie_without/out{tPlot[i]}.png")

    plt.subplot(2,3,i+4)
    plt.imshow(figValues, transform=tr + plt.gca().transData, interpolation='none')
    plt.xlim(55, 695)
    plt.ylim(55, 695)
    plt.xticks([])
    plt.yticks([])
    
    if i==0: plt.ylabel("Without initial volume")
    
    plt.xlabel(f"Time $t={tPlot[i]*5e-3}$")


saveFig(plt.gcf(), "initvolume")