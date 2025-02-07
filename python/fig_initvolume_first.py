import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from proPlot import *

ApplyStyle()


tr = transforms.Affine2D().translate(-750,-750)
tr = tr.rotate_deg(180)

plt.figure(figsize=figSize(0.6,0.6))

figValues = plt.imread(f"comparaison_with_without_init_volume/movie_with/out1.png")

plt.imshow(figValues, transform=tr + plt.gca().transData, interpolation='none')
plt.xlim(55, 695)
plt.ylim(55, 695)
plt.xticks([])
plt.yticks([])


saveFig(plt.gcf(), "initvolume_show")