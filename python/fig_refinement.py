import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from proPlot import *

ApplyStyle()

fig_lvl7 = plt.imread("output3D/MD/movie_md_lvl7.png")
fig_lvl8 = plt.imread("output3D/noMD/lvl8_noMD.png")

plt.figure()

plt.subplot(1,2,1)
plt.imshow(np.flipud(fig_lvl7), interpolation='none')
plt.xlim(20, 730)
plt.ylim(20, 730)
plt.xticks([])
plt.yticks([])
plt.xlabel(r"a) \texttt{level} $= 7$")

plt.subplot(1,2,2)
plt.imshow(np.flipud(fig_lvl8), interpolation='none')
plt.xlim(20, 730)
plt.ylim(20, 730)
plt.xticks([])
plt.yticks([])
plt.xlabel(r"a) \texttt{level} $= 8$")


saveFig(plt.gcf(), "comp_ref")