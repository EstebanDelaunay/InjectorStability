import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from proPlot import *

ApplyStyle()

fig_lvl8_noMD = plt.imread("output3D/noMD/lvl8_noMD.png")
fig_lvl8 = plt.imread("python/images/movie_case1_0.150.png")

plt.figure()

plt.subplot(1,2,1)
plt.imshow(np.flipud(fig_lvl8_noMD), interpolation='none')
plt.xlim(20, 730)
plt.ylim(20, 730)
plt.xticks([])
plt.yticks([])
plt.xlabel("a) Without Manifold Death")

plt.subplot(1,2,2)
plt.imshow(fig_lvl8, interpolation='none')
plt.xlim(20, 730)
plt.ylim(20, 730)
plt.xticks([])
plt.yticks([])
plt.xlabel("b) With Manifold Death")


saveFig(plt.gcf(), "comp_md")