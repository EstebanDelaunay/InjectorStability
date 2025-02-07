import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from proPlot import *

ApplyStyle()

fig_tol5 = plt.imread("output3D/comp_tol/tol5.png")
fig_tol7 = plt.imread("output3D/comp_tol/tol7.png")

plt.figure()

plt.subplot(1,2,1)
plt.imshow(np.flipud(fig_tol5), interpolation='none')
plt.xlim(20, 730)
plt.ylim(20, 730)
plt.xticks([])
plt.yticks([])
plt.xlabel(r"a) \texttt{TOLERANCE} $= 10^{-5}$")

plt.subplot(1,2,2)
plt.imshow(np.flipud(fig_tol7), interpolation='none')
plt.xlim(20, 730)
plt.ylim(20, 730)
plt.xticks([])
plt.yticks([])
plt.xlabel(r"a) \texttt{TOLERANCE} $= 10^{-7}$")

saveFig(plt.gcf(), "comp_tol")