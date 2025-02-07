import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from proPlot import *

ApplyStyle()

fig_f = plt.imread("python/images/multi_200.png")
fig_fc = plt.imread("python/images/adapt_200.png")

tr = transforms.Affine2D().translate(-750,-750)
tr = tr.rotate_deg(180)

plt.figure()

plt.subplot(1,2,1)
plt.imshow(fig_f, transform=tr + plt.gca().transData, interpolation='none')
plt.xlim(20, 730)
plt.ylim(20, 730)
plt.xticks([])
plt.yticks([])
plt.title("a)", loc="left")
plt.title("Fixed Grid", loc="center")

plt.subplot(1,2,2)
plt.imshow(fig_fc, transform=tr + plt.gca().transData, interpolation='none')
plt.xlim(20, 730)
plt.ylim(20, 730)
plt.xticks([])
plt.yticks([])
plt.title("b)", loc="left")
plt.title("Adaptive Grid", loc="center")

saveFig(plt.gcf(), "adaptmulti1")