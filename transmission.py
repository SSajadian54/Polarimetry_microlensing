import speclite.filters
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
from matplotlib import rcParams
from numpy import ma
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
import math
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib as mpl

from pylab import figure, show, legend, ylabel




plt.clf()


bessell = speclite.filters.load_filters('bessell-*')
speclite.filters.plot_filters(bessell, wavelength_limits=(2900, 9300))


print bessell[0]


print "bessel is : ",  rband


fig=plt.gcf()
fig.savefig("hotstar.jpg",dpi=200)
print(">>>>>>>>>>>>> The model light curve was made <<<<<<<<<<<<<<<")



