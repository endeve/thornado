import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({'font.size': 13})#,'family':'monospace'})
mpl.rcParams.update({'legend.fontsize': 11})

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

mpl.rcParams['pgf.texsystem'] = 'pdflatex'
mpl.rcParams.update({'pgf.rcfonts' : False})

mpl.rcParams['lines.linewidth'] = 3.25

mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1

mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 1

mpl.rcParams['xtick.major.pad']='8'
mpl.rcParams['ytick.major.pad']='8'
colors = ["#430067", "#94216a", "#ff004d", "#ff8426", "#ffdd34", "#50e112", "#3fa66f", "#365987", "#000000", "#0033ff", "#29adff", "#00ffcc", "#c2c3c7", "#ab5236", "#5f574f"]
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)
