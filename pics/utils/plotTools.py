#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 13:38:07 2019

@author: oshah
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib.ticker import (
    MultipleLocator,
    FormatStrFormatter,
    AutoMinorLocator,
    LogLocator,
    FixedLocator,
)
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from matplotlib.colors import LinearSegmentedColormap

from pics.runparams import plot_params


colors = [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1]]


def reverse_colourmap(cmap, name="my_cmap_r"):
    """
    In:
    cmap, name
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1 - t[0], t[2], t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k, reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)

    return my_cmap_r


class ColorMap:
    def __init__(self, col=1, row=1, n_bin=3, hspace=1.0, wspace=0.25):

        self.row = row
        self.col = col
        self.hspace = hspace
        self.wspace = wspace
        self.colors = colors
        self.n_bin = n_bin

        # initiate subplots
        self.fig, self.ax = plt.subplots(self.col, self.row)

        # because matplotlib uses a stupid and useless format to store the
        # axis, create here a 2d array that can be accessed consistently
        # regardless of the row x col format.
        if col == 1 and row == 1:
            self.ax = [[self.ax]]

        elif col == 1 and row > 1:
            self.ax = [[axx] for axx in self.ax]

        elif col > 1 and row == 1:
            self.ax = [self.ax]

        self.cmap_name = "my_list"

        # Create customized color list
        self.colors = []

        for n in range(n_bin):
            self.colors.append(
                (1 - 1 / (n + 1), 1 - 1 / (2 * (n + 1)), 1 - 1 / (n + 1))
            )

        # Create the colormap
        self.cm = LinearSegmentedColormap.from_list(
            self.cmap_name, self.colors, N=self.n_bin
        )
        # Fewer bins will result in "coarser" colomap interpolation

        N = 256

        borders = np.linspace(0, 1, n_bin + 1)
        refine = 5

        a = np.linspace(0, 0.01, refine, endpoint=False)
        b = np.linspace(0.01, 1, n_bin - refine)

        borders = np.concatenate((a, b))

        # cols = np.linspace(0, 1, n_bin)

        vals = np.ones((N, 4))
        vals[:, 0] = np.logspace(0, -2, N)  # red
        vals[:, 1] = np.linspace(0, 0, N)  # green
        vals[:, 2] = np.linspace(0, 0, N)  # blue

        # iterate over all pixels in the colormap
        for i in range(N):
            # iterate over all color segments in the colormap and decide
            # which segment the current pixel belongs to
            for j in range(len(borders) - 1):
                if i / N >= borders[j] and i / N < borders[j + 1]:
                    vals[i, 0] = colors[0][j]
                    vals[i, 1] = colors[1][j]
                    vals[i, 2] = colors[2][j]

            if i / N > borders[-2]:
                vals[i, 0] = colors[0][-1]
                vals[i, 1] = colors[1][-1]
                vals[i, 2] = colors[2][-1]

        self.cm = ListedColormap(vals)

        # self.fig.subplots_adjust(left=0.02, bottom=0.06, right=0.75, top=0.94, wspace=0.0)
        self.fig.subplots_adjust(hspace=self.hspace, wspace=self.wspace)
        """
        x = np.arange(0, np.pi, 0.1)
        y = np.arange(0, 2 * np.pi, 0.1)
        X, Y = np.meshgrid(x, y)
        Z = np.cos(X) * np.sin(Y) * 10
        
        for i in range(self.row):
            for j in range(self.col):
                ax = self.ax[i][j]
                # Create the colormap
                cm = LinearSegmentedColormap.from_list(
                    cmap_name, colors, N=n_bin)
                # Fewer bins will result in "coarser" colomap interpolation
                im = ax.imshow(Z, interpolation='nearest', origin='lower', cmap=cm,
                               aspect = 'auto')
                ax.set_title("N bins: %s" % n_bin)
                self.fig.colorbar(im, ax=ax)
                
        """

    def Plot(self, Z=None):
        for i in range(self.row):
            for j in range(self.col):
                ax = self.ax[i][j]
                self.im = ax.imshow(
                    Z,
                    interpolation="nearest",
                    origin="lower",
                    cmap=self.cm,
                    aspect="auto",
                )

                ax.set_title("N bins: %s" % self.n_bin)
                self.fig.colorbar(self.im, ax=ax)


class Plot:
    def __init__(
        self,
        col=1,
        row=1,
        axis_labels=None,
        axis_limits=None,
        minorlocatorx=None,
        majorlocatorx=None,
        minorlocatory=None,
        majorlocatory=None,
        x_minor_labels=None,
        y_minor_labels=None,
        x_major_labels=None,
        y_major_labels=None,
        hspace=1.0,
        wspace=0.25,
        plot_titles=None,
        logx=False,
        logy=False,
        majorlabelsize=14,
        minorlabelsize=8,
        axislabelsize=16,
        axislabelpad=5,
        majortickwidth=1,
        minortickwidth=1,
        minor_labels=False,
        sharex=False,
        sharey=False,
        titlefontsize=16,
        title_pad=10,
        figsize=(20, 10),
        grid=True,
        top=True,
        right=True,
        left=True,
        bottom=True,
        labelright=True,
        labeltop=True,
        labelleft=True,
    ):

        letters = ["a", "b", "c", "d", "e", "f", "g", "h"]

        self.sharex = sharex
        self.sharey = sharey
        self.logx = logx
        self.logy = logy
        self.hspace = hspace
        self.wspace = wspace
        self.col = col
        self.row = row

        self.plot_params = {
            "top": top,
            "bottom": bottom,
            "right": right,
            "left": left,
            "direction": "in",
            "pad": plot_params["pad"],
            "labelright": labelright,
            "labelleft": labelleft,
            "labeltop": labeltop,
        }

        # If no logx and logy are given, set them to False for all axes
        if logx == False:
            self.logx = [[False for j in range(col)] for i in range(row)]

        if logy == False:
            self.logy = [[False for j in range(col)] for i in range(row)]

        # If logx is set to True for all axes
        if logx == True:
            self.logx = [[True for j in range(col)] for i in range(row)]

        if logy == True:
            self.logy = [[True for j in range(col)] for i in range(row)]

        if not plot_titles == None:
            self.plot_titles = plot_titles

        else:
            self.plot_titles = [["" for j in range(col)] for i in range(row)]

        if not axis_limits == None:
            self.axis_limits = axis_limits

        else:
            self.axis_limits = [
                [[[0.0, 1.0], [0.0, 1.0]] for j in range(col)] for i in range(row)
            ]

        if majorlocatorx == None:
            majorlocatorx = [[None for j in range(col)] for i in range(row)]

        if majorlocatory == None:
            majorlocatory = [[None for j in range(col)] for i in range(row)]

        self.majorlocatorx = []
        self.majorlocatory = []
        for i in range(row):
            self.majorlocatorx.append([])
            self.majorlocatory.append([])
            for j in range(col):
                delta_lims_x = (
                    self.axis_limits[i][j][0][1] - self.axis_limits[i][j][0][0]
                )

                n_major_ticks_x = int(delta_lims_x) + 1

                delta_lims_y = (
                    self.axis_limits[i][j][1][1] - self.axis_limits[i][j][1][0]
                )

                n_major_ticks_y = int(delta_lims_y) + 1

                if majorlocatorx[i][j] == None or self.logx[i][j]:
                    if self.logx[i][j]:

                        self.majorlocatorx[i].append(
                            np.linspace(
                                int(self.axis_limits[i][j][0][0]),
                                int(self.axis_limits[i][j][0][1]),
                                n_major_ticks_x,
                            )
                        )

                    else:
                        self.majorlocatorx[i].append(delta_lims_x / n_major_ticks_x)

                else:
                    self.majorlocatorx[i].append(majorlocatorx[i][j])

                if majorlocatory[i][j] == None or self.logy[i][j]:
                    if self.logy[i][j]:

                        self.majorlocatory[i].append(
                            np.linspace(
                                int(self.axis_limits[i][j][1][0]),
                                int(self.axis_limits[i][j][1][1]),
                                n_major_ticks_y,
                            )
                        )

                    else:
                        self.majorlocatory[i].append(delta_lims_y / n_major_ticks_y)

                else:
                    self.majorlocatory[i].append(majorlocatory[i][j])

        if not minorlocatorx == None:
            self.minorlocatorx = minorlocatorx

        else:
            self.minorlocatorx = [[0.5 for j in range(col)] for i in range(row)]

        if not minorlocatory == None:
            self.minorlocatory = minorlocatory

        else:
            self.minorlocatory = [[0.5 for j in range(col)] for i in range(row)]

        if not axis_labels == None:
            self.axis_labels = axis_labels

        else:
            self.axis_labels = [
                [["x [a.u.]", "y [a. u.]"] for j in range(col)] for i in range(row)
            ]

        self.fig, self.ax = plt.subplots(row, col, figsize=figsize)
        self.fig.subplots_adjust(hspace=self.hspace, wspace=self.wspace)

        # because matplotlib uses a stupid and useless format to store the
        # axis, create here a 2d array that can be accessed consistently
        # regardless of the row x col format.
        if col == 1 and row == 1:
            self.ax = [[self.ax]]

        elif col == 1 and row > 1:
            self.ax = [[axx] for axx in self.ax]

        elif col > 1 and row == 1:
            self.ax = [self.ax]

        self.create_ticklabels()

        for i in range(row):
            for j in range(col):
                axx = self.ax[i][j]
                xmin = self.axis_limits[i][j][0][0]
                xmax = self.axis_limits[i][j][0][1]
                ymin = self.axis_limits[i][j][1][0]
                ymax = self.axis_limits[i][j][1][1]

                # set width of frame surrounding the plots
                for pos in ["top", "bottom", "left", "right"]:
                    axx.spines[pos].set_linewidth(1.0)

                # set labels to right for 2 column plots
                if col == 2 and j == 1:
                    axx.yaxis.tick_right()
                    axx.yaxis.set_label_position("right")

                axx.set_xlim(self.axis_limits[i][j][0])
                axx.set_ylim(self.axis_limits[i][j][1])

                # Only set title for upper row
                if i == 0:
                    axx.set_title(
                        self.plot_titles[i][j], fontsize=titlefontsize, pad=title_pad
                    )

                # Set x-axis labels
                axx.set_xlabel(
                    self.axis_labels[i][j][0],
                    fontsize=axislabelsize,
                    labelpad=axislabelpad,
                )

                # Set y-axis labels
                axx.set_ylabel(
                    self.axis_labels[i][j][1],
                    fontsize=axislabelsize,
                    labelpad=axislabelpad,
                )

                axx.tick_params(which="both", **self.plot_params, zorder=1000)

                # If shared x-axis, turn off top tick labels for lower panels and
                # axis labels for all upper panels
                # and lower labels for all upper panels
                if self.sharex:
                    if i < row - 1:
                        axx.xaxis.label.set_visible(False)
                        axx.tick_params(which="both", labelbottom=False)

                    if i > 0:
                        axx.tick_params(which="both", labeltop=False)

                # If shared y-axis, turn off right tick labels for left panels
                # and left axis labels for all right panels
                # and right axis labels for all left panels
                if self.sharey:
                    if j > 0:
                        print("turn off labels")
                        axx.yaxis.label.set_visible(False)
                        axx.tick_params(which="both", labelleft=False)

                    if j < col - 1:
                        axx.tick_params(which="both", labelright=False)

                axx.tick_params(
                    which="major",
                    length=plot_params["majorticklen"],
                    labelsize=majorlabelsize,
                    width=majortickwidth,
                )
                axx.tick_params(
                    which="minor",
                    length=plot_params["minorticklen"],
                    labelsize=minorlabelsize,
                    width=minortickwidth,
                )

                plt.rcParams["axes.linewidth"] = 100

                axx.zorder = 10
                # create log ticks
                if self.logx[i][j]:
                    xticks_minor = []
                    for ll in range(int(xmax) - int(xmin) + 1):
                        left = xmin + ll
                        right = left + 1
                        # print('\n left/right =', left, right)
                        partial_xticks = np.linspace(10**left, 10**right, 10)
                        for lll in range(10):
                            xticks_minor.append(np.log10(partial_xticks[lll]))

                    xticks_major = self.majorlocatorx[i][j]

                    axx.xaxis.set_major_locator(FixedLocator(self.majorlocatorx[i][j]))
                    axx.xaxis.set_minor_locator(FixedLocator(xticks_minor))

                # create linear ticks
                else:
                    xticks_minor = np.arange(
                        xmin, xmax + self.minorlocatorx[i][j], self.minorlocatorx[i][j]
                    )

                    xticks_major = np.arange(
                        xmin, xmax + self.majorlocatorx[i][j], self.majorlocatorx[i][j]
                    )

                    axx.xaxis.set_major_locator(
                        MultipleLocator(self.majorlocatorx[i][j])
                    )

                    axx.xaxis.set_minor_locator(
                        MultipleLocator(self.minorlocatorx[i][j])
                    )

                # create log ticks
                if self.logy[i][j]:
                    yticks_minor = []
                    for ll in range(int(ymax) - int(ymin) + 1):
                        left = ymin + ll
                        right = left + 1
                        # print('\n left/right =', left, right)
                        partial_yticks = np.linspace(10**left, 10**right, 10)
                        for lll in range(10):
                            yticks_minor.append(np.log10(partial_yticks[lll]))

                    yticks_major = self.majorlocatory[i][j]

                    axx.yaxis.set_major_locator(FixedLocator(self.majorlocatory[i][j]))
                    axx.yaxis.set_minor_locator(FixedLocator(yticks_minor))

                # create linear ticks
                else:

                    yticks_minor = np.arange(
                        ymin, ymax + self.minorlocatory[i][j], self.minorlocatory[i][j]
                    )

                    yticks_major = np.arange(
                        ymin, ymax + self.majorlocatory[i][j], self.majorlocatory[i][j]
                    )

                    axx.yaxis.set_major_locator(
                        MultipleLocator(self.majorlocatory[i][j])
                    )

                    axx.yaxis.set_minor_locator(
                        MultipleLocator(self.minorlocatory[i][j])
                    )

                if minor_labels:
                    axx.set_xticklabels(self.x_minor_labels[i][j], minor=True)
                    axx.set_yticklabels(self.y_minor_labels[i][j], minor=True)

                axx.set_xticklabels(self.x_major_labels[i][j], minor=False)
                axx.set_yticklabels(self.y_major_labels[i][j], minor=False)

                axx.set_facecolor(plot_params["backgroundcol"])
                """
                axx.grid(axis='both', which='major', 
                         color=plot_params['gridcol'],
                         alpha=plot_params['gridalpha'],
                         linewidth=plot_params['lwdthgridmajor'],
                         zorder=0)

                axx.grid(axis='both', which='minor', 
                         color=plot_params['gridcol'],
                         alpha=plot_params['gridalpha'],
                         linewidth=plot_params['lwdthgridminor'],
                         zorder=0)
                """
                # construct grid manually because matplotlib builtin methods suck
                # and do not behave as desired

                if grid:
                    # draw minor x-grid lines
                    for x in xticks_minor:
                        axx.plot(
                            [x, x],
                            [ymin, ymax],
                            color=plot_params["gridcolor"],
                            linewidth=1,
                            zorder=-100,
                            alpha=plot_params["gridalpha"],
                        )

                    # draw minor y-grid lines
                    for y in yticks_minor:
                        axx.plot(
                            [xmin, xmax],
                            [y, y],
                            color=plot_params["gridcolor"],
                            linewidth=1,
                            zorder=-100,
                            alpha=plot_params["gridalpha"],
                        )

                    # draw major x-grid lines
                    for x in xticks_major:

                        axx.plot(
                            [x, x],
                            [ymin, ymax],
                            color=plot_params["gridcolor"],
                            linewidth=2,
                            zorder=-100,
                            alpha=plot_params["gridalpha"],
                        )

                    # draw major y-grid lines
                    for y in yticks_major:
                        axx.plot(
                            [xmin, xmax],
                            [y, y],
                            color=plot_params["gridcolor"],
                            linewidth=2,
                            zorder=-100,
                            alpha=plot_params["gridalpha"],
                        )

    def create_ticklabels(self):

        self.x_major_labels = []
        self.y_major_labels = []
        self.x_minor_labels = []
        self.y_minor_labels = []

        for i in range(self.row):
            self.x_major_labels.append([])
            self.y_major_labels.append([])
            self.x_minor_labels.append([])
            self.y_minor_labels.append([])

            for j in range(self.col):
                if self.logx[i][j]:
                    self.x_major_labels[i].append([])

                else:
                    self.x_major_labels[i].append([""])

                if self.logy[i][j]:
                    self.y_major_labels[i].append([])

                else:
                    self.y_major_labels[i].append([""])

                self.x_minor_labels[i].append([""])
                self.y_minor_labels[i].append([""])

                if self.logx[i][j]:

                    # process x axis
                    for ll in range(
                        int(self.axis_limits[i][j][0][1] - self.axis_limits[i][j][0][0])
                        + 1
                    ):
                        xlabel = self.axis_limits[i][j][0][0] + ll
                        self.x_major_labels[i][j].append(
                            r"$10^{{{1:d}}}$".format(1, xlabel)
                        )

                else:
                    delta_x = (
                        self.axis_limits[i][j][0][1] - self.axis_limits[i][j][0][0]
                    )

                    n_major_ticks_x = int(
                        round(delta_x / self.majorlocatorx[i][j] + 1, 0)
                    )

                    n_minor_ticks_x = int(
                        round(delta_x / self.minorlocatorx[i][j] + 1, 0)
                    )

                    x_major_labels_dummy = np.arange(
                        self.axis_limits[i][j][0][0],
                        self.axis_limits[i][j][0][1] + self.majorlocatorx[i][j],
                        self.majorlocatorx[i][j],
                    )

                    x_minor_labels_dummy = np.linspace(
                        self.axis_limits[i][j][0][0],
                        self.axis_limits[i][j][0][1],
                        n_minor_ticks_x,
                    )

                    # create customized major and minor tick labels
                    roundx = 3

                    # minor ticks
                    # generate x labels
                    for n in range(len(x_minor_labels_dummy)):
                        number = x_minor_labels_dummy[n]
                        if (
                            round(number / self.majorlocatorx[i][j], 0)
                            == round(number / self.majorlocatorx[i][j], roundx)
                            or n == len(x_minor_labels_dummy) - 1
                        ):
                            self.x_minor_labels[i][j].append("")

                        else:
                            self.x_minor_labels[i][j].append(str(round(number, roundx)))

                    # major ticks
                    # generate x labels
                    for n in range(len(x_major_labels_dummy)):
                        number = x_major_labels_dummy[n]
                        if abs((int(number) - number) / number) < 1.0e-6:
                            self.x_major_labels[i][j].append(str(int(number)))

                        else:
                            self.x_major_labels[i][j].append(str(round(number, roundx)))

                if self.logy[i][j]:
                    for ll in range(
                        int(self.axis_limits[i][j][1][1] - self.axis_limits[i][j][1][0])
                        + 1
                    ):
                        ylabel = self.axis_limits[i][j][1][0] + ll
                        self.y_major_labels[i][j].append(
                            r"$10^{{{1:d}}}$".format(1, ylabel)
                        )

                else:

                    delta_y = (
                        self.axis_limits[i][j][1][1] - self.axis_limits[i][j][1][0]
                    )

                    n_major_ticks_y = int(
                        round(delta_y / self.majorlocatory[i][j] + 1, 0)
                    )

                    n_minor_ticks_y = int(
                        round(delta_y / self.minorlocatory[i][j] + 1, 0)
                    )

                    y_major_labels_dummy = np.arange(
                        self.axis_limits[i][j][1][0],
                        self.axis_limits[i][j][1][1] + self.majorlocatory[i][j],
                        self.majorlocatory[i][j],
                    )

                    y_minor_labels_dummy = np.linspace(
                        self.axis_limits[i][j][1][0],
                        self.axis_limits[i][j][1][1],
                        n_minor_ticks_y,
                    )

                    # print ('minor locator y = ', self.minorlocatory[i][j])
                    # print ('delta y =', delta_y)
                    # print ('n minor ticks y =', n_minor_ticks_y)
                    # print ('x minor dummy =', x_minor_labels_dummy)
                    # print (y_minor_labels_dummy)

                    roundy = 3
                    # minor ticks
                    # generate y labels
                    for n in range(len(y_minor_labels_dummy)):
                        number = y_minor_labels_dummy[n]
                        if (
                            round(number / self.majorlocatory[i][j], 0)
                            == round(number / self.majorlocatory[i][j], roundy)
                            or n == len(y_minor_labels_dummy) - 1
                        ):
                            self.y_minor_labels[i][j].append("")

                        else:
                            self.y_minor_labels[i][j].append(str(round(number, roundy)))

                    # major ticks
                    # generate y labels
                    for n in range(len(y_major_labels_dummy)):
                        number = y_major_labels_dummy[n]
                        if abs((int(number) - number) / number) < 1.0e-6:
                            self.y_major_labels[i][j].append(str(int(number)))
                        else:
                            self.y_major_labels[i][j].append(str(round(number, roundy)))
        # print ('y minors =', self.y_minor_labels)
        # print ('x minors =', self.x_minor_labels)

        # print ('y majors =', self.y_major_labels)
        # print ('x majors =', self.x_major_labels)

    def Save(self, format="png", loc="cwd", file_name="fig.png", dpi=120):
        print("writing plot to file:", loc + file_name)
        self.fig.savefig(loc + file_name, dpi=dpi)
