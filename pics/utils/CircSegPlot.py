# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 11:54:09 2022

@author: oliver
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc
import matplotlib as mpl
import numpy as np
import matplotlib
import copy
import functionTools as ftool
from PIMPrunparams import layerColors, layerCodes
import matplotlib.transforms as transforms
from PIMPphysicalparams import m_earth

matplotlib.rcParams["text.usetex"] = True
# This mothafucka speeds up saving graphics with mpl like crazy
matplotlib.use("Agg")

orderedLayers = [
    "inner core",
    "outer core",
    "mantle",
    "solid",
    "supercritical",
    "liquid",
]


class PlanetOverviewChart:
    def __init__(self, structures):
        self.colors = layerColors
        self.structures = structures

        # Array containing the layer keys for all structures in the chart
        self.data = None

    def set_data(self):
        self.data = ["0", "1", "2", "02", "210"]

    def plot_advanced(self, showLegend=True, fnts=16, auto_shape=True, **kwargs):

        # Determines best fit layout for given number of structures
        if auto_shape:
            pass
        else:
            try:
                shape = kwargs["shape"]
            except KeyError:
                print(
                    "WARNING: For auto_shape = False the key argument 'shape' needs to be specified."
                )

    def plot_simple(self, showLegend=True, fnts=16):
        nCols = len(self.structures)

        fig, ax = plt.subplots(1, nCols, figsize=(10, 4))
        alphas = ["a", "b", "c", "d"]
        # Add pie slices for each structure
        for i in range(len(self.structures)):
            cols = [self.colors[st] for st in self.structures[i]]
            pc = PlanetChart(nLayers=len(self.structures[i]), colors=cols)
            pc.draw(ax=ax[i])
            ax[i].text(0.2, -1, alphas[i], fontsize=fnts)

        # Add layer color legend
        if showLegend:
            axx = ax[-1].inset_axes([1.0, 0, 1.0, 1.0])
            # Gather unique structure components
            uniqueStructures = list(set.union(*map(set, self.structures)))

            dy = 0.2
            for i in range(len(uniqueStructures)):
                yPos = 0.5 + dy * 2 - i * dy
                axx.text(0.3, yPos, orderedLayers[i], va="center", fontsize=fnts)
                axx.scatter(0.2, yPos, s=50, color=self.colors[orderedLayers[i]])
                axx.set_xlim(0.0, 1)
            axx.axis("off")

        loc = "/home/os18o068/Documents/PHD/Abbildungen/"
        fig.savefig(loc + "hydro_structures.pdf", format="pdf", bbox_inches="tight")
        plt.close()


class PlanetChart:
    def __init__(
        self,
        data,
        radius=1.0,
        label_offset=0.0,
        on_scale=True,
        d1=0.02,
        d2=0.05,
        **kwargs
    ):
        self.radius = radius
        self.ghosts = []
        self.slices = []
        self.data = data
        self.label_offset = label_offset
        self.structure = [int(a) for a in list(self.data.keys())[0]]
        self.nLayers = len(self.structure)
        self.connection_lines = [None for i in range(self.nLayers + 1)]
        self.d1 = d1
        self.d2 = d2 + d1

        # Take upper limit of radius as layer size
        if on_scale:
            self.layerSizes = [a["R_outer"][1] for a in list(self.data.values())[0]]

        # Use uniform layer sizes
        else:
            self.layerSizes = [
                self.radius * (a + 1) / self.nLayers for a in range(self.nLayers)
            ]

        # Take lower limit of radius as gost sizes
        self.ghostSizes = [a["R_outer"][0] for a in list(self.data.values())[0]]
        for i in range(self.nLayers):
            try:
                self.ghostSizes[i] *= radius / self.layerSizes[-1]
            except TypeError:
                print("not gonna happen")

        self.layerSizes = [a / self.layerSizes[-1] * radius for a in self.layerSizes]

        try:
            self.phi1 = kwargs["phi1"]
        except KeyError:
            self.phi1 = np.pi * 1 / 3

        try:
            self.phi2 = kwargs["phi2"]
        except KeyError:
            self.phi2 = np.pi - self.phi1

        try:
            self.colors = kwargs["colors"]
        except KeyError:
            self.colors = [layerColors[layerCodes[int(st)]] for st in self.structure]

    def get_slices(self, res=5):
        self.slices = []
        for i in range(self.nLayers):
            ps = PieSlice(
                r=self.layerSizes[i],
                res=res,
                phi1=self.phi1,
                phi2=self.phi2,
                zorder=self.nLayers - i,
                color=self.colors[i],
                offset=[0.0, -(self.radius)],
            )

            self.slices.append(ps)

    def get_ghosts(self, res=5):
        self.ghosts = []
        for i in range(self.nLayers):
            dcolor = 0.5 * np.min(np.ma.masked_equal(np.array(self.colors[i]), 0))
            print("dcolor =", dcolor)
            if not self.ghostSizes[i] == None:
                gh = PieGhost(
                    r=self.ghostSizes[i],
                    res=res,
                    phi1=self.phi1,
                    phi2=self.phi2,
                    zorder=100,
                    color=[max(c - dcolor, 0.0) for c in self.colors[i]],
                    offset=[0.0, 0],
                )

            else:
                gh = None

            self.ghosts.append(gh)

    def draw(self, res=5, title="Your \ Title", titleFontSize=14, **kwargs):
        try:
            self.ax = kwargs["ax"]
        except KeyError:
            self.fig, self.ax = plt.subplots()

        for ps in self.slices:
            ps.draw(self.ax)

        for gh in self.ghosts:
            try:
                gh.draw(self.ax)
            except AttributeError:
                pass

        """
        for i in range(self.nLayers):
            ps = PieSlice(r = self.layerSizes[i], 
                          res = res, 
                          phi1 = self.phi1,
                          phi2 = self.phi2,
                          zorder = self.nLayers - i,
                          color = self.colors[i],
                          offset = [0., - (self.radius - self.layerSizes[i])])
            
            ps.draw(self.ax)
        """
        self.ax.set_aspect("equal")
        self.ax.axis("off")
        # self.ax.set_title(r'${}$'.format(title), fontsize = titleFontSize)

    def get_label_positions(self, fac=1.1):
        self.all_yPos = []
        self.all_xPos = []
        self.all_yPos_center = [0.0]
        for ls in self.layerSizes:
            self.all_yPos.append(ls * np.sin(self.phi1))
            self.all_xPos.append(ls * np.cos(self.phi1))
            self.all_yPos_center.append(ls)

        self.all_yPos = np.array(self.all_yPos)
        self.all_xPos = np.array(self.all_xPos)
        self.all_raw_xPos = copy.deepcopy(self.all_xPos)
        self.all_raw_yPos = copy.deepcopy(self.all_yPos)
        self.all_yPos -= self.label_offset * self.radius * np.cos(self.phi1)
        self.all_xPos += self.label_offset * self.radius * np.sin(self.phi1)

        # Gather label positions
        for i in range(self.nLayers):

            self.connection_lines[i] = [
                [
                    self.all_raw_xPos[i],
                    self.all_xPos[i] + self.d1 * np.cos(np.pi / 2 - self.phi1),
                ],
                [
                    -1 + self.all_raw_yPos[i],
                    -1 + self.all_yPos[i] - self.d1 * np.sin(np.pi / 2 - self.phi1),
                ],
            ]

        # Add central connection line
        self.connection_lines[-1] = [
            [0, 0 + self.d1 * np.cos(np.pi / 2 - self.phi1)],
            [-1, -1 - self.d1 * np.sin(np.pi / 2 - self.phi1)],
        ]

        adjusted_indices = []
        # Automatically adjust label positions to avoid overlaps
        for i in range(self.nLayers - 1):
            dx_sq = (self.all_xPos[i] - self.all_xPos[i + 1]) ** 2
            dy_sq = (self.all_yPos[i] - self.all_yPos[i + 1]) ** 2
            dpos = np.sqrt(dx_sq + dy_sq)

            drel = dpos / self.radius

            if drel < 0.1:  # or self.all_yPos[i] > self.all_yPos[i + 1]:
                print("adjusting label {}".format(i))
                adjusted_indices.append(i + 1)
                self.all_yPos[i + 1] *= fac
                self.all_xPos[i + 1] *= fac

                self.connection_lines[i + 1] = [
                    [
                        self.all_raw_xPos[i + 1],
                        self.all_xPos[i + 1] + self.d1 * np.cos(np.pi / 2 - self.phi1),
                    ],
                    [
                        -1 + self.all_raw_yPos[i + 1],
                        -1
                        + self.all_yPos[i + 1]
                        - self.d1 * np.sin(np.pi / 2 - self.phi1),
                    ],
                ]

            if i in adjusted_indices:
                print("adjusting label {}".format(i))
                self.all_yPos[i + 1] *= fac
                self.all_xPos[i + 1] *= fac

                self.connection_lines[i + 1] = [
                    [
                        self.all_raw_xPos[i + 1],
                        self.all_xPos[i + 1] + self.d1 * np.cos(np.pi / 2 - self.phi1),
                    ],
                    [
                        -1 + self.all_raw_yPos[i + 1],
                        -1
                        + self.all_yPos[i + 1]
                        - self.d1 * np.sin(np.pi / 2 - self.phi1),
                    ],
                ]

    def add_labels(
        self,
        title="Structure Label",
        labelFontSize=12,
        titleFontSize=16,
        dataMult={
            "P_outer": 1e-9,
            "T_outer": 1e0,
            "R_outer": 1e-3,
            "P_inner": 1e-9,
            "T_inner": 1e0,
            "R_inner": 1e-3,
            "indigenous_mass": 1.0 / m_earth,
            "mass_fraction": 1e2,
            "x_MgO": 1e2,
            "x_FeO": 1e2,
            "x_SiO2": 1e2,
            "Si_core": 1e2,
            "S_core": 1e2,
            "O_core": 1e2,
            "logfO2": 1,
        },
        middleKeys=None,
        outer=True,
    ):

        data_labels = {
            "T_outer": "K",
            "P_outer": "GPa",
            "R_outer": "km",
            "indigenous_mass": r"$\ M_\oplus$",
            "mass_fraction": r"$\ \rm wt\%$",
            "S_core": r"$\rm \ wt\% \ S$",
            "Si_core": r"$\rm \ wt\%  \ Si$",
            "O_core": r"$\rm \ wt\%  \ O$",
            "x_MgO": r"$\rm \% \ MgO$",
            "x_SiO2": r"$\rm \% \ SiO_2$",
            "x_FeO": r"$\rm \% \ FeO$",
            "logfO2": "",
        }

        trafo = transforms.blended_transform_factory(
            self.ax.transData, self.ax.transData
        )

        self.get_label_positions()
        data = list(self.data.values())[0]

        for i in range(self.nLayers):
            left = r"\begin{eqnarray*} "
            right = r"\end{eqnarray*} "
            labelRight = r"&&{a:g}-{b:g} \ \rm K \\ &&{c:g}-{d:g} \ \rm GPa".format(
                a=ftool.fancyround(
                    data[i]["T_outer"][0] * dataMult["T_outer"], digits=3
                ),
                b=ftool.fancyround(
                    data[i]["T_outer"][1] * dataMult["T_outer"], digits=3
                ),
                c=ftool.fancyround(
                    data[i]["P_outer"][0] * dataMult["P_outer"], digits=3
                ),
                d=ftool.fancyround(
                    data[i]["P_outer"][1] * dataMult["P_outer"], digits=3
                ),
            )

            labelLeft = r"${a:g}-{b:g} \ \rm km$".format(
                a=ftool.fancyround(
                    data[i]["R_outer"][0] * dataMult["R_outer"], digits=3
                ),
                b=ftool.fancyround(
                    data[i]["R_outer"][1] * dataMult["R_outer"], digits=3
                ),
            )

            if not middleKeys == None:
                labelMiddle = r""
                for j in range(len(middleKeys[i])):
                    key = middleKeys[i][j]
                    partial_label_middle = r"${a:g}-{b:g}$".format(
                        a=ftool.fancyround(data[i][key][0] * dataMult[key], digits=2),
                        b=ftool.fancyround(data[i][key][1] * dataMult[key], digits=2),
                    )

                    if middleKeys[i][j] == "logfO2_mantle":
                        partial_label_middle = (
                            r"$-\log(f_{O_2})= \ $" + partial_label_middle
                        )

                    labelMiddle = labelMiddle + partial_label_middle
                    labelMiddle = labelMiddle + data_labels[middleKeys[i][j]]
                    if j < len(middleKeys[i]) - 1:
                        labelMiddle = labelMiddle + "\n"

            # label2 =  r"$D_{it} =\begin{cases}\end{cases}$"
            yPos = self.all_yPos[i] - self.radius
            xPos = (self.radius + yPos) * np.cos(self.phi1)
            xPos = self.all_xPos[i]
            yPos = -1 + self.all_yPos[i]

            rotAngle = -(np.pi / 2 - self.phi1) / np.pi * 180

            if not outer and i == self.nLayers - 1:
                pass
            else:
                # Add temperature and pressure ranges at the right
                self.ax.annotate(
                    r"$\}$",
                    fontsize=24,
                    xy=(
                        xPos + self.d1 * np.cos(np.pi / 2 - self.phi1),
                        yPos - self.d1 * np.sin(np.pi / 2 - self.phi1),
                    ),
                    # xycoords='figure fraction',
                    transform=trafo,
                    rotation=rotAngle + 180,
                    va="center",
                    ha="right",
                    rotation_mode="anchor",
                    # bbox=dict(facecolor='none', edgecolor='blue', pad=0.0)
                )

                self.ax.text(
                    xPos + self.d2 * np.cos(np.pi / 2 - self.phi1),
                    yPos - self.d2 * np.sin(np.pi / 2 - self.phi1),
                    left + labelRight + right,
                    rotation=rotAngle,
                    zorder=100,
                    va="center",
                    ha="left",
                    rotation_mode="anchor",
                    # bbox=dict(facecolor='none', edgecolor='blue', pad=0.),
                    transform=trafo,
                    fontsize=labelFontSize,
                )

                # Add radius range at the left
                self.ax.text(
                    -(xPos + self.d1 * np.cos(np.pi / 2 - self.phi1)),
                    yPos - self.d1 * np.sin(np.pi / 2 - self.phi1),
                    labelLeft,
                    # rotation = rotAngle,
                    zorder=100,
                    va="center",
                    ha="right",
                    rotation_mode="anchor",
                    # bbox=dict(facecolor='none', edgecolor='blue', pad=0.),
                    transform=trafo,
                    fontsize=labelFontSize,
                )

            # Add layer mass fraction in the center of each layer
            self.ax.text(
                0,
                -1 + (self.all_yPos_center[i + 1] + self.all_yPos_center[i]) / 2,
                labelMiddle,
                zorder=100,
                va="center",
                ha="center",
                rotation_mode="anchor",
                # bbox=dict(facecolor='none', edgecolor='blue', pad=0.),
                transform=trafo,
                fontsize=labelFontSize,
            )

            # Add central labels
            if i == 0:
                # Add temperature and pressure ranges at the right
                label = r"&&{a:g}-{b:g} \ \rm K \\ &&{c:g}-{d:g} \ \rm GPa".format(
                    a=ftool.fancyround(
                        data[i]["T_inner"][0] * dataMult["T_inner"], digits=3
                    ),
                    b=ftool.fancyround(
                        data[i]["T_inner"][1] * dataMult["T_inner"], digits=3
                    ),
                    c=ftool.fancyround(
                        data[i]["P_inner"][0] * dataMult["P_inner"], digits=3
                    ),
                    d=ftool.fancyround(
                        data[i]["P_inner"][1] * dataMult["P_inner"], digits=3
                    ),
                )
                x_pos_center = self.label_offset * self.radius * np.cos(self.phi1)
                y_pos_center = -1 - self.label_offset * self.radius * np.sin(self.phi1)

                self.ax.annotate(
                    r"$\}$",
                    fontsize=24,
                    xy=(
                        x_pos_center + self.d1 * np.cos(np.pi / 2 - self.phi1),
                        y_pos_center - self.d1 * np.sin(np.pi / 2 - self.phi1),
                    ),
                    # xycoords='figure fraction',
                    transform=trafo,
                    rotation=rotAngle + 180,
                    va="center",
                    ha="right",
                    rotation_mode="anchor",
                    # bbox=dict(facecolor='none', edgecolor='blue', pad=0.0)
                )

                self.ax.text(
                    x_pos_center + self.d2 * np.cos(np.pi / 2 - self.phi1),
                    y_pos_center - self.d2 * np.sin(np.pi / 2 - self.phi1),
                    left + label + right,
                    rotation=rotAngle,
                    zorder=100,
                    va="center",
                    ha="left",
                    rotation_mode="anchor",
                    # bbox=dict(facecolor='none', edgecolor='blue', pad=0.),
                    transform=trafo,
                    fontsize=labelFontSize,
                )

        # Add connection lines for offset labels
        for i in range(self.nLayers + 1):
            print(self.connection_lines)
            if not self.connection_lines[i] == None:
                # Right side
                self.ax.plot(
                    self.connection_lines[i][0], self.connection_lines[i][1], color="k"
                )
                # Left side
                if i < self.nLayers:
                    self.ax.plot(
                        [-cl for cl in self.connection_lines[i][0]],
                        self.connection_lines[i][1],
                        color="k",
                    )

        self.ax.set_ylim(-1.1, 0.1)
        # self.ax.set_title(r'\rm${}$'.format(title), fontsize = titleFontSize)

    def addData(
        self,
        data,
        title="Structure Label",
        labelFontSize=12,
        titleFontSize=16,
        dataMult={
            "layerPressures": 1e-9,
            "layerTemperatures": 1e0,
            "layerSizes": 1e-3,
            "layerMasses": 1e2,
            "MgO_mantle": 1e2,
            "FeO_mantle": 1e2,
            "SiO2_mantle": 1e2,
            "Si_core": 1e2,
            "S_core": 1e2,
            "O_core": 1e2,
            "logfO2_mantle": 1,
        },
        centralKeys=None,
        middleKeys=None,
        central_data=None,
    ):

        data_labels = {
            "layerTemperatures": "K",
            "layerPressures": "GPa",
            "layerSizes": "km",
            "layerMasses": "wt\%",
            "S_core": r"$\rm \ wt\% \ S$",
            "Si_core": r"$\rm \ wt\%  \ Si$",
            "O_core": r"$\rm \ wt\%  \ O$",
            "MgO_mantle": r"$\rm \% \ MgO$",
            "SiO2_mantle": r"$\rm \% \ SiO_2$",
            "FeO_mantle": r"$\rm \% \ FeO$",
            "logfO2_mantle": "",
        }

        trafo = transforms.blended_transform_factory(
            self.ax.transData, self.ax.transData
        )

        all_yPos = []
        all_xPos = []
        all_yPos_center = [0.0]
        for ls in self.layerSizes:
            all_yPos.append(ls * np.sin(self.phi1))
            all_xPos.append(ls * np.cos(self.phi1))
            all_yPos_center.append(ls)

        all_yPos = np.array(all_yPos)
        # self.ax.set_xlim(-1,1)
        # self.ax.set_ylim(-1,1)
        for i in range(self.nLayers):
            left = r"\begin{eqnarray*} "
            right = r"\end{eqnarray*} "
            labelRight = r"&&{a:g}-{b:g} \ \rm K \\ &&{c:g}-{d:g} \ \rm GPa".format(
                a=ftool.fancyround(
                    data["layerTemperatures"][i][0] * dataMult["layerTemperatures"],
                    digits=3,
                ),
                b=ftool.fancyround(
                    data["layerTemperatures"][i][1] * dataMult["layerTemperatures"],
                    digits=3,
                ),
                c=ftool.fancyround(
                    data["layerPressures"][i][0] * dataMult["layerPressures"], digits=3
                ),
                d=ftool.fancyround(
                    data["layerPressures"][i][1] * dataMult["layerPressures"], digits=3
                ),
            )

            labelLeft = r"${a:g}-{b:g} \ \rm km$".format(
                a=ftool.fancyround(
                    data["layerSizes"][i][0] * dataMult["layerSizes"], digits=3
                ),
                b=ftool.fancyround(
                    data["layerSizes"][i][1] * dataMult["layerSizes"], digits=3
                ),
            )
            """
            labelCenter1 = r'${a:g}-{b:g} \ \rm wt\%$'.format(
                a = ftool.fancyround(i+data['layerMasses'][i][0] * \
                                     dataMult['layerMasses'], digits = 2),
                b = ftool.fancyround(data['layerMasses'][i][1] * \
                                     dataMult['layerMasses'], digits = 2))
            
            labelCenter = f'{labelCenter1} \n {labelCenter1}  \n {labelCenter1}'
            """
            if not middleKeys == None:
                labelMiddle = r""

                for j in range(len(middleKeys[i])):
                    dat = data[middleKeys[i][j]][i]

                    partial_label_middle = r"${a:g}-{b:g}$".format(
                        a=ftool.fancyround(
                            dat[0] * dataMult[middleKeys[i][j]], digits=2
                        ),
                        b=ftool.fancyround(
                            dat[1] * dataMult[middleKeys[i][j]], digits=2
                        ),
                    )

                    if middleKeys[i][j] == "logfO2_mantle":
                        partial_label_middle = (
                            r"$-\log(f_{O_2})= \ $" + partial_label_middle
                        )

                    labelMiddle = labelMiddle + partial_label_middle
                    labelMiddle = labelMiddle + data_labels[middleKeys[i][j]]
                    if j < len(middleKeys[i]) - 1:
                        labelMiddle = labelMiddle + "\n"

            # label2 =  r"$D_{it} =\begin{cases}\end{cases}$"
            yPos = all_yPos[i] - self.radius
            xPos = (self.radius + yPos) * np.cos(self.phi1)
            xPos = all_xPos[i]
            yPos = -1 + all_yPos[i]

            d1 = 0.05
            rotAngle = -(np.pi / 2 - self.phi1) / np.pi * 180
            self.ax.annotate(
                r"$\}$",
                fontsize=24,
                xy=(xPos, yPos),
                # xycoords='figure fraction',
                transform=trafo,
                rotation=rotAngle + 180,
                va="center",
                ha="right",
                rotation_mode="anchor",
                # bbox=dict(facecolor='none', edgecolor='blue', pad=0.0)
            )

            # Add temperature and pressure ranges at the right
            self.ax.text(
                xPos + d1 * np.cos(np.pi / 2 - self.phi1),
                yPos - d1 * np.sin(np.pi / 2 - self.phi1),
                left + labelRight + right,
                rotation=rotAngle,
                zorder=100,
                va="center",
                ha="left",
                rotation_mode="anchor",
                # bbox=dict(facecolor='none', edgecolor='blue', pad=0.),
                transform=trafo,
                fontsize=labelFontSize,
            )

            # Add radius range at the left
            self.ax.text(
                -(xPos + 0.01),
                yPos,
                labelLeft,
                # rotation = rotAngle,
                zorder=100,
                va="center",
                ha="right",
                rotation_mode="anchor",
                # bbox=dict(facecolor='none', edgecolor='blue', pad=0.),
                transform=trafo,
                fontsize=labelFontSize,
            )

            if i >= 0 and not central_data == None:
                # Add layer mass fraction in the center of each layer
                self.ax.text(
                    0,
                    -1 + (all_yPos_center[i + 1] + all_yPos_center[i]) / 2,
                    labelMiddle,
                    zorder=100,
                    va="center",
                    ha="center",
                    rotation_mode="anchor",
                    # bbox=dict(facecolor='none', edgecolor='blue', pad=0.),
                    transform=trafo,
                    fontsize=labelFontSize,
                )
            # Add central labels
            if i == 0:
                # Add temperature and pressure ranges at the right
                label = r"&&{a:g}-{b:g} \ \rm K \\ &&{c:g}-{d:g} \ \rm GPa".format(
                    a=ftool.fancyround(
                        central_data["Temperature"][0] * dataMult["layerTemperatures"],
                        digits=3,
                    ),
                    b=ftool.fancyround(
                        central_data["Temperature"][1] * dataMult["layerTemperatures"],
                        digits=3,
                    ),
                    c=ftool.fancyround(
                        central_data["Pressure"][0] * dataMult["layerPressures"],
                        digits=3,
                    ),
                    d=ftool.fancyround(
                        central_data["Pressure"][1] * dataMult["layerPressures"],
                        digits=3,
                    ),
                )

                self.ax.annotate(
                    r"$\}$",
                    fontsize=24,
                    xy=(0, -1),
                    # xycoords='figure fraction',
                    transform=trafo,
                    rotation=rotAngle + 180,
                    va="center",
                    ha="right",
                    rotation_mode="anchor",
                    # bbox=dict(facecolor='none', edgecolor='blue', pad=0.0)
                )
                self.ax.text(
                    d1 * np.cos(np.pi / 2 - self.phi1),
                    -1 - d1 * np.sin(np.pi / 2 - self.phi1),
                    left + label + right,
                    rotation=rotAngle,
                    zorder=100,
                    va="center",
                    ha="left",
                    rotation_mode="anchor",
                    # bbox=dict(facecolor='none', edgecolor='blue', pad=0.),
                    transform=trafo,
                    fontsize=labelFontSize,
                )

        self.ax.set_ylim(-1.1, 0.1)
        # self.ax.set_title(r'\rm${}$'.format(title), fontsize = titleFontSize)


class PieGhost:
    def __init__(self, r=1.0, res=5, zorder=0, color="b", offset=[0, 0], **kwargs):
        try:
            self.phi1 = kwargs["phi1"]
        except KeyError:
            self.phi1 = np.pi * 1 / 3

        try:
            self.phi2 = kwargs["phi2"]
        except KeyError:
            self.phi2 = np.pi * 2 / 3

        self.r = r  # Radius of circle
        self.angles = np.linspace(self.phi1, self.phi2, 2 ** res)
        self.offset = offset

        # Create circle segment
        self.x = self.r * np.cos(self.angles)
        self.y = self.r * np.sin(self.angles)

        self.color = color
        self.zorder = zorder

    def draw(self, axis, loc=[0, 0], orientation="normal", zorder=0):

        # Draw cirle segment
        axis.plot(
            self.x, self.y - 1.0, color=self.color, linestyle="--", zorder=self.zorder
        )


class PieSlice:
    def __init__(self, r=1.0, res=5, zorder=0, color="b", offset=[0, 0], **kwargs):
        try:
            self.phi1 = kwargs["phi1"]
        except KeyError:
            self.phi1 = np.pi * 1 / 3

        try:
            self.phi2 = kwargs["phi2"]
        except KeyError:
            self.phi2 = np.pi * 2 / 3

        self.r = r  # Radius of circle
        self.angles = np.linspace(self.phi1, self.phi2, 2 ** res)
        self.offset = offset

        # Create circle segment
        self.x = self.r * np.cos(self.angles)
        self.y = self.r * np.sin(self.angles) - self.r * np.cos(np.pi / 2 - self.phi1)

        self.color = color
        self.zorder = zorder
        self.dy = self.r * np.cos(np.pi / 2 - self.phi1)
        # Create triangle
        self.polygon = Polygon(
            [
                (self.x[-1] + self.offset[0], -1 + self.dy),
                (0 + self.offset[0], -1),
                (self.x[0] + self.offset[0], -1 + self.dy),
            ],
            color=self.color,
            zorder=self.zorder,
        )

    def draw(self, axis, loc=[0, 0], orientation="normal", zorder=0):
        # Draw triangle
        axis.add_patch(self.polygon)

        # Draw cirle segment
        axis.fill_between(
            self.x,
            self.offset[1] + self.dy,
            self.y + self.offset[1] + self.dy,
            color=self.color,
            zorder=self.zorder,
        )
