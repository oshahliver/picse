#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 10:04:59 2019

@author: oshah
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt
from PIMPphysicalparams import molar_mass_list, rho0_coeffs_Mg2SiO4,\
                                K0_coeffs_Mg2SiO4, saturation_coeffs_Mg2SiO4,\
                                phase_transition_coeffs_Mg2SiO4, \
                                delta_rho_MgO_coeffs, mFe, mOl, mMg,\
                                G0_coeffs_Mg2SiO4, Rgas

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap     
from PIMPrunparams import color_list, grid_color, plot_params, background_color
import functionTools as ftool
import plotTools


#Mantle geotherm from Stixrude 2009
#[P (GPa), T (K)]

P_mantle_geotherm = [
0.28282828282828376,
2.326797385620921,
4.727272727272727,
7.484254307783729,
11.165775401069524,
14.279263220439695,
14.913844325609041,
16.56803327391563,
18.935234699940594,
20.461081402257875,
22.616755793226396,
24.432560903149138,
24.850861556743922,
26.894830659536545,
29.830065359477118,
33.65656565656566,
38.55258467023174,
43.09209744503862
]

T_mantle_geotherm = [
426.6666666666679,
560,
693.3333333333358,
826.6666666666679,
1013.3333333333321,
1146.6666666666679,
1306.6666666666679,
1386.6666666666679,
1466.6666666666679,
1626.6666666666679,
1653.3333333333321,
1706.6666666666679,
1520.0000000000036,
1653.3333333333321,
1786.6666666666679,
1920,
2053.3333333333358,
2186.666666666668]

#effect of Mg# number on water saturation content in Mg2SiO4
#conventions are: temp (°C), pres (GPa), Mg# (mol %), X_H2O (wt %), phase
#source
data_unsaturated_Mg2SiO4 = np.array([
                                    [[1400, 14.1, 100, 0., 0],
                                     [1400, 14.8, 100, 0., 1],
                                     [1400, 14.8, 100, 0., 1],
                                     [1400, 15.6, 100, 0., 1],
                                     [1400, 19.5, 100, 0., 2],
                                     [1400, 20.0, 100, 0., 2], 'Wu Y 2012'],
                                     [[1600, 20.8, 100, 0., 2],
                                      [1600, 19.2, 100, 0., 1],
                                      [1600, 19.6, 100, 0., 1],
                                      [1600, 12.3, 90, 0., 0],
                                      [1600, 16., 90, 0., 1],
                                      [1600, 17.8, 100, 1],
                                      [1600, 12.6, 100, 0],
                                      [1600, 14.1, 80, 0., 1],
                                      [1600, 19.9, 100, 0., 2],
                                      [1600, 12.2, 45, 0., 2],
                                      [1600, 12., 45, 0., 2],
                                      [1600, 12.9, 55, 0., 2],
                                      [1600, 11.1, 45, 0., 2],
                                      [1600, 14.6, 100, 0., 1],
                                      [1600, 15.0, 96, 0., 1],
                                      [1600, 14.1, 85, 0., 1],
                                      [1600, 13.7, 96, 0., 1],
                                      [1600, 14.6, 96, 0., 1],
                                      [1600, 15.7, 100, 0., 1],
                                      [1600, 15., 100, 0., 1],
                                      [1600, 14.1, 100, 0., 0],
                                      [1600, 14.6, 100, 0., 0],
                                      [1600, 15., 100, 0., 1],
                                      [1600, 13.7, 91, 0., 0],
                                      [1600, 14.1, 91, 0., 0],
                                      [1600, 15.3, 100, 0., 1],
                                      [1600, 15., 100, 0., 0],
                                      [1600, 13.2, 80, 0., 0],
                                      [1600, 13.2, 86, 0., 0],
                                      [1600, 13.7, 86, 0., 0],
                                      [1600, 14.6, 91, 0., 1],
                                      [1200, 14.4, 100, 0., 1],
                                      [1200, 14.4, 80, 0., 2],
                                      [1200, 13.5, 100, 0., 0],
                                      [1200, 13.5, 91, 0., 1],
                                      [1200, 13.9, 100, 0., 0],
                                      [1200, 12.2, 85, 0., 0],
                                      [1200, 17.3, 100, 0., 1],
                                      [1200, 18.6, 100, 0., 1],
                                      [1200, 19.2, 100, 0., 2],
                                      [1200, 19.2, 96, 0., 2],
                                      [1600, 14.4, 90, 0., 3],
                                      [1600, 14.5, 70, 0., 4],
                                      [1600, 19.2, 90, 0., 4], 'Katsura 1989'],
                                     [[1100, 14.87, 100, 0., 3],
                                      [900, 14.88, 100, 0., 3],
                                      [650, 14.83, 100, 0., 5],
                                      [875, 14.92, 100, 0., 3],
                                      [650, 15.16, 100, 0., 5],
                                      [750, 15.2, 100, 0., 4],
                                      [1100, 15.82, 100, 0., 4],
                                      [985, 15.83, 100, 0., 3],
                                      [800, 16., 18, 100, 0., 5],
                                      [880, 16.2, 100, 0., 4],
                                      [700, 16.2, 100, 0., 5],
                                      [1200, 16.62, 100, 0., 4],
                                      [1050, 17.07, 100, 0., 4],
                                      [900, 17.18, 100, 0., 5],
                                      [1000, 17.32, 100, 0., 5],
                                      [1200, 17.91, 100, 0., 4],
                                      [1500, 18.2, 100, 0., 4],
                                      [800, 18.13, 100, 0., 4],
                                      [930, 18.85, 100, 0., 5],
                                      [700, 19.34, 100, 0., 5], 'Suzuki 2000'],
                                      [[1200, 12.62, 81.1, 0., 1],
                                       [1200, 12.62, 68.4, 0., 2],
                                       [1200, 14.28, 87.2, 0., 1],
                                       [1200, 14.28, 77.8, 0., 2],
                                       [1200, 15.68, 93.2, 0., 1],
                                       [1200, 15.68, 88.3, 0., 2],
                                       [1400, 14.27, 83.2, 0., 1],
                                       [1400, 14.27, 75.5, 0., 2],
                                       [1400, 14.69, 85.1, 0., 1],
                                       [1400, 14.69, 78.1, 0., 2],
                                       [1400, 14.98, 85.3, 0., 1],
                                       [1400, 14.98, 77.2, 0., 2],
                                       [1400, 15.13, 85.5, 0., 1],
                                       [1400, 15.3, 77.6, 0., 2],
                                       [1400, 16.46, 92.9, 0., 1],
                                       [1400, 16.46, 88.4, 0., 2],
                                       [1600, 13.52, 75.1, 0., 1],
                                       [1600, 13.52, 65.5, 0., 2],
                                       [1600, 14.62, 80.6, 0., 1],
                                       [1600, 14.62, 71.3, 0., 2],
                                       [1600, 14.98, 81.2, 0., 1],
                                       [1600, 14.98, 73.0, 0., 2],
                                       [1600, 15.76, 84.5, 0., 1],
                                       [1600, 15.76, 77.3, 0., 2],
                                       [1600, 17.98, 93.3, 0., 1],
                                       [1600, 17.98, 89.6, 0., 2], 'Tsujino 2019'],
                                       [[1400, 14.6, 100-17.3, 0., 2],
                                        [1400, 15.6, 100-16.7, 0., 2],
                                        [1400, 15.7, 100-16.1, 0., 2],
                                        [1400, 15.9, 100-15.3, 0.2], 'Inoue 2010']
                                    
                                    ])


data_saturation_Mg2SiO4 = np.array([[#alpha phase
                        [[1100, 2.5, 100, 0.0135],
                         [1100, 5., 100, 0.0496],
                         [1100, 6.5, 100, 0.026],
                         [1100, 8., 100, 0.0867],
                         [1000, 8., 100, 0.086],
                         [1100, 8., 100, 0.098],
                         [1100, 8., 100, 0.049],
                         [1100, 9., 100, 0.0984],
                         [1100, 10., 100, 0.107],
                         [1100, 12., 100, 0.151],
                         [1100, 13., 100, 0.109], 'Kohlstedt 1996'],
                        [[1100, 13.7, 100, 0.76],
                         [1100, 12.95, 100-3, 0.29],
                         [1100, 13., 100-4.1, .37],
                         [1100, 13.05, 100-3.7, .71],
                         [1100, 13.1, 100-3.2, .64],
                         [1100, 13.1, 100-2.9, .55], 'Chen 2002'],
                        [[1100, 6., 100, 0.0548],
                         [1100, 9., 100, 0.1223], 
                         [1100, 2.5, 100, 0.0051],
                         [1250, 2.5, 100, 0.0051],
                         [1325, 2.5, 100, 0.0111],
                         [1400, 2.5, 100, 0.0173], 
                         [1100, 6., 100, 0.1189],
                         [1175, 6., 100, 0.078],
                         [1250, 6., 100, 0.0607],
                         [1325, 6., 100, 0.0551],
                         [1400, 6., 100, 0.0489], 
                         [1175, 9., 100, 0.1567],
                         [1250, 9., 100, 0.1984],
                         [1325, 9., 100, 0.1154],
                         [1400, 9., 100, 0.0386], 'Bali 2008'],
                        [[1100, 13., 100, 0.2], 'Smyth 2003'],
                        [[1200, 12., 100, 0.15], 'Chen 1998'],
                        [[1400, 5, 100, 0.0018],
                         [1400, 13.5, 93.8, 0.168],
                         [1175, 2.5, 92.6, 0.0126],
                         [1250, 2.5, 92.6, 0.0054],
                         [1325, 2.5, 93.2, 0.0032],
                         [1175, 2.5, 91.7, 0.0235],
                         [1250, 2.5, 91.9, 0.0184],
                         [1325, 2.5, 88.8, 0.0065],
                         [1400, 2.5, 93.9, 0.00015], 
                         [1175, 5.0, 92.9, 0.0263],
                         [1250, 5.0, 92.2, 0.0286],
                         [1250, 5.0, 92, 0.0552],
                         [1325, 5.0, 91, 0.0158],
                         [1400, 5.0, 94.7, 0.0129], 
                         [1175, 7.5, 92.7, 0.2072],
                         [1250, 7.5, 92.3, 0.1262],
                         [1325, 7.5, 93.4, 0.0964],
                         [1400, 7.5, 92.7, 0.083],
                         [1175, 7.5, 94, 0.0904],
                         [1250, 7.5, 92, 0.0876],
                         [1325, 7.5, 90.3, 0.0593],
                         [1400, 7.5, 87.2, 0.0284], 
                         [1175, 9., 90.6, 0.142],
                         [1175, 9., 91.8, 0.4690],
                         [1250, 9., 91.7, 0.2805],
                         [1400, 9., 93.2, 0.0046], 'Ferot 2012'],
                        [[1000, 1.5, 100, 6.4e-4],
                         [1000., .2, 100, 2.8e-4],
                         [1060., .2, 100, 5.7e-4],
                         [1060., .2, 100, 5.5e-4],
                         [1100., .2, 100, 2.3e-4],
                         [1110., .2, 100, 5.7e-4],
                         [900., .2, 100, 1.3e-4],
                         [950., .2, 100, 4.2e-4],
                         [900., .2, 100, 3.8e-4], 'Demouchy 2003'],
                        [[1600, 6.3, 100, 0.0915],
                         [1600, 6.3, 100, 0.095],
                         [1400, 6.3, 100, 0.1],
                         [1200, 6.3, 100, 0.104], 'Sokol 2010'],
           
                        [[1000, 0.3, 100.0, 0.000717],
                         [1000, 0.3, 100-8.5, 0.00224],
                         [1000, 0.3, 100-12.0, 0.00323],
                         [1000, 0.3, 100-14.9, 0.00381],
                         [1000, 0.3, 100-15.3, 0.00386],
                         [1000, 0.3, 100-16.8, 0.00439],
                         [1000, 0.3, 100-16.9, 0.00417],
                         
                         [1050, 0.3, 100-0.0, 0.000941],
                         [1050, 0.3, 100-8.5, 0.00256],
                         [1050, 0.3, 100-12.0, 0.00435],
                         [1050, 0.3, 100-14.9, 0.00502],
                         [1050, 0.3, 100-15.3, 0.00538],
                         [1050, 0.3, 100-16.8, 0.00574],
                         [1050, 0.3, 100-16.9, 0.00583],
                         
                         [1100, 0.3, 100-0.0, 0.000986],
                         [1100, 0.3, 100-8.5, 0.00318],
                         [1100, 0.3, 100-12.0, 0.00453],
                         [1100, 0.3, 100-14.9, 0.00511],
                         [1100, 0.3, 100-15.3, 0.00547],
                         [1100, 0.3, 100-16.8, 0.0065],
                         
                         [1150, 0.3, 100-0.0, 0.00126],
                         [1150, 0.3, 100-8.5, 0.00368],
                         [1150, 0.3, 100-12.0, 0.0048],
                         [1150, 0.3, 100-14.9, 0.00673],
                         [1150, 0.3, 100-15.3, 0.00664],
                         [1150, 0.3, 100-16.8, 0.00803],
                         [1150, 0.3, 100-16.9, 0.00771],
                         
                         [1200, 0.3, 100-0.0, 0.0013],
                         [1200, 0.3, 100-8.5, 0.00439],
                         [1200, 0.3, 100-12.0, 0.00592],
                         [1200, 0.3, 100-14.9, 0.0069],
                         [1200, 0.3, 100-15.3, 0.00731],
                         [1200, 0.3, 100-16.8, 0.00883],
                         [1200, 0.3, 100-16.9, 0.00915],
                         
                         [1250, 0.3, 100-0.0, 0.00157],
                         [1250, 0.3, 100-8.5, 0.00448],
                         [1250, 0.3, 100-12.0, 0.00682],
                         [1250, 0.3, 100-14.9, 0.00771],
                         [1250, 0.3, 100-15.3, 0.00771],
                         [1250, 0.3, 100-16.8, 0.00933],
                         [1250, 0.3, 100-16.9, 0.0096],
                         
                         [1300, 0.3, 100-0.0, 0.00184],
                         [1300, 0.3, 100-8.5, 0.00551],
                         [1300, 0.3, 100-12.0, 0.00731],
                         [1300, 0.3, 100-15.3, 0.00892],
                         [1300, 0.3, 100-16.8, 0.0103],
                         [1300, 0.3, 100-16.9, 0.0105], 'Zhao 2004']
                        ],
                        #beta phase
                       [[[1100, 12.5, 88.2, 2.03], 
                        [1100, 13.5, 91.2, 2.12],
                        [1200, 12.5, 88.2, 1.65],
                        [1200, 13.5, 90., 1.98],
                        [1400, 14.5, 86.7, 1.60],
                        [1400, 15.2, 94., 1.37],
                        [1400, 17.0, 94.1, 0.63],
                        [1600, 14.5, 96.9, 0.33],
                        [1600, 15.2, 96.0, 0.31],
                        [1800, 20.0, 93.0, 0.15],
                        [1900, 20.0, 94.5, 0.12],
                        [2000, 16.5, 94.9, 0.07],'Litasov 2008'],
                        [[1200, 17.9, 96.0, 1.8],
                        [1200, 16.8, 88.0, 1.9],
                        [1200, 16.2, 86.0, 1.2],
                        [1200, 17.3, 89.0, 1.7],
                        [1200, 16.8, 89.0, 1.0],
                        [1200, 19.0, 96.0, 1.3],
                        [1200, 15.2, 78.0, 1.4],
                        [1200, 17.9, 89.0, 1.6], 'Mrosko 2015'],
                        [[1400, 17.0, 100-14.5, 3.72],
                         [1400, 16.5, 100-13.9, 2.28],
                         [1400, 16.5, 100-13.9, 2.24],
                         [1400, 16.5, 100-14.0, 1.88],
                         [1400, 16.5, 100-14.0, 1.79], 'Inoue 2010'],
                        [[900, 15., 100, 2.23],
                         [1000, 15., 100, 2.13],
                         [1100, 15., 100, 2.41],
                         [1100, 15., 100, 2.42],
                         [1200, 15., 100, 2.24],
                         [1300, 15., 100, 1.66],
                         [1400, 15., 100, 0.93], 
                         [1200, 14., 100, 2.4],
                         [1200, 14., 100, 2.68], 
                         [1200, 16, 100, 2.6],
                         [1200, 17, 100, 2.43],
                         [1200, 18, 100, 1.24], 'Demouchy 2005'],
                        [[1500, 15.5, 100, 1.53],
                         [1360., 15.5, 100, 3.13],
                         [1450, 15.5, 100, 1.95],
                         [1300, 15.5, 100, 1.83], 
                         [1400, 13.5, 100, 2.89],
                         [1400, 16.5, 100, 2.61],
                         [1400, 15.5, 100, 3.13], 'Kawamoto 1996'],
                        [[1400, 13.5, 100, 0.79], 'Ferot 2012'],
                        [[1300., 15., 100, 1.5], 'Kohn 2002'],
                        [[1100, 15., 100, 2.32],
                         [1100, 15., 100, 2.41], 
                         [1100, 14., 100, 2.13], 'Kohlstedt 1996'],
                        [[1100, 14.2, 100, 3.8],
                         [1100, 14.7, 100, 3.5], 
                         [1100, 13., 100, 1.9],
                         [1100, 13.05, 100, 3.4],
                         [1100, 13.1, 100, 3.4],
                         [1100, 13.1, 100, 2.3],
                         [1100, 13.25, 100, 3.],
                         [1100, 13.6, 100, 2.6],
                         [1100, 13.7, 100, 2.8],
                         [1100, 13.9, 100, 2.2],
                         [1100, 14.2, 100, 2.6], 'Chen 2002'],
                        [[1150., 14.25, 100, 2.9],
                         [1150., 14.25, 100, 3.2], 'Griffin 2014'],
                        [[1400., 15., 100, 0.9],
                         [1400, 16., 100, 0.84], 'Mao 2008'],
                        [[2100, 18., 100, 0.005],
                         [1400, 17., 100, 0.015],
                         [1200, 16., 100, 0.32],
                         [1300, 16., 100, 0.6],
                         [1400, 16., 100, 0.96],
                         [1400, 16., 100, 1.06], 'Jacobsen 2005'],
                        [[1200, 13.3, 100, .8],
                         [1150, 13.5, 100, 1.6], 'Deon 2010'],
                        [[1200, 13.5, 100, 3.8], 'Chen 1998'],
                        [[1300, 15., 100, 2.4], 'Inoue 2004'],
                        [[1300, 15., 100, 0.2212], 'Bolfan 2000'],
                        [[1300, 14., 100, 1.66], 'Holl 2008'],
                        [[1300, 15.5, 100, 2.5], 'Yusa 1997'],
                        
                        ],
                         
                         #gamma phase
                       [[[1400, 17.0, 100-22.4, 1.69],
                         [1400, 16.5, 100-21.3, 1.11],
                         [1400, 16.5, 100-21.3, 1.25],
                         [1400, 16.5, 100-20.8, 1.0],
                         [1400, 16.5, 100-20.8, 1.1],
                         [1600, 23.0, 100-15.9, 0.76],
                         [1600, 32.1, 100-13.5, 0.71],
                         [1600, 23.2, 100-8.7, 0.63], 'Inoue 2010'],
                        [[1300, 19., 100, 2.2], 'Inoue 1998'],
                        [[1100., 19., 100, 2.7], 'Kohlstedt 1995'],
                        [[1100., 19.5, 100, 2.62], 'Kohlstedt 1996'],
                        [[1200., 19., 100, 0.7817],
                         [1300., 19., 100, 0.75], 'Bolfan 2000'],
                        [[1250., 20., 100, 1.77],
                         [1250., 20., 100, 1.65],
                         [1250., 20., 100, 2.5],
                         [1400., 20., 100, 1.13],
                         [1400., 20., 100, 0.85], 
                         [1400, 18., 100, 0.94],
                         [1400, 18., 100, 1.11], 'Thomas 2015'],
                        [[1400., 20., 100, 0.7892],
                         [1400., 20., 100, 0.9263], 
                         [1500, 21., 100, 0.1988],
                         [1500, 22., 100, 0.7358],
                         [1500, 21.5, 100, 0.7354],
                         [1400, 19., 100, 0.8557],
                         [1400, 18., 100, 1.0661], 'Smyth 2003'],
                        [[1300., 20., 100, 2.2], 
                         [1300, 20., 100, 2.0], 'Kudoh 2000'],
                        [[1300, 19., 100, 2.8], 'Yusa 2000'],
                        [[1300, 19., 100, 2.2], 'Chen 1998']
                       ]])



#Define plot symbols for each data set
plot_symbols = {'Chen 1990':['s','col'], 
                'Yusa 2000':['o','col'], 
                'Kudoh':['v','col'], 
                'Smyth 2003':['<','col'],
                'Thomas 2015':['>','col'], 
                'Kohlstedt 1995':['p','col'],
                'Bolfan 2000':['x','col'], 
                'Inoue 1998':['P','col'], 
                'Inoue 2010':['*','col'], 
                'Yusa 1997':['h','col'],
                'Holl 2008':['+','col'], 
                'Inoue 2004':['D','col'], 
                'Chen 1998':['s',''],
                'Deon 2010':['o',''],
                'Jacobsen 2005':['v',''],
                'Mao 2008':['<',''],
                'Griffin 2014':['>',''],
                'Chen 2002':['p',''],
                'Kohlstedt 1996':['P',''],
                'Kohn 2002':['*',''],
                'Ferot 2012':['h',''],
                'Kawamoto 1996':['D',''],
                'Demouchy 2005':['^','col'],
                'Mrosko 2015':['^',''],
                'Litasov 2008':['d','col'],
                'Sokol 2010':['d',''],
                'Bali 2008':['X','col']}

#-----
#phase transitions
#------
#temperature in °C and pressure in GPa
#transition from olivine to wadsleyite
#data from Katsura 1989
P_alpha_beta_Katsura1989 = np.array([15.0, 
                                     14.0,
                                     11.58])

T_alpha_beta_Katsura1989 = np.array([1600., 
                                     1200.,
                                     15.])

#data from Suzuki et al. 2000
P_alpha_beta_Suzuki2000 = np.array([14.88,
                                    14.87,
                                    14.95,
                                    14.88,
                                    15.01,
                                    14.92,
                                    15.88,
                                    15.83])

T_alpha_beta_Suzuki2000 = np.array([1100,
                                    1100,
                                    900,
                                    900,
                                    875,
                                    875,
                                    985,
                                    985])



P_beta_gamma_Suzuki2000 = np.array([15.32,
                                    15.2,
                                    15.83,
                                    15.82,
                                    16.3,
                                    16.2,
                                    16.59,
                                    16.62,
                                    ])

T_beta_gamma_Suzuki2000 = np.array([750.,
                                    750,
                                    1100,
                                    1100,
                                    800,
                                    800,
                                    1200,
                                    1200])



P_alpha_gamma_Suzuki2000 = np.array([14.98,
                                     14.83,
                                     15.32,
                                     15.16,
                                     16.3,
                                     16.18,
                                     16.37,
                                     16.2,
                                     17.26,
                                     17.18,
                                     17.37,
                                     17.32,
                                     18.91,
                                     18.85,
                                     19.5,
                                     19.34])

T_alpha_gamma_Suzuki2000 = np.array([650,
                                     650,
                                     650,
                                     650,
                                     800,
                                     800,
                                     700,
                                     700,
                                     900,
                                     900,
                                     1000,
                                     1000,
                                     930,
                                     930,
                                     700,
                                     700])



#data from Katsura 2004
P_alpha_beta_Katsura2004 = np.array([14.2, 15.4])
T_alpha_beta_Katsura2004 = np.array([1600., 1900.])


#data from Liu 2010
P_alpha_beta_Liu2010 = np.array([11.58])
T_alpha_beta_Liu2010 = np.array([0.])

#data from Wu Yao 2011
P_alpha_beta_WuYao2011 = np.array([14.8])
T_alpha_beta_WuYao2011 = np.array([1400.])

#transition from wadsleyite to ringwoodite
#data from Tsujino 2019
P_beta_gamma_Tsujino2019  = np.array([])*1.0
T_beta_gamma_Tsujino2019  = np.array([])

#data from Wu Yao 2011
P_beta_gamma_WuYao2011 = np.array([19.5])
T_beta_gamma_WuYao2011 = np.array([1400.])

#data from Katsura 1989
P_beta_gamma_Katsura1989 = np.array([14.5, 
                                     19.2, 
                                     16.5, 
                                     13.2, 
                                     13.9, 
                                     12.6, 
                                     17.3, 
                                     18.6, 
                                     21.0, 
                                     19.0, 
                                     13.5])

T_beta_gamma_Katsura1989 = np.array([1600., 
                                     1600., 
                                     1600., 
                                     1600., 
                                     1200., 
                                     1200., 
                                     1200.,
                                     1200.,
                                     1600., 
                                     1200., 
                                     1600.])

P_beta_gamma_Katsura1989 = np.array([21.0, 19.0])
T_beta_gamma_Katsura1989 = np.array([1600., 1200.])

#data from Inoue 2006
P_beta_gamma_inoue2006 = np.array([18.9, 19.4, 18.5, 18.8])
T_beta_gamma_inoue2006 = np.array([1400., 1500., 1300., 1450.])

#gather all transition data points in one array
P_alpha_beta_data = np.concatenate((P_alpha_beta_Katsura1989, 
                               P_alpha_beta_Katsura2004, 
                               P_alpha_beta_WuYao2011,
                               P_alpha_beta_Liu2010))


T_alpha_beta_data = np.concatenate((T_alpha_beta_Katsura1989, 
                               T_alpha_beta_Katsura2004, 
                               T_alpha_beta_WuYao2011,
                               T_alpha_beta_Liu2010))


P_beta_gamma_data = np.concatenate((P_beta_gamma_Katsura1989,
                               P_beta_gamma_Tsujino2019, 
                               P_beta_gamma_WuYao2011,
                               P_beta_gamma_inoue2006,
                               P_beta_gamma_Suzuki2000))

T_beta_gamma_data = np.concatenate((T_beta_gamma_Katsura1989,
                               T_beta_gamma_Tsujino2019, 
                               T_beta_gamma_WuYao2011,
                               T_beta_gamma_inoue2006,
                               T_beta_gamma_Suzuki2000))


'''
#perform linear fit to transition data in °C and Pa
coeffs_alpha_beta = np.polyfit(T_alpha_beta_data, P_alpha_beta_data, 1)
fit_alpha_beta = np.poly1d(coeffs_alpha_beta)

coeffs_beta_gamma = np.polyfit(T_beta_gamma_data, P_beta_gamma_data, 1)
fit_beta_gamma = np.poly1d(coeffs_beta_gamma)
'''

def P_beta_gamma_fit(T):
    """linear fit to beta-gamma phase transition data from various authors.
    Returns transition pressure in GPa as function of temperature in °C
    """
    c1, c2 = phase_transition_coeffs_Mg2SiO4[1]
    
    return (c1 + c2*T)

def P_alpha_beta_fit(T):
    """linear fit to alpha-beta phase transition data from various authors.
    Returns transition pressure in GPa as function of temperature in °C
    """
    c1, c2 = phase_transition_coeffs_Mg2SiO4[0]
    
    return (c1 + c2*T)

def T_beta_gamma_fit(P):
    """input pressure in GPa, output temperature in °C
    """
    c1, c2 = phase_transition_coeffs_Mg2SiO4[1]
    
    return max((P -  c1)/c2, 0.)

def T_alpha_beta_fit(P):
    """input pressure in GPa, output temperature in °C
    """
    c1, c2 = phase_transition_coeffs_Mg2SiO4[0]
    
    return max((P -  c1)/c2, 0.)


def phaseCheck(T=None, P=None, Fe_number=0.0):
    """Takes pressure in Pa and temperature in K as input and returns phase of 
    (Mg,Fe)2SiO4 according to the linear fits to experimental data. Note, that for
    (Mg,Fe)2SiO4 there are regions of divariant and univariant coexistence as 
    described e.g. in  Katsura 1989. This regions are neglected here.
    The effect of Fe on the phase transitions is estimated from Inoue 2010 for
    relevant iron contents of Fe/Mg < 0.2. The effect of water on the phase 
    diagram is neglegted for the relevant range of water contents.
    """

    #print ('intercept temperature =', round(T_intercept,1), '°C')
                    
    #compute transition pressures at given temperature in Pa
    #Note that the transition fits take T in °C
    trans1 = P_alpha_beta_fit(T-273.15)*1.0e9
    trans2 = P_beta_gamma_fit(T-273.15)*1.0e9
    print ('1,2 =', trans1, trans2)
    #check in which regime given pressure is
    
    #alpha phase
    if P < trans1:
        return 0
    
    #beta phase
    elif P >= trans1:
        if trans1 < trans2:
            if P <= trans2:
                return 1
            
            else:
                return 2
            
        elif trans1 >= trans2:
            return 2


def P_phase_range(T=None, ph=None):
    
    #alpha phase
    if ph == 0:
        P_min = 0.
        P_max = P_alpha_beta_fit(T)
        
    #beta phase
    elif ph == 1:
        P_min = P_alpha_beta_fit(T)
        P_max = P_beta_gamma_fit(T)
        
    #gamma phase
    elif ph == 2:
        P_min = P_beta_gamma_fit(T)
        P_max = 23.
        
    return np.array([P_min, P_max])


def T_phase_range(P=None, ph=None):
    """input pressure in GPa, output temperature in °C
    """
    
    #compute intercept of alpha-beta and beta-gamma line
    a1 = phase_transition_coeffs_Mg2SiO4[0][0]
    a2 = phase_transition_coeffs_Mg2SiO4[1][0]
    b1 = phase_transition_coeffs_Mg2SiO4[0][1]
    b2 = phase_transition_coeffs_Mg2SiO4[1][1]
    
    T_intercept = (a1-a2)/(b2-b1)
                    
    #alpha phase
    if ph == 0:
        T_min = T_alpha_beta_fit(P)
        T_max = 2000.
        
    #beta phase
    elif ph == 1:
        T_min = max(T_beta_gamma_fit(P), T_intercept)
        T_max = max(T_alpha_beta_fit(P), T_intercept)
        
    #gamma phase
    elif ph == 2:
        T_min = 800.
        T_max = T_beta_gamma_fit(P)
        
    return np.array([T_min, T_max])

def rho0_forsterit(C_H2O):
    """Compute ambient density of forsterit as function of water wt% H2O
    according to Jacobsen et al. 2008
    """
    return 3225. *(1.-0.014*C_H2O/0.0089)


def fit_FeX():
    pres = np.array([.0, 18.4, 18.4, 29. , 29. , 41. , 41. , 64. , 64. , 82. , 82. ])
    X = np.array([0., 0.63, 0.59, 0.98, 0.94, 0.97, 0.95, 1.05, 1.06, 1.04, 1.06])
    rho = np.array([8.331, 8.14, 8.14, 7.94, 7.93, 8.28, 8.28, 8.83, 8.83, 9.21, 9.21])
    K = np.array([174.1, 232., 212., 277., 250., 326., 300., 416., 390., 484., 460.])
    
    rho0 = np.array([8.33,  9.06896049,  9.06896049,  9.40527703,  9.40527703,  9.7393202 ,
        9.7393202 , 10.28553569, 10.28553569, 10.65321301, 10.65321301])
    K0 = np.array([174., 260.56106042, 260.56106042, 307.38650002, 307.38650002,
       358.70550151, 358.70550151, 453.55646214, 453.55646214,
       525.48401269, 525.48401269])
    
    delta_rho = rho-rho0
    delta_K = K-K0
    print ('delta_rho =', delta_rho)
    data = {'p':pres, 'x':X, 'px':pres*X, 'drho':delta_rho}
    dataf = pd.DataFrame(data)
    dataf.to_csv()
    
    x = dataf[['p', 'x']].values.reshape(-1,2)
    
    z = dataf['drho']
    
    lg = LinearRegression()
    
    model = lg.fit(x, z)
                    
    print ('\n', model.intercept_)
    
    for co in model.coef_:
        print (co)    
    
    
    
def fit_data():
    """Fit XH2O as function of pressure (GPa), temperature (K) and magnesium
    content (mol%) for hydrated Olivine
    """
    
    #collect data axes
    
    pres = []
    temp =[]
    XH2O = []
    Mg = []
    
    models = []
    
    #iterate over phases    
    for i in range(len(data_saturation_Mg2SiO4)):
        pres.append([])
        temp.append([])
        XH2O.append([])
        Mg.append([])
        
        dat = data_saturation_Mg2SiO4[i]
        
        #iterate over different data sources
        #ommit last entry as it contains the source reference
        for j in range(len(dat)):
            for k in range(len(dat[j])-1):
                #convert temperature to K
                temp[i].append(dat[j][k][0]+273.)
                
                #pressure in GPa
                pres[i].append(dat[j][k][1])
                
                #Convert Mg# to Fe# in mol fraction
                Mg[i].append(1-dat[j][k][2]/100)

                #Water content in wt fraction
                XH2O[i].append(dat[j][k][3]/100)
        
    print ('the fitted coefficients are:')
    
    for i in range(len(data_saturation_Mg2SiO4)):
    
        data = {'pres': np.array(pres[i]),
                'temp': np.array(temp[i]),
                'temp_pres': np.array(pres[i])*np.array(temp[i]),
                'XH2O': np.array(XH2O[i]),
                'Mg_temp':np.array(Mg[i])*np.array(temp[i]),
                'Mg_pres':np.array(Mg[i])*np.array(pres[i]),
                'Mg_pres_temp':np.array(Mg[i])*np.array(temp[i])*np.array(pres[i]),
                'Mg': np.array(Mg[i])}
        
        dataf = pd.DataFrame(data)
        dataf.to_csv()
        
        thresh = 1
        
        if i < thresh:
            x = dataf[['Mg', 
                       'pres', 
                       'Mg_pres', 
                       'temp', 
                       'Mg_temp',
                       'temp_pres',
                       'Mg_pres_temp']].values.reshape(-1, 7)               

        else:
            x = dataf[['Mg', 
                       'pres',  
                       'temp', 
                       'temp_pres']].values.reshape(-1, 4)  
            
        z = dataf['XH2O']
        
        lg = LinearRegression()
        
        model = lg.fit(x, z)
                        
        print ('\n', model.intercept_)
        
        for co in model.coef_:
            print (co)
        
        
        if i < thresh:
            models.append([[[model.intercept_, model.coef_[0]], 
                 [model.coef_[1], model.coef_[2]]], 
                [[model.coef_[3], model.coef_[4]], 
                 [model.coef_[5], 0.]]])


        else:
            models.append([[[model.intercept_, model.coef_[0]], 
                 [model.coef_[1], 0.]], 
                [[model.coef_[2], 0.], 
                 [model.coef_[3], 0.]]])            
            

            
    params = ['temp', 'pres', 'Fe']
        
    for i in range(2):
        for j in range(2):
            for k in range(2):
                print ('')
                if i == 1:
                    print (params[0])
                if j == 1:
                    print (params[1])
                    
                if k == 1:
                    print (params[2])
            
    '''
    for i in range(len(data_saturation_Mg2SiO4)):
        data = {'pres': np.array(pres[i]),
                'temp': np.array(temp[i]),
                'temp_pres': np.array(pres[i])*np.array(temp[i]),
                'XH2O': np.array(XH2O[i]),
                'Mg': np.array(Mg[i])}
                #'pres_Mg': np.array(pres[i])*np.array(Mg[i]),
                #'temp_Mg': np.array(temp[i])*np.array(Mg[i]),
                #'pres_temp_Mg': np.array(pres[i])**np.array(temp[i])*np.array(Mg[i])}
        
        dataf = pd.DataFrame(data)
        dataf.to_csv()
        
        x = dataf[['pres', 
                   'temp', 
                   'temp_pres', 
                   'Mg']].values.reshape(-1, 4)
        
        z = dataf['XH2O']
        
        lg = LinearRegression()
        
        model = lg.fit(x, z)
        
        models.append(model)
                
        print ('\n', model.intercept_)
        for co in model.coef_:
            print (co)
    '''
    return models


def fit(P, T, Mg, mod):
    c = mod.coef_
    return mod.intercept_ + P*c[0] + T*c[1] + P*T*c[2] + Mg*c[3]


def rho0(X_H2O, Fe_number, ph):
    """Compute ambient density in kg m^-3 of hydrated Olivine containing 
    variable amounts of iron.
    
    X_H2O = water content in wt%
    Fe_number = [Fe]/([Mg]+[Fe]) iron number in mol%
    ph = region in phase diagram (0=alpha, 1=beta, 2=gamma)
    
    Coefficients to experimental data are taken from Mao 2016
    """
    c1, c2, c3 = rho0_coeffs_Mg2SiO4[ph]
    
    return c1 + c2*Fe_number + c3*X_H2O


def K0(X_H2O, Fe_number, ph):
    """Compute isothermal bulk modulus at ambient conditions in GPa for 
    hydrated Olivine containing variable amounts of iron.    
    
    X_H2O = water content in wt%
    Fe_number = [Fe]/([Mg]+[Fe]) iron number in mol%
    ph = region in phase diagram (0=alpha, 1=beta, 2=gamma)
    
    Coefficients to experimental data are taken from Mao 2016
    Note that in Mao 2016 coefficients were fitted to the adiabatic bulk
    modulus. We assume that the slopes for the isothermal bulk modulus are the
    same and use these coefficients and use values for the isothermal bulk
    modulus at ambient conditions as the intercept (first coefficient) instead
    of the values for the adiabatic bulk modulus.
    """

    c1, c2, c3, c4 = K0_coeffs_Mg2SiO4[ph]
    
    return c1 + c2*Fe_number + c3*X_H2O + c4*Fe_number*X_H2O


def G0(X_H2O=0.0, Fe_number=0.0, ph=None):
    """Compute shear modulus at ambient conditions in GPa for
    hydrated Olivine containing variable amounts of iron.
    
    X_H2O = water content in wt%
    Fe_number = [Fe]/([Mg]+[Fe]) iron number in mol%
    ph = region in phase diagram (0=alpha, 1=beta, 2=gamma)
    
    Coefficients to experimental data are taken from Mao 2016
    """
    c1, c2, c3, c4 = G0_coeffs_Mg2SiO4[ph]
    
    return c1 + c2*Fe_number + c3*X_H2O + c4*Fe_number*X_H2O


def xi(eta, m1, m2):
    """Compute mole fraction from mass fraction of species 2 where species
    1 has a molar abundance of (1-xi)
    """
    return m1*eta/(eta*m1-eta*m2+m2)


def eta(xi, m1, m2):
    """Compute mass fraction from mole fraction of species 2 where species
    1 has a weight abundance of (1-eta)
    """
    return xi*m2/((1.-xi)*m1 + xi*m2)


def X_H2O_sat(P=None, T=None, Fe_number=0., phase=None):
    """Compute saturation water content in wt% of Mg2SiO4 at given pressure,
    temperature and iron content.
    
    P(Pa), T(K), Fe_number(mol%)
    
    Default is Fe_number=0.0
    
    """    
    #convert temperature and pressure in GPa and °C as this is what is used
    #in the fit
    
    try:
        P = P*1.0e-9
    
    except TypeError:
        print ('WARNING: invalid pressure input given:', P)
    
    #determine phase region
    #note that phaseCheck takes P in Pa and T in K
    if phase == None:
        phase = phaseCheck(P=P*1.0e9, T=T, Fe_number=Fe_number)
        
    else:
        pass

    #use the fitted coefficients for the phase at hand
    c = saturation_coeffs_Mg2SiO4[phase]
    
    #return saturation water content and avoid negative values as outputs
    #from the simple linear relation
    #fit has been performed in GPa, K, mol fraction and wt fraction
    XH2O =  0.
    
    for i in range(2):
        for j in range(2):
            for k in range(2):
                XH2O += T**i*P**j*(Fe_number/100)**k*c[i][j][k]
    
    return np.maximum(XH2O*100, 0.001)
        
def XH2O_Stv(P=None, T=None, AlSi=0.):
    c = [-3.765, 6.583e-2, -8.916e-4, -6.454e2, 3.544e1, 1.015e-1]
    XH2O = c[0] + c[1]*P*1.0e-9 + c[2]*T + c[3]*AlSi + c[4]*AlSi*P*1.0e-9 + c[5]*AlSi*T
    
    return 10.**XH2O


def fit_MgO(pi=9, pf=12, N=10, temps=[300], deg=3, log=True):
    import eos
    from scipy.optimize import curve_fit
    
    MgO_dens = np.zeros([len(temps), N])
    MgO2H2_dens = np.zeros([len(temps), N])
    pres = np.logspace(8, 12, 10)
    diffs = np.zeros([len(temps), N])
    legend_list = []
    plot_list = []
    xxx = np.logspace(pi, pf, N)
    
    x = []
    y = []
    z = []
    
    for t in range(len(temps)):
        #compute data for fit
        print ('processing temperature', t,': ',temps[t],' K')
        for p in range(len(pres)):
            MgO_dens[t][p] = eos.Compute(what='dens', ll=5, T=temps[t], 
                    P=pres[p])[0]
            
            MgO2H2_dens[t][p] = eos.Compute(what='dens', ll=11, T=temps[t], 
                       P=pres[p])[0]
            
            x.append(pres[p])
            y.append(temps[t])
            diff = (MgO_dens[t][p]-MgO2H2_dens[t][p])/MgO_dens[t][p]
            z.append(diff)

        #compute data for plot
        for p in range(len(xxx)):
            MgO_dens[t][p] = eos.Compute(what='dens', ll=5, T=temps[t], 
                    P=xxx[p])[0]
            
            MgO2H2_dens[t][p] = eos.Compute(what='dens', ll=11, T=temps[t], 
                       P=xxx[p])[0]
            

            diffs[t][p] = (MgO_dens[t][p]-MgO2H2_dens[t][p])/MgO_dens[t][p]
        
    print (diffs[0])
    
    
    def fit_model(x, y, model):
        coefs = model.coef_
        inter = model.intercept_
        p=x
        
        for pp in range(len(p)):
            p[pp] = max(p[pp], 8)
        
        return p*coefs[0] + p**2*coefs[1] + p**3*coefs[2] + p**4*coefs[3] + \
            y*coefs[4] + p*y*coefs[5] + p**2*y*coefs[6] + p**3*y*coefs[7] + \
            inter
        
    #fit in log(P)
    x = np.log10(np.array(x))
    y = np.array(y)
    z = np.array(z)
    
    data = {'x': x, 'x2':x**2, 'x3':x**3, 'x4': x**4, 'y': y, 
            'xy': x*y, 'x2y': x**2*y, 'x3y': x**3*y, 'z':z}
    dataf = pd.DataFrame(data)

    dataf.to_csv()

    xx = dataf[['x', 'x2', 'x3', 'x4', 'y', 'xy', 'x2y', 'x3y']].values.reshape(-1, 8)
    zz = dataf['z']
    
    lg = LinearRegression()
    
    model = lg.fit(xx, zz)
    
    xx = np.log10(pres)

    
    fig, ax = plt.subplots()
    if log:
        ax.set_xscale('log')
        
    for i in range(len(temps)):

        ax.scatter(xxx, diffs[i]*100, marker='s', color = color_list[i])
        
        pl, = ax.plot(xxx, 100*fit_model(np.log10(xxx), temps[i], model), 
                color = color_list[i], linestyle = '--')
        
        legend_list.append(str(temps[i])+' K')
        plot_list.append(pl)
        
    legend=ax.legend(plot_list, legend_list, loc=1, 
                      fontsize=10, framealpha = plot_params['legendalpha'],
                      )
    ax.tick_params(which='both', direction='in')
    ax.add_artist(legend)
    ax.set_xlim(10**pi, 10**pf)
    ax.set_ylim(0, 60)
    ax.set_xlabel(r'$Pressure \ [Pa]$')
    ax.set_ylabel(r'$ \Delta_{max} \ [ $'+'%'+r'$]$')
    coeffs = [c for c in model.coef_] 
    coeffs.append(model.intercept_)
    return coeffs


def delta_rho_max_MgO(P=None, T=None):
    """Computes relative density change at maximum water content in MgO as
    function of temperature in K and pressure in Pa. The fit coefficients from 
    the fit_MgO() function are used.
    """
    
    coefs = delta_rho_MgO_coeffs
    
    #for P < 1.0e8 Pa use constant value at P = 1.0e8
    x = max(np.log10(P), 8)
    y = T
    
    return (x*coefs[0] + x**2*coefs[1] + x**3*coefs[2] + x**4*coefs[3] + \
        y*coefs[4] + x*y*coefs[5] + x**2*y*coefs[6] + x**3*y*coefs[7] + \
        coefs[8])*100
    
            
def delta_rho(P=None, T=None, ll=1, X_H2O=0.):
    """Computes relative density change upon hydration for a given material as
    a function of water content in wt%, pressure in Pa and temperature in K.
    The value that is returned multiplied by the anhydrous density gives
    the hydrous density of the material.
    
    Input is:
        
        material ll (int)
        pressure P (Pa)
        temperature T (K)
        water content X_H2O (wt%)
        
    Output is:
        
        relative density change (0.0 - 1.0)
    """
    if ll == 5:
        #maximum amount of wt% H2O in MgO if all MgO -> Mg(OH)2
        X_H2O_max = 30.089
        if X_H2O > 0. and X_H2O <= X_H2O_max:
            return (1.-delta_rho_max_MgO(P=P, T=T)*X_H2O/X_H2O_max)
        
        elif X_H2O > X_H2O_max:
            print ('WARNING: maximum hydration level exceeded')
        
        else:
            return 1.
    
    else:
        return 1.


def Sievert(T, P_partial):
    P0 = 1.0e5
    xi_H = np.sqrt(P_partial/P0)*np.exp((-31.8e3-T*38.1)/(T*Rgas))
    return xi_H


def PlotRelDev(temps=[1000, 1500, 2000, 2500, 3000], P_min=0, P_max=30,
               N=25):
    import eos
    
    temps = np.array(temps)
    
    pres = np.linspace(P_min, P_max, N)
    
    dens = np.zeros([len(temps), N])

    for i in range(len(temps)):
        for j in range(N):
            dhyd = eos.Compute(ll=1, T=temps[i], P=pres[j]*1.0e9, saturation=True)[0]
            ddry = eos.Compute(ll=1, T=temps[i], P=pres[j]*1.0e9, saturation=False)[0]
            dens[i][j] = abs(dhyd-ddry)/ddry
    
    
    fig, ax = plt.subplots()
    
    for d in range(len(dens)):
        
        ax.semilogy(pres, dens[d], label='T = '+str(temps[d])+' K')
        
    ax.legend()
    
    ax.set_xlabel(r'$P \ [\rmGPa]$')
    ax.set_ylabel(r'$(\rho_{hyd} - \rho_{dry})/\rho_{dry}$')    


def PlotWaterContent(Tmin = 300, Tmax=2000, Pmin=0, Pmax=25, res = 4,
                     cmap = 'inferno', log=False, Fe_number=0.0, **kwargs):
    
    N = 2**res
    x = np.linspace(Tmin, Tmax, N)
    y = np.linspace(Pmin, Pmax, N)
    
    xx, yy = np.meshgrid(x, y)
    
    zz = np.zeros([N, N])
    
    newcmap = cmap
    
    
    #Customize color map
    cmap = cm.get_cmap(cmap, 256)
    newcmap = cmap(np.linspace(0, 1, 256))
    black = np.array([0., 0., 0., 1.])
    dark = np.array([0.1, 0.05, 0.0, 1.])
    middark = np.array([0.2, 0.075, 0.0, 1.])
    lightdark = np.array([0.3, 0.1, 0.0, 1.])
    light = np.array([.375, .125, .0, 1.])
    redlight = np.array([.45, .175, .01, 1.])
    end = np.array([.5, .225, .01, 1.])
    endend = np.array([.53, .25, .02, 1.])
    newcmap[:1, :] = black
    newcmap[1:2, :] = dark
    newcmap[2:3, :] = middark
    newcmap[3:4, :] = lightdark
    newcmap[4:5, :] = light
    newcmap[5:6, :] = redlight
    newcmap[6:7, :] = end
    newcmap[7:8, :] = endend
    newcmap = ListedColormap(newcmap)
    
    
    for i in range(N):
        for j in range(N):
            
            ph = phaseCheck(T=xx[i][j], P=yy[i][j]*1.0e9)
            zz[i][j] = X_H2O_sat(T=xx[i][j], P=yy[i][j]*1.0e9, phase = ph,
              Fe_number=Fe_number)
    
    #plot = plotTools.ColorMap(**kwargs)
    
    #plot.Plot(Z=zz)

    
    fig, ax = plt.subplots()

    geotherm_color = (.25, .5, 1)

    ax.plot(T_mantle_geotherm, P_mantle_geotherm, linewidth=2,
            linestyle ='--', color = geotherm_color)
        
    if log:
        im = ax.imshow(np.log10(zz), origin='lower', 
                   extent=[Tmin, Tmax, Pmin, Pmax],
              aspect='auto', cmap = newcmap,
              vmin=-3, vmax=1)
    
    else:
        im = ax.imshow(zz, origin='lower', 
                   extent=[Tmin, Tmax, Pmin, Pmax],
              aspect='auto', cmap = newcmap,
              vmin=0, vmax=10)
    
    
    ax.set_xlabel(r'$T \ \rm [K]$')
    ax.set_ylabel(r'$P \ \rm [GPa]$')
    cbar = fig.colorbar(im, ax=ax)
    #cbar.ax.set_yticklabels(['0','1','2','>3'])
    cbar.set_label(r'$X_{H2O} \ \rm [wt \%]$', rotation=90, labelpad=10)

    ax.text(600, 10, 'mantle geotherm \n(Stixrude 2009)',
            color = geotherm_color, rotation=50, fontsize = 8)
    ax.text(1400, 2, 'Fe# = '+str(int(Fe_number))+' mol%',
            color = 'white', fontsize=12)


def Compare(N=10, temps=[1000], P_start=1.0e7, P_end=30.0e9, 
            background_col=background_color, fnts1=8, fnts2=10, fnts3=12,
            param = 'dens',
            **kwargs):
    import eos
    
    legend_digits = 2
    dens_scale = 1.0e-3
    zorder_plot = 1
    zorder_grid = 0
    
    m1=molar_mass_list[1]
    m2=molar_mass_list[0]
    
    pres=np.linspace(P_start, P_end, N)
    dens=np.empty([N])
    #water mass fractions
    massfractions = np.array([.0, .01, 0.02, 0.03, 0.04, 0.05, 0.06])
    
    #water mole fractions
    molefractions = xi(massfractions, m1, m2)
    Fe_numbers = np.array([25., 30., 35., 40., 45., 50.])
    
    panel_labels = [['a', 'b'], ['c', 'd'], ['e', 'f'], ['g', 'h'], ['i', 'j']]
    
    fnts=10
    lwd=2

    axis_limits = [[[], []], [[], [], [], []]]

    majorlocatorx = [[[], []], [[], [], [], []]]
                        
    majorlocatory = [[[], []], [[], [], [], []]]

    minorlocatorx = [[[], []], [[], [], [], []]]
    
    minorlocatory = [[[], []], [[], [], [], []]]

    axis_labels = [[[], []], [[], [], [], []]]

    plot_titles = [[[], []], [[], [], [], []]]
    
    
    if param == 'dens':
        x1 = [0, P_end*1.0e-9]
        y1 = [2.5, 4.0]
        scale = 1.0e-3

        x2 = [0, P_end*1.0e-9]
        y2 = [3.0, 4.5]
        
        dxmaj = 5.
        dymaj = .5
        
        dxmin = 1.
        dymin = .1
        
        yaxis_label = r'$\rmDensity \ [10^3 \ kg \ m^{-3}]$'
        
    elif param == 'KT':
        x1 = [0, P_end*1.0e-9]
        y1 = [100, 300]
        scale = 1.0e-9

        x2 = [0, P_end*1.0e-9]
        y2 = [100, 300]
        
        dxmaj = 5.
        dymaj = 50.
        
        dxmin = 1.
        dymin = 25.
        
        yaxis_label = r'$K_T \rm \ [GPa]$'
        
    elif param == 'alpha':
        x1 = [0, P_end*1.0e-9]
        y1 = [0, 20]
        scale = 1.0e5

        x2 = [0, P_end*1.0e-9]
        y2 = [0, 12]
        
        dxmaj = 5.
        dymaj = 5.
        
        dxmin = 1.
        dymin = 1.
        
        yaxis_label = r'$\alpha_{th} \rm \ [10^{-5} \ K^{-1}]$'
        
    elif param == 'adgrad':
        x1 = [0, P_end*1.0e-9]
        y1 = [0.025, .1]
        scale = 1

        x2 = [0, P_end*1.0e-9]
        y2 = [0.025, .1]
        
        dxmaj = 5.
        dymaj = .025
        
        dxmin = 1.
        dymin = .005
        
        yaxis_label = r'$\nabla_{ad}$'
                
    
    for col in range(len(temps)):
        #Prepare plot params for primary plot
        axis_limits[0][0].append([x1, y1])
        axis_limits[0][1].append([x2, y2])
        
        majorlocatorx[0][0].append(dxmaj)
        majorlocatorx[0][1].append(dxmaj)

        majorlocatory[0][0].append(dymaj)
        majorlocatory[0][1].append(dymaj)

        minorlocatorx[0][0].append(dxmin)
        minorlocatorx[0][1].append(dxmin)

        minorlocatory[0][0].append(dymin)
        minorlocatory[0][1].append(dymin)
        
        plot_titles[0][0].append(str(temps[col])+' K')
        plot_titles[0][1].append(str(temps[col])+' K')
        
        axis_labels[0][0].append(['Pressure [GPa]', 
                         yaxis_label])
        axis_labels[0][1].append([r'$\rm{Pressure \ [GPa]}$', 
                         yaxis_label])
        
        #Prepare plot params for relative deviation plot
        axis_limits[1][0].append([x1, [-10, 2]])
        axis_limits[1][1].append([x2, [-30, 2]])
        axis_limits[1][2].append([x2, [0, 400]])
        axis_limits[1][3].append([x2, [-10, 20]])
        
        majorlocatorx[1][0].append(dxmaj)
        majorlocatorx[1][1].append(dxmaj)
        majorlocatorx[1][2].append(dxmaj)
        majorlocatorx[1][3].append(dxmaj)
        
        majorlocatory[1][0].append(2)
        majorlocatory[1][1].append(5)
        majorlocatory[1][2].append(100)
        majorlocatory[1][3].append(5)
        
        minorlocatorx[1][0].append(dxmin)
        minorlocatorx[1][1].append(dxmin)
        minorlocatorx[1][2].append(dxmin)
        minorlocatorx[1][3].append(dxmin)
        
        minorlocatory[1][0].append(1)
        minorlocatory[1][1].append(1)
        minorlocatory[1][2].append(25)
        minorlocatory[1][3].append(1)
        
        plot_titles[1][0].append(str(temps[col])+' K')
        plot_titles[1][1].append(str(temps[col])+' K')
        plot_titles[1][2].append(str(temps[col])+' K')
        plot_titles[1][3].append(str(temps[col])+' K')

        axis_labels[1][0].append(['Pressure [GPa]', 
                                    r'$\delta \rho \ \rm [\%]$'])
        axis_labels[1][1].append([r'$\rm{Pressure \ [GPa]}$', 
                                    r'$\delta K_T \ \rm [\%] $'])
        axis_labels[1][2].append([r'$\rm{Pressure \ [GPa]}$', 
                                    r'$\delta \alpha_{th} \ \rm [\%] $'])    
        axis_labels[1][3].append([r'$\rm{Pressure \ [GPa]}$', 
                                    r'$\delta \nabla_{ad} \ \rm [\%] $'])                

    plots = []
    rows = [2, 4]
    hspaces = [.05, .1]
    
    
    for i in range(2):
        pl = plotTools.Plot(col=len(temps), row=rows[i], 
                            axis_limits=axis_limits[i], 
                            majorlocatorx=majorlocatorx[i],
                            majorlocatory=majorlocatory[i], 
                            minorlocatorx=minorlocatorx[i], 
                            minorlocatory=minorlocatory[i], 
                            axis_labels = axis_labels[i],
                            plot_titles=plot_titles[i],
                            sharex=[[True]], 
                            sharey=[True for t in range(len(temps))],
                            hspace=hspaces[i], 
                            wspace=.1,
                            majorlabelsize=fnts2, 
                            title_pad=30,
                            axislabelsize=fnts3,
                            figsize=(10, 5),
                            **kwargs)
    
        plots.append(pl)
    
    for t in range(len(temps)):
        plot=plots[0]
        T = temps[t]
        ax = [plot.ax[0][t], plot.ax[1][t]]
        
        #ax[0].text(19., 3.81, str(T)+' K', fontsize=fnts3)
        #ax[1].text(19., 4.31, str(T)+' K', fontsize=fnts3)
        #ax[0].text(22., 2.6, '{}'.format(panel_labels[t][0]), fontsize=fnts3,
         # weight='bold')
        #ax[1].text(22., 2.6 + 0.5, '{}'.format(panel_labels[t][1]), fontsize=fnts3,
         # weight='bold')
         
        plot_list2=[]
        plot_list3=[]
        legend_list1=[]
        legend_list3=[]
    
        #plot different water mass fractions
        for c in range(len(massfractions)):
            xi2=xi(massfractions[c], m1, m2)
            cmax = len(massfractions)-1
            color=(.2+.5*(1.-c/cmax), .4, .25+.75*(c/cmax)**2)#color_list[-c-7]
            
            dens0 = np.empty([N])
            dens_max_hyd = np.empty([N])
            dens0_max_hyd = np.empty([N])
                    
            for i in range(N):
                #Compute dens at given hydration using linear hydration model
                d0, T, P, dPdrho, phase, X_H2O, alpha = eos.Compute(T=T, 
                                P=pres[i], ll=1, X_H2O = massfractions[c]*100)
                
                if param == 'dens':
                    dens0[i] = d0
                    
                elif param == 'KT':
                    dens0[i] = d0 * dPdrho
                    
                elif param =='alpha':
                    dens0[i] = alpha
                    
                elif param == 'adgrad':
                    dens0[i] = P/((1+alpha*T)*d0*dPdrho)
                        
            p2,=ax[0].plot(pres*1.0e-9, dens0*scale, color=color, 
                  linestyle='--', linewidth=lwd, zorder=zorder_plot)
            
            if not molefractions[c] == 0.0:
                mole = str(ftool.fancyround(molefractions[c]*100, legend_digits))
    
            else:
                mole = (1+legend_digits)*' '+'0.0'
            
            legend_list1.append(
                                r'$X_{H2O} \rm \ = \ $'+\
                                str(int(round(100*massfractions[c], 
                                          legend_digits)))+r'$\rm\ wt\%$')
            
            plot_list2.append(p2)
        
        max_hyd_list = []
        xi_hyd_list = []
        
        #plot density of water saturated, iron free Olivine
        for i in range(N):
    
            max_hyd = X_H2O_sat(P=pres[i], T=T)
            max_hyd_list.append(max_hyd)
            
            #discard negative hydration
            max_hyd = max(max_hyd, 0.)
            
            #compute water mole fraction from given mass fraction
            xi2=xi(max_hyd/100., m1, m2)
            xi_hyd_list.append(xi2)

            #Compute density with max hydration using bilinear hydration model
            dens, T, P, dPdrho, phase, X_H2O, alpha = eos.Compute(T=T, 
                                            P=pres[i], ll=1, saturation=True)

            if param == 'dens':
                dens0_max_hyd[i] = dens
                
            elif param == 'KT':
                dens0_max_hyd[i] = dens * dPdrho

            elif param == 'alpha':
                dens0_max_hyd[i] = alpha
                
            elif param == 'adgrad':
                dens0_max_hyd[i] = P/((1+alpha*T)*dens*dPdrho)
            
        #plot density of water saturated Olivine for different iron contents
        for c in range(len(Fe_numbers)):
            params_hyd = np.zeros([N, 6])
            params_dry = np.zeros([N, 6])            
            delta = np.zeros([len(Fe_numbers), N, 6])
            dens_dry_Fe = np.zeros([N])
            dens0_max_hyd_Fe = np.zeros([N])
            
            for i in range(N):

                dens_dry, T, P, dPdrho_dry, phase, X_H2O_dry, alpha_dry = eos.Compute(T=T, 
                                        P=pres[i], ll=1, Fe_number=Fe_numbers[c])  

                dens_max_hyd, T, P, dPdrho_max_hyd, phase, X_H2O_max_hyd,  \
                alpha_max_hyd = eos.Compute(T=T, 
                                P=pres[i], ll=1, Fe_number=Fe_numbers[c], 
                                saturation=True)  

                params_hyd[i][:] = dens_max_hyd, dPdrho_max_hyd, X_H2O_max_hyd, \
                                alpha_max_hyd, dens_max_hyd*dPdrho_max_hyd,\
                                P/((1+alpha_max_hyd*T)*dens_max_hyd*dPdrho_max_hyd)
                params_dry[i][:] = dens_dry, dPdrho_dry, X_H2O_dry, alpha_dry,\
                                dens_dry*dPdrho_dry, \
                                P/((1+alpha_dry*T)*dens_dry*dPdrho_dry)
                
                if param == 'dens':
                    dens_dry_Fe[i] = dens_dry
                    dens0_max_hyd_Fe[i] = dens_max_hyd

                elif param == 'KT':
                    dens_dry_Fe[i] = dens_dry*dPdrho_dry
                    dens0_max_hyd_Fe[i] = dens_max_hyd*dPdrho_max_hyd

                elif param == 'alpha':
                    dens_dry_Fe[i] = alpha_dry
                    dens0_max_hyd_Fe[i] = alpha_max_hyd
                    
                elif param == 'adgrad':
                    dens_dry_Fe[i] = P/((1+alpha_dry*T)*dens_dry*dPdrho_dry)
                    dens0_max_hyd_Fe[i] = T/((1+alpha_max_hyd*T)*\
                                            dens_max_hyd*dPdrho_max_hyd)

                max_hyd = X_H2O_sat(P=pres[i], T=T, Fe_number = Fe_numbers[c])
                
                #discard negative hydration
                max_hyd = max(max_hyd, 0.)
                
                #Compute relative difference between dry and hydrous
                delta[c][i] = (params_hyd[i] - params_dry[i])/params_dry[i]
                
            #print ('\n', delta[c])
            #print (np.transpose(delta[c])[0])
            #delta.reshape(len(Fe_numbers), 4, N)       
            #print (delta[0][0])
            
            plots[1].ax[0][t].plot(pres*1.0e-9, 100*np.transpose(delta[c])[0], 
                 color=color_list[-c-1])

            plots[1].ax[1][t].plot(pres*1.0e-9, 100*np.transpose(delta[c])[4], 
                 color=color_list[-c-1])
            
            plots[1].ax[2][t].plot(pres*1.0e-9, 100*np.transpose(delta[c])[3], 
                 color=color_list[-c-1])

            plots[1].ax[3][t].plot(pres*1.0e-9, 100*np.transpose(delta[c])[-1], 
                 color=color_list[-c-1])
        
            if not Fe_numbers[c] < 10.:
                mole = str(ftool.fancyround(Fe_numbers[c], legend_digits))
    
            else:
                mole = legend_digits*' '+str(ftool.fancyround(Fe_numbers[c], legend_digits))
    
            legend_list3.append(r'$\xi_{\rm Fe} = $'+str(int(Fe_numbers[c]))+r'$\rm\ mol\%$')
            
            #Plot iron bearig, anhydrous (Mao 2016)
            p5, = ax[1].plot(pres*1.0e-9, dens_dry_Fe*scale, 
                    color=color_list[-c-1], linestyle='-', zorder=zorder_plot,
                    linewidth=lwd)
            
            #Plot iron bearing hydrous
            p6, = ax[1].plot(pres*1.0e-9, dens0_max_hyd_Fe*scale, 
                    linestyle='--',  color=color_list[-c-1], zorder=zorder_plot,
                    linewidth=lwd)
            
            plot_list3.append(p5)
    
        #Plot water saturation model
        p4,=ax[0].plot(pres*1.0e-9, dens0_max_hyd*scale, color='k', 
              linestyle=':', linewidth=lwd)
    
        if t == 0:
            legend1=ax[0].legend(plot_list2, legend_list1, loc=4, fontsize=fnts1,
                              edgecolor=plot_params['legendedgecolor'],
                              framealpha = plot_params['legendalpha'],
                              facecolor = background_color)
            
            legend2=ax[0].legend([plot_list2[0], p4], 
                              [r'$\mathbf{linear \ hydration}$'+'\n'+ r'$\mathbf{ (Mao \ 2016)}$',
                               r'$\mathbf{saturated}$'+'\n'+\
                               r'$\mathbf{(this \ work)}$'], 
                              loc=2, 
                              fontsize=fnts1,
                              edgecolor=plot_params['legendedgecolor'],
                              framealpha = plot_params['legendalpha'],
                              facecolor = background_color)
            
            legend3=ax[1].legend(plot_list3, legend_list3, loc=4, fontsize=fnts1,
                              edgecolor=plot_params['legendedgecolor'],
                              framealpha = plot_params['legendalpha'],
                              facecolor = background_color)
            
            legend4=ax[1].legend([p5, p6], [r'$\mathbf{anhydrous}$'+\
                      r'$\mathbf{ \ (Mao \ 2016)}$', 
                      r'$\mathbf{saturated \ (this \ work)}$'], loc=2, 
                              fontsize=fnts1, framealpha = plot_params['legendalpha'],
                              edgecolor=plot_params['legendedgecolor'], 
                              facecolor = background_color)
        
            ax[0].add_artist(legend1)
            ax[0].add_artist(legend2)
        
            ax[1].add_artist(legend3)  
            ax[1].add_artist(legend4)              
    
    plots[0].fig.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/rho_hydr_iron_compare.pdf',
           format = 'pdf', bbox_inches='tight')
    return plot


#Here all the data is plotted along side with the linear hydration model
#from the fits
'''
#gather pres data
pres_data = [[], [], []]
temp_data = [[], [], []]
XH2O_data = [[], [], []]
Mg_data = [[], [], []]

#rearange data for plotting
for ph in range(3):
    for i in range(len(data_saturation_Mg2SiO4[ph])):
        for j in range(len(data_saturation_Mg2SiO4[ph][i])-1):
            temp_data[ph].append(data_saturation_Mg2SiO4[ph][i][j][0])
            pres_data[ph].append(data_saturation_Mg2SiO4[ph][i][j][1])
            Mg_data[ph].append(data_saturation_Mg2SiO4[ph][i][j][2])
            XH2O_data[ph].append(data_saturation_Mg2SiO4[ph][i][j][3])
        



mod = []
for i in range(3):
    mo = fit_data()
    
    mod.append(mo)

mod_alpha, mod_beta, mod_gamma = mod



ymax_list = [.5, 5., 4.]
ytick_steps = [6, 6, 5]


pres_list = [[.2, 1.5, 2.5, 5., 6., 6.3, 7.5, 8., 9., 12.], 
             [13.3, 13.5, 14.2, 14.25, 14.7, 15., 15.5, 16., 17., 18.], 
             [19., 20.]]

temp_list = [[1100, 1200, 1400],
             [1100, 1200, 1300, 1400], 
             [1300, 1400, 1500]]

pres_range_list = np.array([[0., 16.],
                   [12., 22.],
                   [17., 23.]])

temp_range_list = np.array([[800, 1800],
                   [0., 3000],
                   [500., 2000]])


fig, ax = plt.subplots(3, 3)

#turn off axis at the right for showing scatter legend
for axx in ax:
    axx[-1].axis('off')

pos = 0.
label_list = []
scatter_list = []
for symb in plot_symbols:
    
    lb = plot_symbols[symb][0]
    fc = plot_symbols[symb][1]
    
    if fc == '':
        facecolor = ''
        
    else:
        facecolor = 'k'
   
    sc = ax[0][2].scatter(0, 0, marker=lb, facecolor=facecolor, color='k')
    
    scatter_list.append(sc)
    label_list.append(symb)
    pos -= 10.

legend=ax[0][2].legend(scatter_list, label_list, loc=2, framealpha=1.,
         frameon=False)

ax[0][2].add_artist(legend)

for i in range(len(pres_list)):
    
    ax[i][0].set_ylabel(r'$wt \% \ H_2O$')
    
    for j in range(len(pres_list[i])):
        p = pres_list[i][j]
        T_alpha_beta = T_alpha_beta_fit(p)
        T_beta_gamma = T_beta_gamma_fit(p)
        
        #define phase range for plotting
        plot_range_temp = T_phase_range(P=p, ph=i)
        
        #Compute saturated water content from bilinear fit
        X_H2O = fit(p, plot_range_temp, 100.,  mod[i])
        
        ax[i][0].plot(plot_range_temp, X_H2O,
          label = str(p)+' GPa', color = color_list[j])
        
        ax[i][0].plot([T_alpha_beta, T_alpha_beta], [0., ymax_list[i]],
          linestyle = '--', color = color_list[j])

        ax[i][0].plot([T_beta_gamma, T_beta_gamma], [0., ymax_list[i]],
          linestyle = ':', color = color_list[j])
        
        ax[i][0].set_ylim(0, ymax_list[i])
        ax[i][0].yaxis.set_ticks(np.linspace(0., ymax_list[i], ytick_steps[i]))
        
        ax[i][0].set_xlim(temp_range_list[i][0], temp_range_list[i][1])
        ax[i][0].tick_params(labelright='off', labelleft='on', top=True, 
          right=True)
        
        
        #extract and plot experimental data for the given phase
        for k in range(len(data_saturation_Mg2SiO4[i])):
            dat = data_saturation_Mg2SiO4[i][k]
            for l in range(len(dat)-1):
                author = dat[-1]
                if dat[l][1] == p and dat[l][2] == 100:
                    try:
                        symb = plot_symbols[author]
                    except KeyError:
                        symb = ['','']
                        
                        
                    if symb[1] == 'col':
                        fc = color_list[j]
                        
                    else:
                        fc = symb[1]
                        
                    ax[i][0].scatter(dat[l][0], dat[l][3], color=color_list[j],
                      marker=symb[0], facecolor = fc)

    
    for j in range(len(temp_list[i])):
        t = temp_list[i][j]
        P_alpha_beta = P_alpha_beta_fit(t)
        P_beta_gamma = P_beta_gamma_fit(t)

        #define phase range for plotting
        plot_range_pres = P_phase_range(T=t, ph=i)
        
        
        ax[i][1].plot(plot_range_pres, 
          fit(plot_range_pres, t, 100., mod[i]), 
          label = str(t)+' °C', color = color_list[j])

        ax[i][1].plot([P_alpha_beta, P_alpha_beta], [0., ymax_list[i]],
          linestyle = '--', color = color_list[j])

        ax[i][1].plot([P_beta_gamma, P_beta_gamma], [0., ymax_list[i]],
          linestyle = ':', color = color_list[j])
        
        ax[i][1].set_ylim(0., ymax_list[i])
        ax[i][1].set_xlim(pres_range_list[i][0], pres_range_list[i][1])
        ax[i][1].yaxis.set_ticks(np.linspace(0., ymax_list[i], ytick_steps[i]))
        
        ax[i][1].tick_params(labelright='on', labelleft='off', top=True, 
          right=True)
        
        #extract and plot experimental data for the given phase
        for k in range(len(data_saturation_Mg2SiO4[i])):
            dat = data_saturation_Mg2SiO4[i][k]
            for l in range(len(dat)-1):
                author = dat[-1]
                if dat[l][0] == t and dat[l][2] == 100:
                    try:
                        symb = plot_symbols[author]
                    except KeyError:
                        symb = ['','']
                        
                    if symb[1] == 'col':
                        fc = color_list[j]
                        
                    else:
                        fc = symb[1]
                        
                    ax[i][1].scatter(dat[l][1], dat[l][3], color=color_list[j],
                      marker=symb[0], facecolor = fc )
        
for axx in ax:
    for axxx in axx:
        axxx.legend()

ax[2][0].set_xlabel(r'$Temperature \ [°C]$')
ax[2][1].set_xlabel(r'$Pressure  \ [GPa]$')
'''

#temperatures ar in °C, pressures in GPa
from matplotlib.lines import Line2D
temp=np.linspace(10, 2500, 2)
pres=P_beta_gamma_fit(temp)

fig, ax = plt.subplots()
ax.tick_params(which='major', direction='in', right='on', top='on',
               size=10)

ax.tick_params(which='minor', direction='in', right='on', top='on')

ax.scatter(T_alpha_beta_data, P_alpha_beta_data, color=(.6, .6, .6),
             linestyle='None', marker='s')

ax.scatter(T_beta_gamma_data, P_beta_gamma_data, color=(.6, .6, .6),
             linestyle='None', marker='s')


trans_color=color_list[7]
phase_color='grey'
ax.set_xlim(10, 2500)
ax.set_xlabel(r'$\rm Temperature \ (°C)$')
ax.set_ylabel(r'$\rm Pressure \  (GPa)$')


#ax.set_title(r'$phase \ diagram \ of \ (Mg)_2SiO_4$')
ax.text(750, 6.5, r'$Forsterit \ (\alpha)$', color=phase_color, 
        rotation=0, fontsize=10)
ax.text(750, 14., r'$Wadsleyite \ (\beta)$', color=phase_color, 
        rotation=0, fontsize=10)
ax.text(750, 21.5, r'$Ringwoodite \ (\gamma)$', color=phase_color, 
        rotation=0, fontsize=10)

#add data points to phase diagram
#loop over phases
for i in range(len(data_saturation_Mg2SiO4)):
    #loop over data sets
    for j in range(len(data_saturation_Mg2SiO4[i])):
        dat = data_saturation_Mg2SiO4[i][j]
        
        #loop over all data points
        for d in range(len(dat)-1):
        
            if dat[d][2] == 100:
                ax.scatter(dat[d][0], dat[d][1], color = color_list[i], 
                           facecolor='', s=30)
            
#add unsaturated data
plot_list = []
label_list = []
for dat in data_unsaturated_Mg2SiO4:
    for i in range(len(dat)-1):
        ph = dat[i][-1]
        if ph > 2:
            ph = ph+1
        
        # Plot only zero Fe content
        if dat[i][2] == 100:
            pl = ax.scatter(dat[i][0], dat[i][1], color=color_list[ph], s=30)
        
    plot_list.append(pl)
    label_list.append(dat[-1])
   
phs = [0, 1, 2, 4, 5, 6]
labs = [r'$\alpha$', r'$\beta$', r'$\gamma$', r'$\alpha + \beta$', r'$\beta + \gamma$' ,r'$\alpha + \gamma$']
legend_elements = [Line2D([0, 1000], [0, 20], 
                          color = color_list[phs[i]],
                          label = "k",
                          ls = '',
                          marker = "o") for i in range(6)]
ax.legend(legend_elements, labs, loc = 4)
     
ax.plot([0, 2000], [P_alpha_beta_fit(0), P_alpha_beta_fit(2000)], color=(.8,.8,.8),
        linestyle = '--', linewidth=2)
ax.plot([0, 2000], [P_beta_gamma_fit(0), P_beta_gamma_fit(2000)], color=(.8,.8,.8),
        linestyle = '--', linewidth=2)
#ax.legend()
ax.set_xlim([600, 2000])
fig.savefig("Mg2SiO2_phase_diagram.pdf", format = "pdf", bbox_inches = "tight")

'''
molar mass of mixture
m=(1-xi)m1+xi mH2O

mass fraction of water
eta=xi mH2O/m = xi mH2O/((1-xi)m1+xi mH2O)
eta((1-xi)m1+xi mH2O) = xi mH2O
eta m1 - eta xi m1 + eta xi mH2O = xi mH2O
eta m1 = xi(eta m1 - eta mH2O +mH2O)
xi = m1*eta/(eta m1 - eta mH2O + mH2O)
'''