#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:40:45 2019

@author: oshah
"""

import inspect
import time
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import sys, re, os
from contextlib import contextmanager
from copy import deepcopy
import pandas as pd
import matplotlib as mpl
from sklearn.linear_model import LinearRegression
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pics import physicalparams
from pics.runparams import color_list
import pickle


def clean_names(folder):
    for root, dirs, filenames in os.walk(folder):
        for dirname in dirs:

            newdirname = dirname.replace(" ", "_")
            newdirname = newdirname.replace("-", "_")
            src = os.path.join(root, dirname)
            dst = os.path.join(root, newdirname)
            print(dirname + "-->", newdirname)
            print("dst =", dst)
            os.rename(src, dst)

    for root, dirs, filenames in os.walk(folder):
        for filename in filenames:

            newfilename = filename.replace(" ", "_")
            newfilename = newfilename.replace("-", "_")
            src = os.path.join(root, filename)
            dst = os.path.join(root, newfilename)
            print(filename + "-->", newfilename)
            print("dst =", dst)
            os.rename(src, dst)


def save_objects(objs, filename):
    with open(filename, "wb") as output:
        pickle.dump(objs, output, pickle.HIGHEST_PROTOCOL)


def load_objects(filename, noisy=False):
    with open(filename, "rb") as f:
        objs = pickle.load(f)

    return objs


def std_corr(arr):
    """Compute corrected standard deviation of array"""
    std = 0.0
    for i in range(len(arr)):
        std += (arr[i] - np.mean(arr)) ** 2

    return np.sqrt(std / (len(arr) - 1))


# return np.sqrt(sum([(arr-np.mean(arr))**2])/(len(arr)-1))


def my_cmap(res=3, start=0, reverse=False):

    border_colors = list(
        reversed(
            [
                (0.7, 0.2, 0.6),
                (0.6, 0.4, 1.0),
                (0.2, 0.4, 1.0),
                (0.2, 0.6, 0.9),
                (0.2, 0.8, 0.8),
                (0.2, 0.8, 0.4),
                (0.6, 0.8, 0.4),
                (0.6, 0.6, 0.2),
                (0.8, 0.4, 0.2),
                (1.0, 0.2, 0.2),
                (1.0, 0.5, 0.5),
                (0.7, 0.7, 0.7),
                (0.5, 0.5, 0.5),
                (0.2, 0.2, 0.2),
                (0.0, 0.0, 0.0),
            ]
        )
    )

    if reverse:
        border_colors = list(reversed(border_colors))

    n_additional_bins = 2**res

    colors = []

    for i in range(len(border_colors) - 1 - start):
        colors.append(border_colors[i + start])
        for r in range(n_additional_bins):
            colors.append(
                np.asarray(border_colors[i + start])
                + (r + 1)
                / (n_additional_bins + 1)
                * (
                    np.asarray(border_colors[i + 1 + start])
                    - np.asarray(border_colors[i + start])
                )
            )

    colors.append(border_colors[-1])

    newcmap = mpl.colors.ListedColormap(colors)

    return newcmap


def extract_col_from_cmap(cmap, val_range, val):
    """Extracts the closest color in cmap for val within a range val_range"""
    N = len(cmap.colors) - 1
    ind = int(N * (val - val_range[0]) / (val_range[1] - val_range[0]))
    return cmap.colors[ind]


def my_dev_cmap(
    N_col=81, border_colors=np.array([color_list[2], [1, 1, 1], color_list[0]])
):

    colors = np.empty([N_col, 3])

    colors[0][:] = border_colors[0][:]

    for i in range(N_col - 1):
        for j in range(3):

            if i < int(N_col / 2):
                ind = [0, 1]

            else:
                ind = [1, 2]

            grad = np.array(
                2 * (border_colors[ind[1]][j] - border_colors[ind[0]][j]) / (N_col - 1)
            )

            colors[i + 1][j] = max(colors[i][j] + grad, 0.0)
            colors[i + 1][j] = min(colors[i + 1][j], 1.0)

    newcmap = mpl.colors.ListedColormap(colors)
    return newcmap


def my_colorbar(mappable, size="5%", pad=0.1, label="a colorbar", fnts=12):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=size, pad=pad)
    cbar = fig.colorbar(mappable, cax=cax)
    cbar.set_label(label, fontsize=fnts)
    cbar.ax.tick_params(labelsize=fnts)
    return cbar


def extract_value_range(data=None, scalings=None, hard_mins=None, hard_maxs=None):
    """
    Extracts the min and max of a data array. The data is assumed to have
    shape [n_rows, n_columns, x, y, n_params]
    """
    N_params = len(data.T)

    if hard_mins == None:
        hard_mins = [np.nan for i in range(N_params)]

    if hard_maxs == None:
        hard_maxs = [np.nan for i in range(N_params)]

    if scalings == None:
        scalings = [1.0 for i in range(N_params)]

    vmins = np.empty([N_params])
    vmaxs = np.empty([N_params])
    print("shape =", np.shape(data))

    for i in range(N_params):

        vmins[i] = np.nanmin(data.T[i]) * scalings[i]
        vmaxs[i] = np.nanmax(data.T[i]) * scalings[i]

        extrem = max(abs(vmins[i]), abs(vmaxs[i]))

        if vmins[i] / abs(vmins[i]) != vmaxs[i] / abs(vmaxs[i]):
            vmins[i] = -extrem
            vmaxs[i] = extrem

        else:
            pass

    return vmins, vmaxs


def read_data_file(file):
    data = []

    with open(file, "r", encoding="utf-8") as infile:
        lines = infile.readlines()
        print(lines[0])
        # Get shape of the original data array
        for c in range(len(lines[0])):
            char = lines[0][c]

            if char == "(":
                shape = lines[0][c + 1 :]
                shape = shape.split(",")
                shape[-1] = shape[-1][0:-2]
                # print ("shape = ", shape)
                for s in range(len(shape)):
                    shape[s] = int(shape[s])

        for line in lines:
            if line.startswith("#"):
                pass

            else:
                data.append(float(line))

    return np.asarray(data), shape


def label_to_param(label=None, types=[], labels=[]):
    params = []

    for i in range(len(labels)):
        pattern = labels[i]
        type = types[i]

        remainder = label[label.index(pattern) + len(pattern) + 1 :]
        r = remainder.split("_")
        print(r)
        if type == "float":

            exp = len(r[1])

            p = float(r[0]) + float(r[1]) * 10 ** (-exp)

        elif type == "int":
            p = int(r[0])

        params.append(p)

    return params


def param_to_label(params=[], types=[], labels=[], rounds=None):
    """


    Parameters
    ----------
    params : TYPE, optional
        Array containing the parameter values for which the label is to be
        generated. The default is [].
    types : TYPE, optional
        Array containing the types of the parameter values. The default is [].
    labels : TYPE, optional
        Array containing the labels of the individual parameters. The default is [].
    rounds : TYPE, optional
        Array containing the number of digits to which each parameter is to be
        rounded for the label. The default is None.

    Returns
    -------
    final_label : TYPE
        DESCRIPTION.

    """

    # By default all parameters are rounded to the first digit
    if rounds == None:
        rounds = [1 for i in range(len(params))]

    final_label = "_"

    for i in range(len(params)):
        type = types[i]
        param = params[i]
        label = labels[i]

        final_label = final_label + label + "_"

        if type == "float":
            # Float numbers are divided into digtis before and after the comma
            p1 = int(round(param, 1))

            p2 = 10 ** rounds[i] * round(param - p1, rounds[i])

            # Compute how many digits need to be printed for p2
            if not p2 == 0.0:
                try:
                    log = int(np.log10(p2))

                except OverflowError:
                    log = rounds[i] - 1

            else:
                log = rounds[i] - 1

            p2 = int(p2)

            # Determine how many zeros there are before the first significant digit
            if log + 1 < rounds[i]:
                if log == 0:
                    zeros = rounds[i] - 1

                else:
                    zeros = min((rounds[i] - (log + 1)), log + 1)

                str_p2 = "0" * zeros + str(p2)

            else:
                str_p2 = str(p2)

            p = str(p1) + "_" + str_p2

        elif type == "int":
            p = str(int(param))

        else:
            print("WARNING: invalid param type")

        if i < len(params) - 1:
            p = p + "_"

        final_label = final_label + p

    return final_label


def gather_data(file):

    data_in = np.load(file)

    data_out = np.empty([12, 14])

    count = 0
    for i in range(2):
        for j in range(2):
            for k in range(3):
                data_out[count][0] = data_in[0][0][i][j][k]  # M
                data_out[count][1] = data_in[1][0][i][j][k]  # R
                data_out[count][2] = data_in[2][0][i][j][k]  # Mg#
                data_out[count][3] = data_in[5][0][i][j][k]  # TS
                data_out[count][4] = data_in[4][0][i][j][k]  # M_Core
                data_out[count][5] = data_in[10][0][i][j][k]  # M_H2O, Core
                data_out[count][6] = (
                    data_in[8][0][i][j][k] - data_in[10][0][i][j][k]
                )  # MH2O, Mantle
                data_out[count][7] = data_in[8][0][i][j][k] * 0.0  # M_Ocean
                data_out[count][8] = data_in[14][0][i][j][k]  # xi_H
                data_out[count][9] = data_in[21][0][i][j][k]  # P_H2
                data_out[count][10] = data_in[13][0][i][j][k]  # R_Ocean
                data_out[count][11] = data_in[12][0][i][j][k]  # PC
                data_out[count][12] = data_in[11][0][i][j][k]  # TC
                data_out[count][13] = 200.0 + i * 500.0

                count += 1

    return data_out


def to_latex_table(data, file_name, digits=3):
    """Takes numpy array as input and converts it to txt file in the standard
    latex table format.
    """
    with open(file_name, "w") as the_file:

        # Create line string from data row

        lines = []

        for i in range(len(data)):
            line_str = ""
            for j in range(len(data[i])):
                val = scinot(data[i][j], digits=digits)
                line_str += str(val)

                if j < len(data[i]) - 1:
                    line_str += " & "

                else:
                    line_str += "\\\\" + "\n"

            lines.append(line_str)
            the_file.write(line_str)

    return lines


def for_recursive(ranges=[], n=None, c=0, iter_list=[], f=None, **kwargs):

    if n == None:
        n = len(ranges)

    if iter_list == []:
        iter_list = [0] * n

    if c == n - 1:
        for iter_list[c] in ranges[c]:
            # print (iter_list)
            f(iter_list=iter_list, **kwargs)

    else:
        for iter_list[c] in ranges[c]:
            for_recursive(
                ranges=ranges, n=n, c=c + 1, iter_list=iter_list, f=f, **kwargs
            )


def row(pnt, order=2):
    """Compute row of coefficient matrix 'M' for 2nd order 2d interpolation"""
    x, y = pnt
    r = []
    for i in range(order + 1):
        for j in range(order + 1):
            r.append(x**i * y**j)

    return r


def gather_pairs(x, y, order=2):
    res = []
    for i in range(len(x)):
        for j in range(len(y)):
            res.append([x[i][j], y[i][j]])
    return res


def construct_matrix(x, y, order=2):
    """Construct coefficient matrix for the points (x[i], y[i])"""
    xx, yy = np.meshgrid(x, y)
    pairs = gather_pairs(xx, yy, order=order)
    print("pairs =", pairs)
    matrix = [row(p, order=order) for p in pairs]

    return matrix


def bilinear_model(x, y, z):
    """x and y are 2d arrays containing the x and y values at which the
    data is given. z is an array containing the 4 results
    """

    # consturct matrix
    matrix = construct_matrix(x, y, order=1)
    coeffs = np.linalg.solve(matrix, z)

    return coeffs


def fit_data(x, y, param, z):
    """Perform bilinear fit of the form:
    z = ax + by + cxy + d*param
    """
    data = {"x": x, "y": y, "xy": x * y, "z": z, "param": param}

    dataf = pd.DataFrame(data)
    dataf.to_csv()

    xx = dataf[["x", "y", "xy", "param"]].values.reshape(-1, 4)
    zz = dataf["z"]

    lg = LinearRegression()

    model = lg.fit(xx, zz)

    # print ('the fitted intercept is:', model.intercept_)
    # print ('the fitted coefficients are:', model.coef_ )

    return model


def fit(mod, x, y, param):
    c = mod.coef_
    return mod.intercept_ + x * c[0] + y * c[1] + x * y * c[2] + param * c[3]


def sectomin(sec=None, digits=1, ms=True):
    """Converts seconds into hours, minutes and seconds and gives the seconds
    with the desired number of digits
    """

    # Compute number of full hours
    hours = int(sec / 3600 + 0.5)

    # Compute remainding seconds
    seconds = sec - hours * 3600.0

    # Compute number of full minutes
    minutes = int(seconds / 60.0 + 0.5)

    # Compute remainding seconds
    seconds = sec - minutes * 60.0

    if not ms:
        # Round seconds to desired precission
        seconds = round(seconds, digits)
        miliseconds = 0

    else:
        # Compute number of full seconds
        seconds = int(seconds + 0.5)

        # Compute remainding miliseconds
        miliseconds = (sec - hours * 3600 - minutes * 60 - seconds) * 1.0e3

        # Round miliseconds to desired precission
        miliseconds = round(miliseconds, digits)

    return [hours, minutes, seconds, miliseconds]


def printTime(sec=None, digits=1, ms=True, where=""):
    hours, minutes, seconds, miliseconds = sectomin(sec=sec, digits=digits, ms=ms)

    if digits == 0:
        if ms:
            miliseconds = int(miliseconds)

        else:
            seconds = int(seconds)

    if ms:
        print(
            "Elapsed time in",
            where,
            ":",
            hours,
            "h",
            minutes,
            "min",
            seconds,
            "sec",
            miliseconds,
            "ms",
        )

    else:
        print("Elapsed time in", where, ":", hours, "h", minutes, "min", seconds, "sec")


def fancyround(num, digits=2):
    """Converts a given value (num) into a rounded value with a fixed number
    of significant digits (digits)
    """

    # compute order of magnitude of the given value
    if not abs(num) < 1e-10 and not np.isnan(num):
        try:
            log = np.log10(num)
            exponent = int(log + 0.5 * (log / abs(log) - 1.0))

        except ValueError:
            exponent = 0

    else:
        exponent = 0

    # round the value to the desired number of significant digits
    dummy_num = num * 10 ** (-exponent)

    dummy_num = round(dummy_num, digits - 1)

    return dummy_num * 10**exponent


def scinot(num, digits=2):
    """Takes value (num) and number of digits (dig) as input and converts the
    number into scientific notation with the corresponding amount of digits
    """
    try:
        log = np.log10(num)
        exponent = int(log + 0.5 * (log / abs(log) - 1.0))

    except ValueError:
        print("num =", num)
        exponent = -10

    except OverflowError:
        exponent = 0

    return str(round(num / 10**exponent, digits - 1)) + "e" + str(exponent)


def print_table(data, cols, wide):
    """Prints formatted data on columns of given width."""
    n, r = divmod(len(data), cols)
    pat = "{{:{}}}".format(wide)
    line = "\n".join(pat * cols for _ in range(n))
    last_line = pat * r
    print(line.format(*data))
    print(last_line.format(*data[n * cols :]))


class Reprinter:
    def __init__(self):
        self.text = ""

    def moveup(self, lines):
        for _ in range(lines):
            sys.stdout.write("\x1b[A")

    def reprint(self, text):
        # Clear previous text by overwritig non-spaces with spaces
        self.moveup(self.text.count("\n"))
        sys.stdout.write(re.sub(r"[^\s]", " ", self.text))

        # Print new text
        lines = min(self.text.count("\n"), text.count("\n"))
        self.moveup(lines)
        sys.stdout.write(text)
        self.text = text


def checkKey(name="", **kwargs):
    value = None
    for key, val in sorted(kwargs.items()):
        if key == name:
            value = val
    return value


def supressStdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def updateKey(f, c0, whicharg, noisy=False, identity="", **kwargs):
    """This function compares the arguments passed to the function f with the
    arguments passed to the functions that calls updateKey. If this function is
    missing one or more parameters, updateKey will add a key : val to the
    dictionary kwargs and return it
    """
    if noisy:
        print("function <updateKey>; id =", identity)

    external = [arg for arg in inspect.getfullargspec(f).args]
    local = [key for key, val in kwargs.items()]

    diff = set(external) - set(local)

    try:
        missingarg = list(diff)[0]

    except IndexError:
        pass

    try:
        if missingarg == whicharg:
            try:
                if noisy:
                    print("old keys =", kwargs)
            except NameError:
                pass

            kwargs[missingarg] = c0

            try:
                if noisy:
                    print("new keys =", kwargs)
            except NameError:
                pass
        else:
            if noisy:
                print("NOTE: Missing argument and given key do not match!")
                print(
                    "      missingarg =",
                    missingarg,
                    ", given key =",
                    whicharg,
                    "\nYou have either passed an ambiguous key",
                    "or forgotten to pass a passive key",
                )
    except UnboundLocalError:
        pass

    for key, val in sorted(kwargs.items()):
        if key == whicharg:
            whicharg = key
            kwargs[key] = c0

    return kwargs


def bisec(
    *args,
    a=0.0,
    b=0.0,
    y=0.0,
    f=None,
    limit=100,
    whicharg="",
    noisy=False,
    identity="",
    plot=False,
    axis=None,
    acc=1.0e-8,
    brash=False,
    log=False,
    **kwargs
):
    """Find the value c at which f(c) = y +- eps*y by bisecting between a and b"""
    #    try:
    #       acc = kwargs['acc']
    #  except KeyError:
    #     acc=1.0e-5

    if noisy:
        print()
        print("function <bisec>; id =", identity)

    # print ('acc in bisec:', acc)
    counter = 0
    c0 = (a + b) / 2.0
    c = c0

    # add missing parameter to kwargs for which the bisection in f is to be carried
    # out
    kwargs = updateKey(f, c0, whicharg, noisy=noisy, **kwargs)

    if noisy:
        print("kwargs in bisec =", kwargs)
    if plot:
        fig, axis = plt.subplots()
        axis.grid(color=(0.85, 0.85, 0.85))

    # no kwargs contains all relevant parameters, one of which will be updated
    # along the ride. args contains possible positional arguments of the function
    # f which have to be passed to bisec manually for each function specifically
    # NOTE: the bisection can thus only be carried out for keyword arguments

    fc = f(*args, **kwargs)
    if fc == y:
        if noisy:
            print("NOTE: fc = y in bisection")

    else:
        while abs(fc - y) / y > acc:
            counter += 1

            # compute new midpoint
            c = (a + b) / 2.0

            # update argument to evaluate function at new midpoint
            kwargs[whicharg] = c
            # evaluate function at new midpoint
            fc = f(*args, **kwargs)

            # update right value for bisection
            kwargs[whicharg] = b
            # evaluate function at new b value
            fb = f(*args, **kwargs)

            if brash:
                print("\nc=", c)
                print("a, b=", a, b)

            if plot:
                if noisy:
                    print("plotting on axis:", axis)
                    print("point: ", c, fc)

                axis.scatter(c, fc, color="k", s=15, marker="x")
                axis.tick_params(right="on", top="on", direction="in")
                axis.set_xlabel(whicharg)
                axis.set_ylabel(str(f))

                if log:
                    axis.set_yscale("log")

            # check if function value at new midpoint differs from target value
            if not y == fc:
                try:
                    # target value between a and c
                    if (fc - y) / abs(fc - y) == (fb - y) / abs(fb - y):
                        b = c

                    # target vale between b and c
                    else:
                        a = c

                    # check if accuracy is reached to update c for putput
                    if abs(fc - y) / y <= acc:
                        # compute new modpoint
                        c = (a + b) / 2.0

                except BaseException as e:
                    if noisy or brash:
                        print("WARNING: an exception has occured: ", e)
            else:
                break

            if counter > limit:
                if noisy or brash:
                    print(
                        "WARNING: bisection step limit exceeded in <",
                        identity,
                        "> after",
                        limit,
                        "steps !",
                    )
                break
    return c


def deriv(
    *args,
    x0=0.0,
    f=None,
    whicharg="x",
    noisy=False,
    identity="",
    plot=False,
    axis=None,
    brash=False,
    rr=10,
    type="tangent",
    order=4,
    log="",
    pnts=30,
    acc=1.0e-6,
    **kwargs
):
    """Differentiate function <f> with respect to x at x=x0 where <f> must be
    replaced by the name of the function. For multivalued functions the
    argument with respect to which the differntiation should be carried out
    can be set by the <whicharg> key (default is 'x', assuming f=<f>(x), where)
    """
    # try:
    #   acc = kwargs['acc']
    # except KeyError:
    #   acc=1.0e-6

    try:
        eps = kwargs["eps"]

    except KeyError:
        eps = 1.0e-2

    time_init = time.time()
    kwargs = updateKey(f, x0, whicharg, noisy=noisy, **kwargs)
    exit_code = 0

    if type == "poly":
        """Fit polynomial of order <order> to the function f in an interval
        [xl, xr] around x0 and compute its derivative at x0"""
        h = x0 * eps
        xl = x0 - h
        xr = x0 + h
        xx = np.linspace(xl, xr, 1000, endpoint=True)
        yy = []

        for i in range(len(xx)):
            kwargs[whicharg] = xx[i]
            yy.append(f(*args, **kwargs))

        z = np.polyfit(xx, yy, order)
        poly = np.poly1d(z)
        dpoly = np.poly1d.deriv(poly)

        time_delta = time.time() - time_init

        if noisy or brash:
            print()
            print("function <deriv>; id =", identity)
            print("type = poly")
            print("updated kwargs =", kwargs)
            print("elapsed time:", round(time_delta * 1.0e3, 3), "ms")
        return dpoly(x0)

    elif type == "tangent":
        """Compute the slope between two points left and right of x0. In each
        step the interval is decreased until the relative change of the value
        of the slope is < acc
        """
        h = x0 * eps
        xl = x0 - h
        xr = x0 + h
        dx = xr - xl
        kwargs[whicharg] = xr
        yr = f(*args, **kwargs)

        kwargs[whicharg] = xl
        yl = f(*args, **kwargs)
        dy = yr - yl
        newslope = dy / dx

        if noisy or brash:
            print()
            print("function <deriv>; id =", identity)
            print("type = tangent")
            print("updated kwargs =", kwargs)
            print(
                "xl =",
                round(xl, rr),
                "\nxr =",
                round(xr, rr),
                "\nyl =",
                round(yl, rr),
                "\nyr =",
                round(yr, rr),
                "\nh =",
                round(h, rr),
                "\neps =",
                eps,
                "\nnewslope =",
                round(newslope, rr),
            )
            print("starting iteration...")

        if plot:
            if axis == None:
                fig, axis = plt.subplots()

            if log == "x":
                xx = np.logspace(np.log10(xl), np.log10(xr), pnts)

            else:
                xx = np.linspace(xl, xr, pnts)

            yy = []
            for x in xx:
                kwargs[whicharg] = x
                yy.append(f(*args, **kwargs))

            if log == "double":
                axis.loglog(xx, yy, zorder=0, label="numerical")

            elif log == "x":
                axis.semilogx(xx, yy, zorder=0, label="numerical")

            elif log == "y":
                axis.semilogy(xx, yy, zorder=0, label="numerical")

            elif log == "":
                axis.plot(xx, yy, zorder=0, label="numerical")

            axis.set_xlabel("x")
            axis.set_ylabel("f(x)")
            scatter_x = []
            scatter_y = []

        # set initial deviation to 100% as before the first iteration it can not
        # be estimated quantitatively
        dev = 1.0
        count = 0
        while dev > acc:
            # -----------------------------------------------------------------
            count += 1
            oldslope = newslope
            # print ('oldslope=', oldslope)
            h = h / 2.0
            xl = x0 - h
            xr = x0 + h
            dx = xr - xl

            kwargs[whicharg] = xr
            yr = f(*args, **kwargs)

            kwargs[whicharg] = xl
            yl = f(*args, **kwargs)
            dy = yr - yl
            # -----------------------------------------------------------------

            if plot:
                scatter_x.append([xl, xr])
                scatter_y.append([yl, yr])

            try:
                newslope = dy / dx

            # if a singularity is encountered use previous slope as result
            except ZeroDivisionError:
                if noisy or brash:
                    print("WARNING: newslope = 0.0 -> division by zero")
                    print("return oldslope =", oldslope)
                exit_code = 1
                return oldslope
                break

            if brash:
                print("#", count)
                print(
                    "xl =",
                    round(xl, rr),
                    "\nxr =",
                    round(xr, rr),
                    "\nyl =",
                    round(yl, rr),
                    "\nyr =",
                    round(yr, rr),
                    "\nh =",
                    round(h, rr),
                    "\neps =",
                    eps,
                    "\nnewslope =",
                    round(newslope, rr),
                )

            if dx == 0.0:
                if noisy or brash:
                    print("dx = 0")

            try:
                dev = abs((oldslope - newslope) / newslope)
                if brash:
                    print("dev =", round(dev, 10), "\n-----")

            except ZeroDivisionError:
                if noisy or brash:
                    print("WARNING: newslope = 0.0 -> division by zero")
                    print("return oldslope =", oldslope)
                exit_code = 1
                return oldslope
                break

        if plot:

            for i in range(len(scatter_x)):
                axis.scatter(scatter_x[i], scatter_y[i], zorder=1)

        if noisy or brash:
            time_delta = time.time() - time_init
            print(
                "iteration terminated after", count, "steps with exit code", exit_code
            )
            print("elapsed time:", round(time_delta * 1.0e3, 3), "ms")

        return newslope


def integrate(
    start=0.0,
    end=None,
    y=[],
    h=None,
    N=1,
    dydx=None,
    order=4,
    plot=False,
    whicharg=None,
    noisy=False,
    brash=False,
    identity="non passed",
    axis=None,
    eps=0.5,
    oldgrad=None,
    **kwargs
):
    """This method integrates a function dydx starting at x=start with initial
    conditions y(x=start) = y from x to x=end. y is an array of n entries
    corresponding to the number of variables for multivalued functions. In this
    case the integration step h will be ommitted if passed and computed as
    h = (end - start)/end.
    Alternatively one can give an integration step h and number of integration
    steps N. The integration is then carried out from x=start to x=x+N*h in N
    steps. In this case the end value will be omitted if passed as argument.
    """

    if whicharg == None:
        print(
            'WARNING: keyword argument "whicharg" has not been passed to',
            "functionTools.integrate()!",
        )
        print(
            "This will most likely lead to a KeyError or TypeError in the",
            "dydx(**kwargs) function call.",
        )

    kwargs = updateKey(dydx, start, whicharg, noisy=noisy, y=y, **kwargs)

    if noisy or brash:
        print()
        print("function <integrate>; id =", identity)
        print("updated kwargs =", kwargs)
        print("old vector =", y)

        if h != None:
            print("integration by step size <h>. end point <end> is omitted!")
            print("step size =", h)
            print("number of steps =", N)

        if end != None:
            print("integration by end point <end>. step size <h> is omitted!")
            print("end point =", end)
            print("number of steps =", N)

        if start == None:
            print("WARNING: <start> is missing")
        if h == None and end == None:
            print("WARNING: <h> or <end> is missing")
    try:
        h = (end - start) / N
    except TypeError:
        pass

    if order == 4:
        # Butcher-table for the 4th order solver
        alist = [
            [0.0, 0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
        blist = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]
        clist = [0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0]

    x = start
    # set up data arrays for subsequent plotting
    if plot:
        x_list = [x]
        y_list = np.zeros((N + 1, len(y)))
        y_list[0][:] = y[:]

    # --------------------------------------------------------------------------
    y = np.array(y)
    yold = deepcopy(y)
    newgrads = None

    for hh in range(N):
        k_list = []
        for i in range(order):
            # print ("i=", i)
            ki = [y[s] for s in range(len(y))]
            # print ("ki = ", ki)
            for j in range(i):
                for s in range(len(y)):
                    ki[s] += h * k_list[j][s] * alist[i][j]
                # print ("ki = ", ki)

            try:
                kwargs[whicharg] = x + h * clist[i]
                kwargs["y"] = ki
                # print ('y=', kwargs['y'])
                # use updated key arguments to compute new gradients:
                # newgrad = dydx(x = x + h*clist[i], y = ki, ...)
                # print ("grads before f=", newgrads)
                newgrads = dydx(oldgrad=oldgrad, **kwargs)
                # print ("grads after f=", newgrads)

            except KeyError:
                newgrads = None
                print("ERROR: missing key in function <integrate>; id =", identity)
                print("kwargs =", kwargs)
                print("key <", whicharg, "> does not exist")

            k_list.append(newgrads)
        # print ("y_out = ", y)
        for i in range(len(y)):
            for j in range(len(blist)):
                y[i] += h * k_list[j][i] * blist[j]

        if identity == "construct shell":
            if y[3] > yold[3]:
                # print ('i knew it...')
                y[3] = yold[3]  # +h*oldgrad[3]

        # update x variable if integration over multiple grid points is performed
        x += h
        # ----------------------------------------------------------------------

        # gather results for plotting
        if plot:
            yy = y[:]
            x_list.append(x)
            y_list[hh + 1][:] = yy[:]

    # rearange result vector to gather all values for a given parameter
    if plot:
        y_list = y_list.transpose()

        if axis == None:
            fig, axis = plt.subplots(len(y))

        for i in range(len(y_list)):
            axis[i].plot(
                x_list,
                y_list[i],
                marker="o",
                linestyle="-",
                markevery=max(int(N / 10), 1),
            )

    if noisy or brash:
        print("new vector =", y)

    return x, y


def deriv_partial(f=None, **kwargs):
    """Compute partial derivative of function f(v1, v2) with respect to v1 or
    v2. In the actual function, the variables can be called anything else, e.g.
    x and y. Then the variable values at which the derivatives have to be
    evaluated must be passed accordingly to the derivderiv function.

    Example: Differentiate the function fct(x, y)=x**2*y**3 with respect to x
    and evaluate it at x=a, y=b where a=2 and b=4.

    In []: deriv_partial(f=fct, y=b, x0=a, whicharg='x')
    Out[]: 255.99999999999858

    Note that the following would yield the same result (c be an number
    different from a):

    In []: deriv_partial(f=fct, x=c y=b, x0=a, whicharg='x')
    Out[]: 255.99999999999858

    The deriv() function that is called will automatically use the argument
    for x0 as input for x to evaluate the function. The additional passing
    of a value for x itself will be ignored by the deriv() function.
    """
    # print ('kwargs in deriv_partial:', kwargs)
    d = deriv(f=f, **kwargs)
    return d


def derivderiv(f=None, **kwargs):
    """Compute partial derivative of f(v1, v2) with respect to v1 and v2.
    In the actual function, the variables can be called anything else, e.g.
    x and y. Then the variable values at which the derivatives have to be
    evaluated must be passed accordingly to the derivderiv function.

    Example: Differentiate the function fct(x, y)=x**2*y**3 with respect to x
    and y and evaluate it at x=a, y=b where a=2 and b=4.

    In []: derivderiv(f=fct, x=a, y=b)
    Out[]: 192.0000003929472
    """
    # gather arguments passed to the function
    all_args = inspect.getfullargspec(f)[0]

    # by convention the first two arguments are for the x and y variables
    # to pass to the function respectively
    # (here denoted as variable 1 and variable 2 to avoid confusion if the
    # actual function arguments are called x and y in the function itself)
    # print ('kwargs in derivderiv:', kwargs)
    # print ('all args in derivderiv:', all_args)
    v1 = kwargs[all_args[0]]
    v2 = kwargs[all_args[1]]

    # first differentiate with respect to variable 1
    def fx(**kwargs):

        d = deriv_partial(f=f, whicharg=all_args[0], x0=v1, **kwargs)
        return d

    # then differentiate this with respect to variable 2
    dd = deriv(f=fx, whicharg=all_args[1], x0=v2, **kwargs)
    return dd


def row1d(x, order=2):
    """Compute row of coefficient matrix 'M' for 2nd order 2d interpolation"""
    r = []
    for i in range(order + 1):
        r.append(x**i)
    return r


def deriv_x_row(pnt, order=2):
    x, y = pnt
    r = []
    for i in range(order + 1):
        for j in range(order + 1):
            r.append(i * x ** (i - 1) * y ** (j))

    return r


def deriv_x_row1d(x, order=2):
    r = []
    for i in range(order + 1):
        r.append(i * x ** (i - 1))

    return r


def deriv_y_row(pnt, order=2):
    x, y = pnt
    r = []
    for i in range(order + 1):
        for j in range(order + 1):
            r.append(x ** (i) * j * y ** (j - 1))

    return r


def deriv_xy_row(pnt, order=2):
    x, y = pnt
    r = []
    for i in range(order + 1):
        for j in range(order + 1):
            r.append((i) * x ** (i - 1) * j * y ** (j - 1))

    return r


"""
def gather_pairs(x, y, order=2, dim=1):
    
    res=[]
    for i in range(len(x)):
        for j in range(len(y)):
            res.append([x[i][j], y[i][j]])
    return res
"""
'''  
def construct_matrix(grid, order=2, dim=1):
    """Construct coefficient matrix for the points (x[i], y[i])
    """
    if dim == 1:
        pairs, = grid
    
    if dim == 2:
        x, y = grid
        xx, yy = np.meshgrid(x, y)
        pairs = gather_pairs(xx, yy, order=order, dim=dim)
#        print ('pairs=', self.pairs)
    matrix = [row(p, order=order, dim=dim) for p in pairs]
    return matrix
'''


def interpolate_deriv(
    x=None,
    y=None,
    x_list=None,
    y_list=None,
    z_list=None,
    derivatives_x=None,
    derivatives_y=None,
    derivatives_xy=None,
    order=3,
    **kwargs
):
    """ """

    print("x_list:", x_list)
    print("y_list:", y_list)
    # construct coefficient matrix
    matrix = []
    b_vector = []

    # the first 4 rows correspond to the four grid point values
    for a in range(2):
        for b in range(2):
            b_vector.append(z_list[a][b])
            matrix.append(row(pnt=[x_list[a], y_list[b]], order=order))
            # print ('fct(a, b)=',fct(x_list[a], y_list[b]))

    # the next four rows correspond to the x-derivatives on six four points
    for a in range(2):
        for b in range(order - 1):
            b_vector.append(derivatives_x[a][b])
            matrix.append(deriv_x_row(pnt=[x_list[a], y_list[b]], order=order))
            # print ('deriv(a,b)=', deriv(x=x_list[a], y=y_list[b], whicharg='x'))

    # the next four rows correspond to the y-derivatives at four grid points
    for a in range(2):
        for b in range(order - 1):
            b_vector.append(derivatives_y[a][b])
            matrix.append(deriv_y_row(pnt=[x_list[a], y_list[b]], order=order))
            # print ('deriv(a,b)=', deriv(x=x_list[a], y=y_list[b], whicharg='y'))

    # the last four rows correspond to the mixed xy-derivatives at four
    # grid points
    for a in range(order - 1):
        for b in range(order - 1):
            b_vector.append(derivatives_xy[a][b])
            matrix.append(deriv_xy_row(pnt=[x_list[a], y_list[b]], order=order))
            # print ('deriv(a,b)=', deriv(x=x_list[a], y=y_list[b], whicharg='y'))

    # print ('matrix=',self.matrix)
    # print ('b vector=', self.b)
    a_vector = np.linalg.solve(matrix, b_vector)
    # print ('a vecotr=', self.a)
    return a_vector


def interpolate_deriv1d(grid=[], data=[], order=2, plot=False, which="left", **kwargs):
    x_list = grid
    values, derivatives = data

    matrix = []
    b_vector = []

    for i in range(2):
        matrix.append(row1d(x_list[i], order=order))
        b_vector.append(data[0][i])

    if order == 2:
        if which == "left":
            matrix.append(deriv_x_row1d(x_list[0], order=order))
            b_vector.append(data[1][0])

        elif which == "right":
            matrix.append(deriv_x_row1d(x_list[1], order=order))
            b_vector.append(data[1][1])

        elif which == "both":
            matrix.append(deriv_x_row1d(x_list[0], order=order))
            b_vector.append(data[1][0])

    else:
        for i in range(order - 1):
            matrix.append(deriv_x_row1d(x_list[i], order=order))
            b_vector.append(data[1][i])

    # compute the coefficients for the 2nd order polynom fit
    a_vector = np.linalg.solve(matrix, b_vector)

    if which == "both":
        matrix[-1] = deriv_x_row1d(x_list[1], order=order)
        b_vector[-1] = data[1][1]
        a_vector_ = np.linalg.solve(matrix, b_vector)

    def func(x, a, order):
        return sum([x**i * a[i] for i in range(order + 1)])

    if plot:
        xx = np.linspace(min(x_list), max(x_list), 25)
        yy = func(xx, a_vector, order=order)

        if which == "both":
            yy_ = func(xx, a_vector_, order=order)
            # plt.plot(xx, yy_, color='g')
            p1 = data[0][0]
            p2 = data[0][1]
            x1 = x_list[0]
            x2 = x_list[1]
            plt.plot(xx, yy_, linestyle="--", color="grey")
            plt.plot(xx, yy, linestyle="--", color="black")

        plt.scatter(x_list[0], data[0][0], color="b")
        plt.scatter(x_list[1], data[0][1], color="g")

        plt.plot(xx, yy, color="b")

    # print ('a:',a_vector)
    if which == "both":
        # print ('a_:',a_vector_)
        return (a_vector + a_vector_) / 2

    else:
        return a_vector


def interpolate(grid=[], data=[], dim=1, order=2, **kwargs):
    """Takes x and y values as input and gives f(x, y) as output using a
    polynomial nth order interpolation scheme. The scheme solves the
    matrix equation M*a=b where 'M' is the coefficient matrix and 'b'
    the solution vector.

    grid: array of dimension dim containing the grid axis for each dimension
        grid = [[x1,... xN], [y1, ... yN]]
    data: array of dimension dim containing the function values
        example for 2d:
        data = [[f(x1, y1), f(x1, y2)], [f(x2, y1), f(x2, y2)]]
    """

    if dim == 1:
        (x_list,) = grid
        if not len(grid[0]) == len(x_list):
            print("WARNING: grid size does not match data matrix")

    elif dim == 2:
        x_list, y_list = grid

        if not len(grid[0]) == len(x_list):
            print("WARNING: x grid size does not match data matrix")

        if not len(grid[1]) == len(y_list):
            print("WARNING: y grid size does not match data matrix")

    # compute coefficient matrix using the xy-grid points
    matrix = construct_matrix(grid, order=order, dim=dim)

    # compute the solutions at the given grid points
    b_vector = []
    if dim == 1:
        b_vector = data

    elif dim == 2:
        for i in range(order + 1):
            for j in range(order + 1):
                # note that here j and i must be interchanged in comparison
                # to the interpolate_deriv methode, otherwise it doesn't work
                b_vector.append(data[j][i])

    # print ('b=',self.b)
    # print ('dim matrix=', len(self.matrix), len(self.matrix[0]))

    # compute the coefficients for the 2nd order polynom fit
    a_vector = np.linalg.solve(matrix, b_vector)

    return a_vector


# ==============================================================================
# Test functions to check algorithms and give example of usage
# This part is normally commented out


# def f1(x1=None, x2=None, **kwargs):
#     return x2 ** 2


# def f2(x1=None, x2=None, x3=None, **kwargs):
#     return np.exp(x3)


# def f3(x1=None, x2=None, **kwargs):
#     return x1 * x2 + x2 ** 2


# def f4(x1=None, **kwargs):
#     return np.sqrt(x1) * np.sin(x1)


# Test some stuff

# fig, ax = plt.subplots()
# w = 5.
# number = .5/np.sqrt(w)*np.sin(w) + np.sqrt(w)*np.cos(w)
# print (round((deriv(x0 = w, f = f4, whicharg = 'x1', type = 'poly',
#                     acc = 1.0e-3, order = 2) - \
#                     number)/number*100, 8),'%')


# #id = 1
# number1 = (bisec(whicharg = 'x1', f = f1, x2 = 2., y = 4., a = 0.,
#                  b = 10., noisy = True, plot = True, axis = ax, identity = '1'))

# #id = 2
# number2 = bisec(whicharg = 'x3', f = f2, x1 = 5., x2 = 3, y = 3., a = 0.,
#                 b = 10., eps = 1.0e-6, noisy = True, identity = '2')

# #id = 3
# number3 = bisec(whicharg = 'x2', f = f3, x1 = 1., y = 3., a = 0.,
#                 identity = '3', b = 10., eps = 1.0e-6, plot = True, axis = ax,
#                 noisy = True)

# #id = 4
# number4 = deriv(f = f4, x0 = 2., x1 = 0., x2 = 1., whicharg = 'x1',
#                 noisy = True, eps = 1.0e-8, identity = '4')


# x_list = np.linspace(0., 5., 25)

# y_list = f3(x2 = x_list, x1 = 1.)

# ax.plot(x_list, y_list)
# plt.show()

# print ()
# print ('number 1 =', number1)
# print ('number 2 =', number2)
# print ('number 3 =', number3)
# print ('number 4 =', number4)
