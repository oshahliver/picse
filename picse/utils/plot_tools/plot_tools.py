import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from picse.runparams import param_colors
from picse.physicalparams import m_earth, r_earth, G
import pandas as pd
import numpy as np
import os


def plot_structure(
    profiles,
    display=True,
    file_name="structure",
    file_path=None,
    image_extension="pdf",
    write_html=False,
    write_image=False,
):
    """Plots the structure profils of a planetary object"""
    if file_path == None:
        file_path = os.getcwd()

    df = pd.DataFrame(
        profiles.T,
        columns=[
            "radius",
            "temperature",
            "pressure",
            "density",
            "mass",
            "gravity",
            "moi",
            "dgrav",
        ],
    )

    ncol = 3
    nrow = 2
    fig = make_subplots(
        rows=nrow,
        cols=ncol,
        shared_xaxes=True,
        horizontal_spacing=0.15,
        vertical_spacing=0.1,
    )

    fig.update_layout(height=500, width=800, showlegend=False)
    scalings = [[1.0, 1.0, 1.0], [1e-9, 1.0 / m_earth, 1.0]]
    params = [["temperature", "density", "gravity"], ["pressure", "mass", "moi"]]
    labels = [["Temperature", "Density", "Gravity"], ["Pressure", "Mass", "Norm. MoI"]]
    units = [
        ["(K)", "(km \ m^{-3})", "(m \ s^{-2})"],
        ["(GPa)", "(\mathit{M_\oplus})", ""],
    ]

    for i in range(2):
        for j in range(3):
            ind = ncol * i + j
            fig.add_trace(
                go.Scatter(
                    x=df["radius"] * 1e-3,
                    y=df[params[i][j]] * scalings[i][j],
                    line=dict(color="rgb{}".format(param_colors[ind])),
                ),
                row=i + 1,
                col=j + 1,
            )

            if i + 2 * j > 0:
                a = str(ind + 1)
            else:
                a = ""

            fig["layout"]["yaxis{}".format(a)]["title"] = r"$\rm {} \ {}$".format(
                labels[i][j], units[i][j]
            )

    for i in range(ncol):
        a = str(i + ncol + 1)
        fig["layout"]["xaxis{}".format(a)]["title"] = r"$\rm {} \ ({})$".format(
            "Radius", "km"
        )

    # Write to interactive html
    if write_html:
        print("writing html to", "{}/{}.html".format(file_path, file_name))
        fig.write_html("{}/{}.html".format(file_path, file_name), include_mathjax="cdn")

    # Save static image
    if write_image:
        fig.write_image("{}/{}.{}".format(file_path, file_name, image_extension))

    # Display in browser
    if display:
        fig.show()


def plot_sample(data, specs={}):
    if not "x" in specs.keys():
        specs.update({"x": "mass"})

    if not "y" in specs.keys():
        specs.update({"y": "radius"})

    if not "z" in specs.keys():
        specs.update({"z": "mg_number"})

    fig = go.Figure(
        go.Scatter(
            x=data[specs["x"]],
            y=data[specs["y"]],
            mode="markers",
            marker=dict(
                size=10,
                colorscale="Viridis",
                color=data[specs["z"]],
                colorbar=dict(title="Z", thickness=20, tickmode="auto", nticks=10),
            ),
        )
    )

    # Set the layout
    fig.update_layout(
        coloraxis_colorbar=dict(yanchor="top", y=1, x=1.15, xanchor="left")
    )

    # Show the plot
    fig.show()


def plot_sample_3d(data, specs={}):
    if not "x" in specs.keys():
        specs.update({"x": "mass"})

    if not "y" in specs.keys():
        specs.update({"y": "radius"})

    if not "z" in specs.keys():
        specs.update({"z": "pres_center"})

    if not "c" in specs.keys():
        specs.update({"c": "mg_number"})

    fig = go.Figure(
        go.Scatter3d(
            x=data[specs["x"]],
            y=data[specs["y"]],
            z=data[specs["z"]],
            mode="markers",
            marker=dict(
                size=5,
                colorscale="Viridis",
                color=data[specs["c"]],
                colorbar=dict(title="c", thickness=20, tickmode="auto", nticks=10),
            ),
        )
    )

    # Set the layout
    fig.update_layout(
        scene=dict(xaxis_title="X", yaxis_title="Y", zaxis_title="Z"),
        coloraxis_colorbar=dict(yanchor="top", y=1, x=1.15, xanchor="left"),
    )

    # Show the plot
    fig.show()


def plot_mr(data, specs={}):
    fig = go.Figure()

    for label, dat in data.items():
        fig.add_trace(go.Scatter(x=dat[0] / m_earth, y=dat[1] / r_earth, name=label))
    fig.update_layout(
        xaxis_title=r"$\rm{Mass} \ (M_\oplus)$",
        yaxis_title=r"$\rm {Radius} \ (R_\oplus)$",
    )
    fig.show()


def subplot_labeler(fig, labels):
    pass
