import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from pics.runparams import param_colors
from pics.physicalparams import m_earth, r_earth, G
import pandas as pd
import numpy as np


def plot_structure(profiles, file_name="structure"):
    """Plots the structure profils of a planetary object"""

    df = pd.DataFrame(
        np.array(profiles).T,
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

    scalings = [1e-3, 1., 1e-9, 1., 1./m_earth, 1., 1., 1.]

    fig = make_subplots(
        rows=2, cols=2, shared_xaxes=True, horizontal_spacing=0.25, vertical_spacing=0.1
    )

    fig.update_layout(height=500, width=600, showlegend=False)
    params = [["temperature", "density"], ["pressure", "mass"]]
    units = [["K", "km \ m^{-3}"], ["GPa", "\mathit{M_\oplus}"]]

    for i in range(2):
        for j in range(2):
            fig.add_trace(
                go.Scatter(
                    x=df["radius"]*1e-3,
                    y=df[params[i][j]]*scalings[i+2*j+1],
                    line=dict(color="rgb{}".format(param_colors[2*i + j])),
                ),
                row=i + 1,
                col=j + 1,
            )

            if i + 2 * j > 0:
                a = str(2 * i + j + 1)
            else:
                a = ""

            fig["layout"]["yaxis{}".format(a)]["title"] = r"$\rm {} \ ({})$".format(
                params[i][j].capitalize(), units[i][j]
            )

    fig["layout"]["xaxis3"]["title"] = r"$\rm {} \ ({})$".format("Radius", "km")
    fig["layout"]["xaxis4"]["title"] = r"$\rm {} \ ({})$".format("Radius", "km")

    fig.show()


def subplot_labeler(fig, labels):
    pass
