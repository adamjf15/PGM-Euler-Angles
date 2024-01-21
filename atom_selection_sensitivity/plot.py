import pandas as pd
import numpy as np
import plotly.express as px

datafiles = [
    "atomsel1/eulerAngles",
    "atomsel2/eulerAngles",
    "atomsel3/eulerAngles",
]

np.random.seed(13)

dfs = []
for file in datafiles:
    df = pd.read_table(file, delim_whitespace=True)
    df["atomsel"] = int(file.split("/")[0][-1])
    dfs.append(df)
df = pd.concat(dfs, axis=0)
df["atomsel"] = pd.Categorical(df["atomsel"])

random_structures = df["structure"].sample(12).to_list()
df = df.loc[df["structure"].isin(random_structures), :]

fig = px.scatter_3d(
    df,
    x="yaw",
    y="pitch",
    z="roll",
    symbol="atomsel",
    color="structure"
)

fig.update_traces(
    marker_size=8,
    marker_opacity=1,
    marker_line_width=0.5,
)

fig.update_layout(
    template="none",
    scene_xaxis_range=[-12,5],
    scene_yaxis_range=[-10,30],
    scene_zaxis_range=[-25,5],
    scene_aspectratio={"x":1, "y":1, "z":1},
    height=800, 
    width=800,
    showlegend=False,
    margin=dict(l=25, r=25, b=25, t=25),
)

fig.write_image("atom_selection_sensitivity.svg")

fig.write_html("atom_selection_sensitivity.html")