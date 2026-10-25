import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import coolbox.api as cb
from capcruncher.api.plotting import CCTrack, CCFigure
import pyranges as pr


viewpoints = pr.read_bed("pilot_capture_analysis/plotting_viewpoints.bed")
viewpoints

tracks = [
    CCTrack(None, type="scale"),
    CCTrack(
        "pilot_capture_analysis/capcruncher_output/results/DAAP6matureadipocytes_M1/DAAP6matureadipocytes_M1.hdf5::/chr1_9990409_9991833",
        type="heatmap",
        binsize=10000,
        title="Alpha Tile",
        viewpoint="Alpha",
        normalization="ice",
        transform="yes",
        style="triangular",
    ),
    CCTrack(None, type="spacer"),
    CCTrack(None, type="spacer"),
    CCTrack(None, type="xaxis"),
]

fig = CCFigure(tracks, auto_spacing=False)

fig.show()

#fig.save("output.png")

# Plot a specific region if desired
fig.plot("chr8_26839367_26841115")