import os
import csv
import math

import altair

def get_bench_data(path):
    if os.path.isfile(path):
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            row = next(reader)
            return (float(row["s"]), float(row["max_rss"]))
    else:
        return (None, None)

def group_scatter(df, x, y, color=None, shape=None):

    xrange = (math.floor(df[x].min() * 100) / 100, math.ceil(df[x].max() * 100) / 100)
    yrange = (math.floor(df[y].min() * 100) / 100, math.ceil(df[y].max() * 100) / 100)

    fig =  altair.Chart(df).mark_point().encode(
            x=altair.X(x, scale=altair.Scale(domain=xrange)),
            y=altair.Y(y, scale=altair.Scale(domain=yrange))
        )

    if color:
        fig = fig.encode(color=color)
        
    if shape is not None:
        fig = fig.encode(shape=shape)
    
    return fig
