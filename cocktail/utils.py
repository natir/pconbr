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


def get_file_size(dataset):
    path = f"data/{dataset}/reads.fasta"

    if not os.path.isfile(path):
        return None
    
    return os.path.getsize(path)

    
def group_scatter(df, x, y, color=None, shape=None, title="", xtitle="", ytitle=""):

    xrange = (math.floor(df[x].min() * 100) / 100, math.ceil(df[x].max() * 100) / 100)
    yrange = (math.floor(df[y].min() * 100) / 100, math.ceil(df[y].max() * 100) / 100)

    fig = altair.Chart(df, title=title).mark_point().encode(
        x=altair.X(x, scale=altair.Scale(domain=xrange), title=xtitle),
        y=altair.Y(y, scale=altair.Scale(domain=yrange), title=ytitle)
    )

    if color:
        fig = fig.encode(color=color)
        
    if shape is not None:
        fig = fig.encode(shape=shape)

    return fig


def group_line(df, x, y, color=None, shape=None, column=None, title="", xtitle="", ytitle="", point=False):

    xrange = (math.floor(df[x].min() * 100) / 100, math.ceil(df[x].max() * 100) / 100)
    yrange = (math.floor(df[y].min() * 100) / 100, math.ceil(df[y].max() * 100) / 100)

    fig = altair.Chart(df, title=title).mark_line(point=point).encode(
        x=altair.X(x, scale=altair.Scale(domain=xrange), title=xtitle),
        y=altair.Y(y, scale=altair.Scale(domain=yrange), title=ytitle)
    )

    if color:
        fig = fig.encode(color=color)
        
    if shape is not None:
        fig = fig.encode(shape=shape)

    if column is not None:
        fig = fig.encode(column=column)
    
    return fig

def get_data_set(step):
    for entry in os.scandir(step):
        if entry.is_dir():
            yield entry.name
