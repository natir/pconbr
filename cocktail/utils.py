import os
import csv
import math
import collections

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
        
    if shape:
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
        
    if shape:
        fig = fig.encode(shape=shape)

    if column:
        fig = fig.encode(column=column)
    
    return fig

def get_data_set(step):
    return sorted([entry.name for entry in os.scandir(step) if entry.is_dir()])


def get_error_rate_raw(dataset):
    path = f"data/{dataset}/reads.stats"

    return get_error_rate(path)


def get_error_rate(path):
    if not os.path.isfile(path):
        return None
    
    with open(path) as fh:
        for line in fh:
            if line.startswith("SN\terror rate:"):
                return float(line.split("\t")[2])

    return None


def get_quast_info(path):
    if not os.path.isfile(path):
        return None

    QuastResult = collections.namedtuple('QuastResult', ["nb_contigs",
                                                         "nb_misassemblies",
                                                         "unaligned_length",
                                                         "genome_fraction",
                                                         "nb_mismatches",
                                                         "nb_indels",
                                                         "largest_alignment",
                                                         "NGA50"])

    ret = QuastResult(None, None, None, None, None, None, None, None)

    with open(path) as fh:
        reader = csv.reader(fh, delimiter='\t')

        for row in reader:
            if row[0] == "# contigs":
                ret = ret._replace(nb_contigs=int(row[1]))
            elif row[0] == "# misassemblies":
                ret = ret._replace(nb_misassemblies=int(row[1]))
            elif row[0] == "Unaligned length":
                ret = ret._replace(unaligned_length=int(row[1]))
            elif row[0] == "Genome fraction (%)":
                ret = ret._replace(genome_fraction=float(row[1]))
            elif row[0] == "# mismatches per 100 kbp":
                ret = ret._replace(nb_mismatches=float(row[1]))
            elif row[0] == "# indels per 100 kbp":
                ret = ret._replace(nb_indels=float(row[1]))
            elif row[0] == "Largest alignment":
                ret = ret._replace(largest_alignment=int(row[1]))
            elif row[0] == "NGA50":
                ret = ret._replace(NGA50=int(row[1]))

    if any((e is None for e in ret)):
        return None
    else:
        return ret


def fig_layout(figs, cols):
    fig = None
    
    for i in range(0, len(figs), cols):
        end = i + cols
        if end > len(figs):
            end = len(figs)

        line = None
        for f in figs[i:end]:
            if line is None:
                line = f
            else:
                line |= f

        if fig is None:
            fig = line
        else:
            fig &= line
        
    return fig
