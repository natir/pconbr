import os
import csv
import pandas
import altair

def figure(dataset, tool):
    df = dataframe()
    df = df[(df["dataset"] == dataset) & (df["tool"] == tool)]

    readlength = altair.Chart(df).mark_bar().encode(
        x=altair.X("length:Q", bin=altair.Bin(maxbins=100), axis=None),
        y=altair.Y('count()', axis=altair.Axis(title=None))
    ).properties(
        width=1000,
        height=200
    )

    quality = altair.Chart(df).mark_bar().encode(
        y=altair.X("identity:Q", bin=altair.Bin(maxbins=100), axis=None),
        x=altair.Y('count()', axis=altair.Axis(title=None))
    ).properties(
        width=200,
        height=500
    )

    scatter = altair.Chart(df).mark_circle(size=10).encode(
        x=altair.X("length", axis=altair.Axis(title="Read length in base")),
        y=altair.Y("identity", axis=altair.Axis(title="Read indentity in percent"))
    ).properties(
        width=1000,
        height=500
    )

    return altair.vconcat(readlength, scatter | quality) 

def dataframe():
    data = list()

    for dataset in ["bacteria", "yeast", "metagenome"]:
        for kmer_size in range(13, 21, 2):
            for ratio in range(70, 100, 5):
                path = f"filter/{dataset}/identity/kmrf/reads.k{kmer_size}.r{ratio}.tsv"
                if os.path.isfile(path):
                    with open(path) as fh:
                        reader = csv.DictReader(fh, delimiter='\t')
                        for row in reader:
                            data.append((dataset, f"kmrf_k{kmer_size}_r{ratio}", row["name"], int(row["length"]), float(row["identity"])))

    for dataset in ["bacteria", "yeast", "metagenome"]:
        for quality in range(90, 100, 1):
                path = f"filter/{dataset}/identity/filtlong/reads.q{quality}.tsv"
                if os.path.isfile(path):
                    with open(path) as fh:
                        reader = csv.DictReader(fh, delimiter='\t')
                        for row in reader:
                            data.append((dataset, f"filtlong_q{quality}", row["name"], int(row["length"]), float(row["identity"])))
                            
    return pandas.DataFrame(data, columns=['dataset', 'tool', 'name', 'length', 'identity'])
