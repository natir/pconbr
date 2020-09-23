import os
import csv
import pandas
import altair

def figure(df):
    domain_len = (df["length"].min(), df["length"].max())
    domain_ide = (df["identity"].min(), 100)
    
    readlength = altair.Chart(df).mark_bar().encode(
        x=altair.X("length:Q", bin=altair.Bin(maxbins=100), scale=altair.Scale(domain=domain_len), axis=None),
        y=altair.Y('count()', axis=altair.Axis(title=None))
    ).properties(
        width=1000,
        height=200
    )

    quality = altair.Chart(df).mark_bar().encode(
        y=altair.X("identity:Q", bin=altair.Bin(maxbins=100), scale=altair.Scale(domain=domain_ide), axis=None),
        x=altair.Y('count()', axis=altair.Axis(title=None)),
    ).properties(
        width=200,
        height=500
    )

    scatter = altair.Chart(df).mark_circle(size=10).encode(
        x=altair.X("length", axis=altair.Axis(title="Read length in base"), scale=altair.Scale(domain=domain_len)),
        y=altair.Y("identity", axis=altair.Axis(title="Read indentity in percent"), scale=altair.Scale(domain=domain_ide))
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
                    data += get_data(path, dataset, f"kmrf_k{kmer_size}_r{ratio}")
                    
    for dataset in ["bacteria", "yeast", "metagenome"]:
        for quality in range(90, 100, 1):
            path = f"filter/{dataset}/identity/filtlong/reads.q{quality}.tsv"
            if os.path.isfile(path):
                data += get_data(path, dataset, f"filtlong_q{quality}")
                            
    return pandas.DataFrame(data, columns=['dataset', 'tool', 'name', 'length', 'identity'])


def get_data(path, dataset, filter):
    data = list()
    
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            data.append((dataset, filter, row["name"], int(row["length"]), float(row["identity"])))

    return data
