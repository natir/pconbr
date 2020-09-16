import os
import csv
import pandas
import altair

def figure(dataset, kmer_size, ratio):
    df = dataframe()
    df = df[(df["kmer_size"] == str(kmer_size)) & (df["dataset"] == dataset) & (df["ratio"] == ratio)]

    readlength = altair.Chart(df).mark_bar().encode(
        x=altair.X("length:Q", bin=altair.Bin(maxbins=100), axis=None),
        y=altair.Y('count()', axis=None)
    )

    quality = altair.Chart(df).mark_bar().encode(
        y=altair.X("identity:Q", bin=altair.Bin(maxbins=100), axis=None),
        x=altair.Y('count()', axis=None)
    )

    scatter = altair.Chart(df).mark_circle().encode(
        x="length",
        y="identity"
    )

    return altair.vconcat(readlength, scatter | quality) 

def dataframe():
    data = list()

    for dataset in ["bacteria", "yeast", "metagenome"]:
        for kmer_size in range(13, 21, 2):
            for ratio in range(6, 10, 1):
                path = f"filter/{dataset}/identity/reads.k{kmer_size}.r{ratio/10}.tsv"
                if os.path.isfile(path):
                    with open(path) as fh:
                        reader = csv.DictReader(fh, delimiter='\t')
                        for row in reader:
                            data.append((dataset, str(kmer_size), ratio, row["name"], int(row["length"]), float(row["identity"])))

    return pandas.DataFrame(data, columns=['dataset', 'kmer_size', 'ratio', 'name', 'length', 'identity'])
