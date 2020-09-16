import os
import csv
import pandas
import altair


def figure(dataset):
    df = dataframe()

    return altair.Chart(df[df["dataset"] == dataset]).mark_point().encode(
        x="time",
        y="memory",
        color="counter",
        shape="kmer_size"
    )

def dataframe():
    data = list()
    
    for counter in ["pcon", "kmc", "jellyfish"]:
        for dataset in ["bacteria", "yeast", "metagenome"]:
            for kmer_size in range(13, 21, 2):
                (time, memory) = get_data(counter, dataset, "reads", kmer_size)

                if time is not None:
                    data.append((counter, dataset, str(kmer_size), time, memory))

    return  pandas.DataFrame(data, columns=['counter', 'dataset', 'kmer_size', "time", "memory"])

    
def get_data(counter, dataset, prefix, kmer_size):
    path = f"count/bench/{counter}/{dataset}_{prefix}.k{kmer_size}.tsv"

    if os.path.isfile(path):
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            row = next(reader)
            return (int(row["s"]), int(row["max_rss"]))
    else:
        return (None, None)
