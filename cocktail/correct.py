import os
import csv
import math
import pandas
import altair

def figure_pr(dataset):
    df = dataframe_pr()

    xrange = (math.floor(df["recall"].min() * 100) / 100, math.ceil(df["recall"].max() * 100) / 100)
    yrange = (math.floor(df["precision"].min() * 100) / 100, math.ceil(df["precision"].max() * 100) / 100)
    
    return altair.Chart(df[df["dataset"] == dataset]).mark_point().encode(
        x=altair.X("recall", scale=altair.Scale(domain=xrange)),
        y=altair.Y("precision", scale=altair.Scale(domain=yrange)),
        color="corrector",
    )


def dataframe_pr():
    data = list()

    for dataset in ["bacteria", "yeast", "metagenome"]:
        for kmer_size in range(13, 21, 2):
            (precision, recall) = get_data_pr("br", dataset, kmer_size)
            
            if precision is not None and recall is not None:
                data.append((f"br_k{kmer_size}", dataset, precision, recall))
    
    for corrector in ["canu", "consent", "necat"]:
        for dataset in ["bacteria", "yeast", "metagenome"]:
                (precision, recall) = get_data_pr(corrector, dataset, None)
                
                if precision is not None and recall is not None:
                    data.append((corrector, dataset, precision, recall))
                        
    return  pandas.DataFrame(data, columns=['corrector', 'dataset', 'precision', 'recall'])

    
def get_data_pr(corrector, dataset, kmer_size=None):
    if kmer_size is not None:
        path = f"correct/{dataset}/elector/{corrector}/reads.k{kmer_size}/log"
    else:
        path = f"correct/{dataset}/elector/{corrector}/reads/log"
        
    if os.path.isfile(path):
        precision = None
        recall = None
        with open(path) as fh:
            for line in fh:
                if line.startswith("Precision"):
                    precision = float(line.split(":")[-1])

                if line.startswith("Recall"):
                    recall = float(line.split(":")[-1])

        return (precision, recall)
    
    return (None, None)
