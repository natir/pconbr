import os
import csv
import math
import pandas
import altair

from .utils import get_bench_data


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


def dataframe_bench():
    data = list()
    
    for dataset in ["bacteria", "yeast", "metagenome"]:
        for kmer_size in range(13, 21, 2):
            (time, memory) = get_data_bench("br", dataset, f".k{kmer_size}")
            if time is not None:
                data.append((dataset, f"br_k{kmer_size}", time, memory))

        for corrector in ["canu", "consent", "necat"]:
            (time, memory) = get_data_bench(corrector, dataset, "")
            if time is not None:
                data.append((dataset, corrector, time, memory))

    return pandas.DataFrame(data, columns=['dataset', 'corrector', 'time', 'memory'])


def get_data_bench(corrector, dataset, params):
    path = f"correct/bench/{corrector}/{dataset}_reads{params}.tsv"
        
    return get_bench_data(path)
