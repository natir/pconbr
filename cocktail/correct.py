import os
import csv
import math
import pandas
import altair

from . import utils

def dataframe_elector():
    data = list()

    for dataset in utils.get_data_set("data"):
        for kmer_size in range(13, 21, 2):
            (precision, recall) = get_data_elector("br", dataset, kmer_size)
            
            if precision is not None and recall is not None:
                data.append((f"br_k{kmer_size}", dataset, precision, recall))
    
    for corrector in ["canu", "consent", "necat"]:
        for dataset in ["bacteria", "yeast", "metagenome"]:
                (precision, recall) = get_data_elector(corrector, dataset, None)
                
                if precision is not None and recall is not None:
                    data.append((corrector, dataset, precision, recall))
                        
    return  pandas.DataFrame(data, columns=['corrector', 'dataset', 'precision', 'recall'])

    
def get_data_elector(corrector, dataset, kmer_size=None):
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

    d2size = {d: size / 1_000_000 for (d, size) in map(lambda x: (x, utils.get_file_size(x)), utils.get_data_set("data")) if size is not None}
        
    for dataset in d2size.keys():
        for kmer_size in range(13, 21, 2):
            (time, memory) = get_data_bench("br", dataset, f".k{kmer_size}")

            if time is not None:
                data.append((dataset, f"br_k{kmer_size}", time, memory, d2size[dataset]))
                
        for corrector in ["canu", "consent", "necat"]:
            (time, memory) = get_data_bench(corrector, dataset, "")

            if time is not None:
                data.append((dataset, corrector, time, memory, d2size[dataset]))
                
    return pandas.DataFrame(data, columns=['dataset', 'corrector', 'time', 'memory', 'size'])


def get_data_bench(corrector, dataset, params):
    path = f"correct/bench/{corrector}/{dataset}_reads{params}.tsv"
        
    return utils.get_bench_data(path)


def dataframe_stats():
    data = list()
    
    d2error = {d: e for (d, e) in map(lambda x: (x, utils.get_error_rate_raw(x)), utils.get_data_set("data")) if e is not None}
       
    for dataset in d2error.keys():
        for kmer_size in range(13, 21, 2):
            error_rate = get_error_rate("br", dataset, f".k{kmer_size}")

            if error_rate is not None:
                data.append((dataset, f"br_k{kmer_size}", error_rate, d2error[dataset]))
                
        for corrector in ["canu", "consent", "necat"]:
            error_rate = get_error_rate(corrector, dataset, "")

            if error_rate is not None:
                data.append((dataset, corrector, error_rate, d2error[dataset]))
                
    return pandas.DataFrame(data, columns=['dataset', 'corrector', 'corrected', 'raw'])


def get_error_rate(corrector, dataset, params):
    path = f"correct/{dataset}/{corrector}/reads{params}.stats"

    return utils.get_error_rate(path)

