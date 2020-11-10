import os
import csv
import pandas
import altair

from . import utils

def dataframe():
    data = list()
    
    for counter in ["pcon", "pcon_cd", "kmc", "kmc_cd", "kmc_disk", "kmc_disk_cd", "jellyfish", "jellyfish_cd"]:
        for dataset in utils.get_data_set("data"):
            for kmer_size in range(13, 21, 2):
                (time, memory, mean_load) = get_data(counter, dataset, "reads", kmer_size)
                size = utils.get_file_size(dataset)
                
                if time is not None:
                    data.append((counter, dataset, str(kmer_size), time, memory, mean_load, size))

    return pandas.DataFrame(data, columns=['counter', 'dataset', 'kmer_size', 'time', 'memory', 'mean_load','size'])

    
def get_data(counter, dataset, prefix, kmer_size):
    path = f"count/bench/{counter}/{dataset}_{prefix}.k{kmer_size}.tsv"

    return utils.get_bench_data(path)
    

def get_file_size(dataset):
    path = f"data/{dataset}/reads.fasta"

    return os.path.getsize(path)

