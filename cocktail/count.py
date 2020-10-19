import os
import csv
import pandas
import altair

from .utils import get_bench_data

def dataframe():
    data = list()
    
    for counter in ["pcon", "kmc", "jellyfish"]:
        for dataset in ["bacteria", "yeast", "metagenome", "bacteria5", "bacteria7"]:
            for kmer_size in range(13, 21, 2):
                (time, memory) = get_data(counter, dataset, "reads", kmer_size)
                size = get_file_size(dataset)
                
                if time is not None:
                    data.append((counter, dataset, str(kmer_size), time, memory, size))

    return pandas.DataFrame(data, columns=['counter', 'dataset', 'kmer_size', 'time', 'memory', 'size'])

    
def get_data(counter, dataset, prefix, kmer_size):
    path = f"count/bench/{counter}/{dataset}_{prefix}.k{kmer_size}.tsv"

    return get_bench_data(path)
    

def get_file_size(dataset):
    path = f"data/{dataset}/reads.fasta"

    return os.path.getsize(path)

