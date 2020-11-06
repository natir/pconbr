import os
import csv
import pandas

from . import utils

def assembly_data():
    data = list()

    for dataset in utils.get_data_set("data"):
        res = get_quast_info(dataset, "")
        if res is not None:
            data.append(dataset, "raw", "raw", *res)

        for k in range(13, 21, 2):
            for a in range(10, 31):
                res = get_quast_info(dataset, f".k{k}.a{a}")
                if res is not None:
                    data.append((dataset, str(k), str(a), *res))

    return pandas.DataFrame(data, columns=["dataset", "kmer_size", "abundance",
                                           "nb_contigs", "nb_misassemblies",
                                           "unaligned_length", "genome_fraction",
                                           "nb_mismatches", "nb_indels",
                                           "largest_alignment", "NGA50"])

        
def get_quast_info(dataset, params):
    path = f"polish/{dataset}/quast/polish{params}/"

    return utils.get_quast_info(path)
