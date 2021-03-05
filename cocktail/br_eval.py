import os
import csv
import math
import pandas
import altair
import itertools

from . import utils


def get_data_bench(corrector, dataset, params):
    path = f"correct/bench/{corrector}/{dataset}_reads{params}.tsv"

    return utils.get_bench_data(path)


def dataframe_stats():
    data = list()

    d2error = {d: e for (d, e) in map(lambda x: (x, utils.get_error_rate_raw(x)), utils.get_data_set("data")) if e is not None}

    methods = ["one", "two", "greedy", "gap_size", "graph"]
    for method in itertools.permutations(
            ["one", "two", "greedy", "gap_size", "graph"]
    ):
        method = "-".join(method)
        methods.append(method)

    abundances = list()
    abundances.append("first-minimum")
    for method in ["rarefaction", "percent-most", "percent-least"]:
        for p in range(10, 40, 10):
            abundances.append(f"{method}_{p/100}")

    for dataset in d2error.keys():
        for kmer_size in range(15, 17, 2):
            for method in ["one", "two", "greedy", "gap_size", "graph"]:
                for abundance in abundances:
                    error_rate = get_error_rate("br", dataset, f".k{kmer_size}", method, abundance)

                    if error_rate is not None:
                        data.append((dataset, f"br_k{kmer_size}", method, abundance, error_rate, d2error[dataset]))

    return pandas.DataFrame(data, columns=['dataset', 'kmer_size', 'method', 'abundance', 'corrected', 'raw'])


def get_error_rate(corrector, dataset, kmer_size, method, abundance):
    path = f"br_eval/stats/{dataset}/reads.k{kmer_size}.m{method}.a{abundance}.stats"

    return utils.get_error_rate(path)
