#!/usr/bin/env python3

import os
import csv

from collections import defaultdict

import pandas

def __read_info(filename):
    with open(filename) as finput:
        reader = csv.DictReader(finput, delimiter='\t')
        for record in reader:
            return {"time": record["s"], "memory": record["max_rss"]}

def __populate_entry(origin):
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    with os.scandir(origin.path) as it:
        for entry in it:
            if entry.is_file():
                dataset, k, *_ = entry.name.split(".")
                data[dataset][k] = __read_info(entry.path)

    return data

def read_bench_info():
    data = dict()
    
    with os.scandir("./benchmark") as it:
        for entry in it:
            if entry.is_dir():
                data[entry.name] = __populate_entry(entry)

    return data

def get(value):
    data = read_bench_info()

    datasets = sorted({dataset for tool in data.keys() for dataset in data[tool]})
    all_k = sorted({k for dataset in datasets for k in data["pcon"][dataset]})

    index = pandas.MultiIndex.from_tuples([(d, k) for d in datasets for k in all_k], names=("dataset", "k"))
    df = pandas.DataFrame(index=index)

    df["jellyfish"] = [data["jellyfish"][dataset][k][value] for dataset in datasets for k in all_k]
    df["kmc"]       = [data["kmc"][dataset][k][value] for dataset in datasets for k in all_k]
    df["pcon"]      = [data["pcon"][dataset][k][value] for dataset in datasets for k in all_k]
    
    return df

