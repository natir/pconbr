#!/usr/bin/env python3

import os
import csv

from collections import defaultdict

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

    datasets = list(data["pcon"].keys())
    header = "| dataset | k | Jellyfish | Kmc | Pconbr |\n|:-|:-|-:|-:|-:|\n"

    table = ""
    for dataset in datasets:
        for k in sorted(data["pcon"][dataset].keys()):
            table += "| {} | {} | {} | {} | {} |\n".format(dataset, k, data["jellyfish"][dataset][k][value], data["kmc"][dataset][k][value], data["pcon"][dataset][k][value])
    
    return header + table

