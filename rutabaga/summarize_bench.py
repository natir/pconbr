#!/usr/bin/env python3

import os
import csv


def __read_info(filename):
    with open(filename) as finput:
        reader = csv.DictReader(finput, delimiter='\t')
        for record in reader:
            return {"time": record["s"], "memory": record["max_rss"]}

def __populate_entry(origin):
    data = dict()

    with os.scandir(origin.path) as it:
        for entry in it:
            if entry.is_file():
                data[entry.name.split(".")[0]] = __read_info(entry.path)

    return data

def read_bench_info():
    data = dict()
    
    with os.scandir("./benchmark") as it:
        for entry in it:
            if entry.is_dir():
                data[entry.name] = __populate_entry(entry)

    return data

def get_memory():
    data = read_bench_info()

    table = """
| counter | s_pneumoniae | c_vartiovaarae |
|:--------|-------------:|---------------:|
| Jellyfish | {} | {} |
| Kmc | {} | {} |
| Sssik | {} | {} |
""".format(data["jellyfish"]["s_pneumoniae"]["memory"], data["jellyfish"]["c_vartiovaarae"]["memory"], data["kmc"]["s_pneumoniae"]["memory"], data["kmc"]["c_vartiovaarae"]["memory"], data["ssik"]["s_pneumoniae"]["memory"], data["ssik"]["c_vartiovaarae"]["memory"])
    
    return table

def get_time():
    data = read_bench_info()

    table = """
| counter | s_pneumoniae | c_vartiovaarae |
|:--------|-------------:|---------------:|
| Jellyfish | {} | {} |
| Kmc | {} | {} |
| Sssik | {} | {} |
""".format(data["jellyfish"]["s_pneumoniae"]["time"], data["jellyfish"]["c_vartiovaarae"]["time"], data["kmc"]["s_pneumoniae"]["time"], data["kmc"]["c_vartiovaarae"]["time"], data["ssik"]["s_pneumoniae"]["time"], data["ssik"]["c_vartiovaarae"]["time"])
    
    return table
