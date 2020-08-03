#!/usr/bin/env python3

import os
import re
import csv

from collections import defaultdict

import pandas

def get_parameter(filename):
    simulated_reads_96.k11.a18.one-greedy-gap_size.stats
    match = re.match(".+\.k(?P<k>\d+)(\.a(?P<a>\d+))?\.(?P<method>.+)\.stats", filename)

    ret = {k: int(v) for k, v in match.groupdict().items() if v is not None}
    return ret

def get_error_rate(filename):
    with open(filename) as finput:
        reader = csv.reader(finput, delimiter='\t')
        for record in reader:
            if record[0] == "error rate:":
                return record[1]

def genomic_kmer(filename):
    raw_error_rate = float(get_error_rate("reads/{}.stats".format(filename)))
    
    data = defaultdict(lambda: defaultdict(float))
    with os.scandir("./genetic_kmer/") as it:
        for entry in it:
            if entry.is_file() and entry.name.startswith(filename) and entry.name.endswith(".stats"):
                param = get_parameter(entry.name)
                k = param['k']
                m = param['method']
                
                error_rate = get_error_rate(entry.path)
                if error_rate is None:
                    data[k][m] = -1
                else:
                    data[k][m] = float(error_rate) - raw_error_rate

    all_k = sorted(list({k for k in data.keys()}))
    all_m = sorted(list({s for k in data.keys() for m in data[k].keys()}))

    df = pandas.DataFrame(index=[k for k in all_k])
    df.index.name = "k"
    
    for m in all_m:
        tmp_list = list()
        for k in df.index:
            if k not in data:
                tmp_list.append(None)
            elif s not in data[k]:
                tmp_list.append(None)
            else:
                tmp_list.append(data[k][s] * 100)
        
        df["{}".format(m)] = tmp_list

    return df


def read_kmer(filename):
    raw_error_rate = float(get_error_rate("reads/{}.stats".format(filename)))

    data = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    with os.scandir("./read_kmer/") as it:
        for entry in it:
            if entry.is_file() and entry.name.startswith(filename) and entry.name.endswith(".stats"):
                param = get_parameter(entry.name)
                k = param["k"]
                a = param["a"]
                m = param["m"]
                
                error_rate = get_error_rate(entry.path)
                if error_rate is None:
                    data[k][a][m] = -1
                else:
                    data[k][a][m] = float(error_rate) - raw_error_rate
                    
    all_k = sorted(list({k for k in data.keys()}))
    all_a = sorted(list({a for k in data.keys() for a in data[k].keys()}))                
    all_m = sorted(list({s for k in data.keys() for a in data[k].keys() for m in data[k][a].keys()}))

    index = pandas.MultiIndex.from_tuples([(a, k) for a in all_a for k in all_k], names=("a", "k"))
    df = pandas.DataFrame(index=index)

    for m in all_m:
        tmp_list = list()
        for a, k in df.index:
            if k not in data:
                tmp_list.append(None)
            elif a not in data[k]:
                tmp_list.append(None)
            elif s not in data[k][a]:
                tmp_list.append(None)
            else:
                tmp_list.append(data[k][a][m] * 100)

        df["{}".format(m)] = tmp_list
                
    return df
