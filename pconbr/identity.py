#!/usr/bin/env python3

import os
import re
import csv

from collections import defaultdict

def get_parameter(filename):
    match = re.match(".+\.k(\d+)(\.a(\d+))?\.s(\d+).stats", filename)

    if match.group(2) is None:
        return (int(match.group(1)), int(match.group(4)), None)
    else:
        return (int(match.group(1)), int(match.group(3)), int(match.group(4)))

def get_error_rate(filename):
    with open(filename) as finput:
        reader = csv.reader(finput, delimiter='\t')
        for record in reader:
            if record[0] == "error rate:":
                return record[1]

def genomic_kmer():
    raw_error_rate = float(get_error_rate("reads/simulated_reads.stats"))
    
    data = defaultdict(lambda: defaultdict(float))
    with os.scandir("./genetic_kmer/") as it:
        for entry in it:
            if entry.is_file() and entry.name.endswith(".stats"):
                k, s, _ = get_parameter(entry.name)
                error_rate = get_error_rate(entry.path)
                if error_rate is None:
                    data[k][s] = -1
                else:
                    data[k][s] = float(error_rate) - raw_error_rate

    all_k = sorted(list({k for k in data.keys()}))
    all_s = sorted(list({s for k in data.keys() for s in data[k].keys()}))

    table = "| "
    for s in all_s:
        table += "| s{}".format(s)
    table += "|\n|:-"

    for s in all_s:
        table += "|-:"
    table += "|\n"

    for k in all_k:
        table += "| k{}".format(k)
        for s in all_s:
            table += " | {:.6f}".format(data[k][s]*100)
        table += "|\n"

    return table

def read_kmer():
    raw_error_rate = float(get_error_rate("reads/simulated_reads.stats"))
        
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    with os.scandir("./read_kmer/") as it:
        for entry in it:
            if entry.is_file() and entry.name.endswith(".stats"):
                k, a, s = get_parameter(entry.name)
                error_rate = get_error_rate(entry.path)
                if error_rate is None:
                    data[k][a][s] = -1
                else:
                    data[k][a][s] = float(error_rate) - raw_error_rate
                    
    all_k = sorted(list({k for k in data.keys()}))
    all_a = sorted(list({a for k in data.keys() for a in data[k].keys()}))                
    all_s = sorted(list({s for k in data.keys() for a in data[k].keys() for s in data[k][a].keys()}))

    table = "| | "
    for s in all_s:
        table += "| s{}".format(s)
    table += "|\n|:-|:-"

    for s in all_s:
        table += "|-:"
    table += "|\n"

    for k in all_k:
        for a in all_a:
            table += "| k{} | a{}".format(k, a)

            for s in all_s:
                table += " | {:.6f}".format(data[k][a][s]*100)
            table += "|\n"
            
    return table
