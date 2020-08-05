#!/usr/bin/env python3

import struct

import csv
import pandas

from collections import Counter, defaultdict 

def count(filename, columns_name):
    count = dict()
    with open(filename) as fh:
        reader = csv.reader(fh)
        for row in reader:
            count[row[0]] = int(row[1])

    df = pandas.DataFrame.from_dict(count, orient='index', columns=[columns_name])
    df.sort_index(inplace=True)

    return df

def generate_csv_count(inputs, output):
    # Get data
    df = pandas.DataFrame()
    for (k, v) in inputs.items():
        for value in v:
            df = pandas.merge(df, count(value, k), how="outer", left_index=True, right_index=True)

    # Clean up
    df = df.fillna(0)
    for c in df.columns:
        df = df.astype({c: 'int64'})

    # Write result
    df.to_csv(output)

def generator_true_set(filename):
    with open(filename) as fh:
        reader = csv.reader(fh)
        for row in reader:
            if int(row[1]) > 0:
                yield row[0]
    
def count_true_false(filename, true_set, columns_name):

    count = {True: Counter(), False: Counter()}

    with open(filename) as fh:
        reader = csv.reader(fh)
        for row in reader:
            count[row[0] in true_set][int(row[1])] += 1

    df = pandas.DataFrame.from_dict(dict(count[True]), orient='index', columns=[columns_name + "_true"])
    false_df = pandas.DataFrame.from_dict(dict(count[False]), orient='index', columns=[columns_name + "_false"])

    df = pandas.merge(df, false_df, how="outer", left_index=True, right_index=True)

    df.sort_index(inplace=True)
    return df

def generate_csv_true_false(inputs, output):

    k2file = defaultdict(dict)

    for v in inputs:
        k = v.split(".")[1][1:]
        if v.endswith("a0.csv"):
            k2file[k][True] = v
        else:
            k2file[k][False] = v

    df = pandas.DataFrame()
    for (k, v) in k2file.items():
        df = pandas.merge(df, count_true_false(v[False], set(generator_true_set(v[True])), k), how="outer", left_index=True, right_index=True)

    # Clean up
    df = df.fillna(0)
    for c in df.columns:
        df = df.astype({c: 'int64'})

    # Write result
    df.to_csv(output)
