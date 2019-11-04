#!/usr/bin/env python3

import struct

import pandas

from collections import Counter, defaultdict

BUFF_SIZE=2048

def buff_read(reader):
    while True:
        byte = reader.read(BUFF_SIZE)
        if len(byte) == 0:
            break
        
        for b in byte:
            yield b

def read_pcon_file(filename):
    with open(filename, 'rb') as reader:
        k, nb_bit = struct.unpack('BB', reader.read(2))
        for val in buff_read(reader):
            if nb_bit == 4:
                yield val & 0b1111
                yield (val & 0b11110000) >> 4
            elif nb_bit == 8:
                yield val
            else:
                raise StopIteration
    

def count(filename, columns_name):
    count = dict(Counter(read_pcon_file(filename)))

    df = pandas.DataFrame.from_dict(count, orient='index', columns=[columns_name])
    df.sort_index(inplace=True)

    return df

def generate_csv_count(inputs, output):
    # Get data
    df = pandas.DataFrame()
    for (k, v) in inputs.items():
        df = pandas.merge(df, count(v, k), how="outer", left_index=True, right_index=True)

    # Clean up
    df = df.fillna(0)
    for c in df.columns:
        df = df.astype({c: 'int64'})

    # Write result
    df.to_csv(output)

def generate_true_set(filename):
    return {i for i, c in enumerate(read_pcon_file(filename)) if int(c) > 0}
    
def count_true_false(filename, true_set, columns_name):

    count = {True: Counter(), False: Counter()}

    for (i, v) in enumerate(read_pcon_file(filename)):
        count[i in true_set][v] += 1

    df = pandas.DataFrame.from_dict(dict(count[True]), orient='index', columns=[columns_name + "_true"])
    df.sort_index(inplace=True)

    false_df = pandas.DataFrame.from_dict(dict(count[False]), orient='index', columns=[columns_name + "_false"])
    false_df.sort_index(inplace=True)

    df = pandas.merge(df, false_df, how="outer", left_index=True, right_index=True)

    return df

def generate_csv_true_false(inputs, output):

    k2file = defaultdict(dict)

    for (k, v) in inputs.items():
        if k.startswith("true_"):
            k2file[k.split("_")[1]][True] = v
        else:
            k2file[k][False] = v

    df = pandas.DataFrame()
    for (k, v) in k2file.items():
        df = pandas.merge(df, count_true_false(v[False], generate_true_set(v[True]), k), how="outer", left_index=True, right_index=True)

    # Clean up
    df = df.fillna(0)
    for c in df.columns:
        df = df.astype({c: 'int64'})

    # Write result
    df.to_csv(output)
