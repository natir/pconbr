#!/usr/bin/env python3

import struct

import pandas

from collections import Counter

def read_pcon_file(filename):
    with open(filename, 'rb') as reader:
        k, nb_bit = struct.unpack('BB', reader.read(2))
        val = read_a_value(reader)
        while val != -1:
            if nb_bit == 4:
                yield str(val & 0b1111)
                yield str((val & 0b11110000) >> 4)
            elif nb_bit == 8:
                yield str(val)
            else:
                raise StopIteration
                
            val = read_a_value(reader)

            
def read_a_value(reader):
    byte = reader.read(1)
    if byte == b'':
        return -1

    return struct.unpack('B', byte)[0]

def get(filename, columns_name):
    count = dict(Counter(read_pcon_file(filename)))

    df = pandas.DataFrame.from_dict(count, orient='index', columns=[columns_name])
    df.index = df.index.astype(int)
    df.sort_index(inplace=True)

    return df
