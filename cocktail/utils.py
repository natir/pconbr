import os
import csv

def get_bench_data(path):
    if os.path.isfile(path):
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            row = next(reader)
            return (int(row["s"]), int(row["max_rss"]))
    else:
        return (None, None)
