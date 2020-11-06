import os
import csv
import pandas
import altair

from . import utils

def identity_data(dataset, filter, params):
    data = list()

    path = f"filter/{dataset}/identity/{filter}/reads.{params}.tsv"      
    return get_id_data(path, dataset, filter, params)
 

def get_id_data(path, dataset, filter, params):
    data = list()
    
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            data.append((dataset, filter, params, row["name"], int(row["length"]), float(row["identity"])))

    return pandas.DataFrame(data, columns=['dataset', 'tool', 'params', 'name', 'length', 'identity'])


def figure_len_identity(df, len_filter=0.998, color=None):
    max = df.quantile(len_filter).length
    df = df.loc[df.length<max]
    
    domain_len = (df.length.min(), df.length.max())
    domain_ide = (df.identity.min(), 100)
    
    readlength = altair.Chart(df).mark_bar().encode(
        x=altair.X("length:Q", bin=altair.Bin(maxbins=200), scale=altair.Scale(domain=domain_len), axis=None),
        y=altair.Y('count()', axis=altair.Axis(title=None))
    ).properties(
        width=1000,
        height=200
    )

    quality = altair.Chart(df).mark_bar().encode(
        y=altair.X("identity:Q", bin=altair.Bin(maxbins=200), scale=altair.Scale(domain=domain_ide), axis=None),
        x=altair.Y('count()', axis=altair.Axis(title=None)),
    ).properties(
        width=200,
        height=500
    )

    scatter = altair.Chart(df).mark_circle(size=2).encode(
        x=altair.X("length", axis=altair.Axis(title="Read length in base"), scale=altair.Scale(domain=domain_len)),
        y=altair.Y("identity", axis=altair.Axis(title="Read indentity in percent"), scale=altair.Scale(domain=domain_ide)),
    ).properties(
        width=1000,
        height=500
    )

    if color:
        scatter = scatter.encode(color=color)
        readlength = readlength.encode(color=color)
        quality = quality.encode(color=color)

    return readlength & (scatter | quality) 

def figure_filtred(dataframe, filter, params, len_filter=0.998):
    df_raw = get_id_data(f"data/{dataframe}/reads.len_id.tsv", "yeast", "raw", "")

    df_filter = identity_data(dataframe, filter, params)
    filter_name = set(df_filter.name)

    df_raw["filtered"] = df_raw.name.apply(lambda x: x not in filter_name)

    return figure_len_identity(df_raw, len_filter=len_filter, color="filtered")


def dataframe_bench():
    data = list()

    for dataset in ["bacteria", "yeast", "metagenome", "bacteria5", "bacteria7"]:
        for qual in range(90, 100):
            (time, memory) = get_data_bench("filtlong", dataset, f"q{qual}")
            size = utils.get_file_size(dataset)
                
            if time is not None:
                data.append(("filtlong", dataset, f"q{qual}", time, memory, size))

        for kmer_size in range(13, 21, 2):
            for ratio in range(70, 100, 5):
                (time, memory) = get_data_bench("kmrf", dataset, f"k{kmer_size}.r{ratio}")
                size = utils.get_file_size(dataset)
                
                if time is not None:
                    data.append(("kmrf", dataset, f"k{kmer_size}.r{ratio}", time, memory, size))
                    
    return pandas.DataFrame(data, columns=['filter', 'dataset', 'params', 'time', 'memory', 'size'])


def get_data_bench(filter, dataset, params):
    path = f"filter/bench/{filter}/{dataset}_reads.{params}.tsv"

    return utils.get_bench_data(path)
    

def get_error_rate(filter, dataset, params):
    path = f"filter/{dataset}/{filter}/reads{params}.stats"

    return utils.get_error_rate(path)


def dataframe_stats():
    data = list()

    d2error = {d: e for (d, e) in map(lambda x: (x, utils.get_error_rate_raw(x)), utils.get_data_set("data")) if e is not None}
         
    for dataset in d2error.keys():
        for kmer_size in range(13, 21, 2):
            for ratio in range(70, 100, 5):
                error_rate = get_error_rate("kmrf", dataset, f".k{kmer_size}.r{ratio}")

                if error_rate is not None:
                    data.append((dataset, f"kmrf.k{kmer_size}.r{ratio}", error_rate, d2error[dataset]))
                
        for qual in range(90, 100):
            error_rate = get_error_rate("filtlong", dataset, f".q{qual}")

            if error_rate is not None:
                data.append((dataset, f"filtlong.q{qual}", error_rate, d2error[dataset]))
                
    return pandas.DataFrame(data, columns=['dataset', 'filter', 'corrected', 'raw'])


def assembly_data():
    data = list()

    for dataset in utils.get_data_set("data"):
        res = get_quast_info(dataset, "raw", "")
        if res is not None:
            data.append((dataset, "raw", "", *res))

    for filter in ["filtlong", "kmrf"]:
        for k in range(13, 21, 2):
            for ratio in range(70, 100, 5):
                res = get_quast_info(dataset, filter, f".k{k}.r{ratio}")
                if res is not None:
                    data.append((dataset, filter, f".k{k}.r{ratio}", *res))

    return pandas.DataFrame(data, columns=["dataset", "filter", "kmer_size",
                                           "abundance", "nb_contigs",
                                           "nb_misassemblies", "unaligned_length",
                                           "genome_fraction", "nb_mismatches",
                                           "nb_indels", "largest_alignment",
                                           "NGA50"])

def get_quast_info(dataset, filter, params):
    path = f"filter/{dataset}/quast/{filter}/reads{params}/report.tsv"

    return utils.get_quast_info(path)
