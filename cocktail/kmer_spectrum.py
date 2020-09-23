import csv
import pandas

def get_kmer_spectrum(dataset, kmer_size):
    ref_count_kmer = csv.reader(open("count/{dataset}/pcon/reference.k{kmer_size}.csv"))
    read_count_kmer = csv.reader(open("data/{dataset}/pcon/reads.k{kmer_size}.csv"))

    all_kmer = [0 for _ in range(0, 256)]
    true_kmer = [0 for _ in range(0, 256)]
    false_kmer = [0 for _ in range(0, 256)]

    ref_kmer = {kmer for (kmer, count) in ref_count_kmer}

    for (kmer, count) in read_count_kmer:
        count = int(count)

        all_kmer[count] += 1
        if kmer in ref_kmer:
            true_kmer[count] += 1
        else:
            false_kmer[count] += 1

    pandas.DataFrame(data=(all_kmer, true_kmer, false_kmer), columns=["all", "true", "false"])
        
