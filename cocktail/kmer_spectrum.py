import gzip

import altair
import pandas

def get_kmer_spectrum(dataset, kmer_size):
    ref_count_kmer = read_pcon_bin(f"count/{dataset}/pcon/reference.k{kmer_size}.pcon")
    read_count_kmer = read_pcon_bin(f"count/{dataset}/pcon/reads.k{kmer_size}.pcon")

    all_kmer = [0 for _ in range(0, 256)]
    true_kmer = [0 for _ in range(0, 256)]
    false_kmer = [0 for _ in range(0, 256)]

    ref_kmer = {kmer for (kmer, count) in ref_count_kmer if count > 0}

    for (kmer, count) in read_count_kmer:        
        #print(kmer, count)
        count = int(count)

        all_kmer[count] += 1
        if kmer in ref_kmer:
            true_kmer[count] += 1
        else:
            false_kmer[count] += 1

    data = [("all", i, all_kmer[i])for i in range(1, 256)]
    data += [("true", i, true_kmer[i])for i in range(1, 256)]
    data += [("false", i, false_kmer[i])for i in range(1, 256)]

    df = pandas.DataFrame(data=data, columns=["type", "abundance", "count"])

    return df

def figure(df):
    return altair.Chart(df).mark_area(opacity=0.7).encode(
        x=altair.X("abundance", scale=altair.Scale(domain=(1, 255))),
        y=altair.Y("count", scale=altair.Scale(type="log", base=10)),
        #column="type",
        color="type"
    )

def read_pcon_bin(path):
    fh = gzip.open(path)

    k = fh.read(1)
    
    while True:
        buffer = fh.read(8*1_000_000)
        if buffer == b"":
            break
        
        for i in range(0, len(buffer)):
            count = buffer[i]
        
            yield (i, count)
