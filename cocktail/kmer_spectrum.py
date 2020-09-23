import csv
import pandas

def get_kmer_spectrum(dataset, kmer_size):
    ref_count_kmer = csv.reader(open(f"count/{dataset}/pcon/reference.k{kmer_size}.csv"))
    read_count_kmer = csv.reader(open(f"count/{dataset}/pcon/reads.k{kmer_size}.csv"))

    all_kmer = [0 for _ in range(0, 256)]
    true_kmer = [0 for _ in range(0, 256)]
    false_kmer = [0 for _ in range(0, 256)]

    print("begin set of true kmer")
    ref_kmer = {kmer for (kmer, count) in ref_count_kmer}
    print("end set of true kmer")

    print("begin compute histograme")
    for (kmer, count) in read_count_kmer:
        count = int(count)

        all_kmer[count] += 1
        if kmer in ref_kmer:
            true_kmer[count] += 1
        else:
            false_kmer[count] += 1
    print("end compute histograme")
    
    df = pandas.DataFrame(data={"all": all_kmer, "true": true_kmer, "false": false_kmer}, index=range(0, 256))
    return df

def figure(df, column):
    return altair.Chart(df).mark_bar().encode(
        x=altair.X(column, bin=altair.Bin(maxbins=255), scale=altair.Scale(type='log', base=10)),
        y=altair.Y('count()'),
    )
