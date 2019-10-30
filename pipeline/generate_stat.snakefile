import sys
sys.path.append('..')

import pconbr.kmer_count.curve

rule kmer_spectrum_simulated_reads:
    input:
        k9  = "reads/simulated_reads.k9.n8.pcon",
        k11 = "reads/simulated_reads.k11.n8.pcon",
        k13 = "reads/simulated_reads.k13.n8.pcon",
        k15 = "reads/simulated_reads.k15.n8.pcon",
        k17 = "reads/simulated_reads.k17.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads.csv"
    run:
        pconbr.kmer_count.curve.generate_csv(input, output[0])

rule kmer_spectrum_e_coli_pb:
    input:
        k9  = "reads/SRR8494911.k9.n8.pcon",
        k11 = "reads/SRR8494911.k11.n8.pcon",
        k13 = "reads/SRR8494911.k13.n8.pcon",
        k15 = "reads/SRR8494911.k15.n8.pcon",
        k17 = "reads/SRR8494911.k17.n8.pcon",
    output:
        "stats/kmer_spectrum/e_coli_pb.csv"
    run:
        pconbr.kmer_count.curve.generate_csv(input, output[0])

rule kmer_spectrum_e_coli_ont:
    input:
        k9  = "reads/SRR8494940.k9.n8.pcon",
        k11 = "reads/SRR8494940.k11.n8.pcon",
        k13 = "reads/SRR8494940.k13.n8.pcon",
        k15 = "reads/SRR8494940.k15.n8.pcon",
        k17 = "reads/SRR8494940.k17.n8.pcon",
    output:
        "stats/kmer_spectrum/e_coli_ont.csv"
    run:
        pconbr.kmer_count.curve.generate_csv(input, output[0])
        
