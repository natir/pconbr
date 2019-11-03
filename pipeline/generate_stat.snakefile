import os 
import sys
sys.path.append(os.path.abspath(os.getcwd()))

from pconbr.kmer_count import curve

include: "parameter_exploration.snakefile"

rule kmer_spectrum_simulated_reads:
    input:
        k13 = "reads/simulated_reads.k13.n8.pcon",
        k15 = "reads/simulated_reads.k15.n8.pcon",
        k17 = "reads/simulated_reads.k17.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads.csv"
    run:
        curve.generate_csv_count(input, output[0])

rule kmer_spectrum_e_coli_pb:
    input:
        k13 = "reads/SRR8494911.k13.n8.pcon",
        k15 = "reads/SRR8494911.k15.n8.pcon",
        k17 = "reads/SRR8494911.k17.n8.pcon",
    output:
        "stats/kmer_spectrum/e_coli_pb.csv"
    run:
        curve.generate_csv_count(input, output[0])

rule kmer_spectrum_e_coli_ont:
    input:
        k13 = "reads/SRR8494940.k13.n8.pcon",
        k15 = "reads/SRR8494940.k15.n8.pcon",
        k17 = "reads/SRR8494940.k17.n8.pcon",
    output:
        "stats/kmer_spectrum/e_coli_ont.csv"
    run:
        curve.generate_csv_count(input, output[0])
        
rule kmer_spectrum_s_pneumoniae:
    input:
        k13 = "reads/SRR8556426.k13.n8.pcon",
        k15 = "reads/SRR8556426.k15.n8.pcon",
        k17 = "reads/SRR8556426.k17.n8.pcon",
    output:
        "stats/kmer_spectrum/s_pneumoniae.csv"
    run:
        curve.generate_csv_count(input, output[0])

rule kmer_spectrum:
    input:
        "stats/kmer_spectrum/simulated_reads.csv",
        "stats/kmer_spectrum/e_coli_ont.csv",
        "stats/kmer_spectrum/e_coli_pb.csv",
        "stats/kmer_spectrum/s_pneumoniae.csv",

rule kmer_spectrum_true_false_simulated_reads:
    input:
        true_k13 = "references/CP028309.k13.n8.pcon",
        k13 = "reads/simulated_reads.k13.n8.pcon",
        true_k15 = "references/CP028309.k15.n8.pcon",
        k15 = "reads/simulated_reads.k15.n8.pcon",
        true_k17 = "references/CP028309.k17.n8.pcon",
        k17 = "reads/simulated_reads.k17.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads_true_false.csv"
    run:
        curve.generate_csv_true_false(input, output)
        
rule kmer_spectrum_true_false:
    input:
        "stats/kmer_spectrum/simulated_reads_true_false.csv",
