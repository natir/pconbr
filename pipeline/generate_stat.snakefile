import os 
import sys
sys.path.append(os.path.abspath(os.getcwd()))

from pconbr.kmer_count import curve

include: "parameter_exploration.snakefile"

rule kmer_spectrum_simulated_reads_85:
    input:
        k13 = "reads/simulated_reads_85.k13.n8.pcon",
        k15 = "reads/simulated_reads_85.k15.n8.pcon",
        k17 = "reads/simulated_reads_85.k17.n8.pcon",
        k19 = "reads/simulated_reads_85.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads_85.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_count(input, output[0])

rule kmer_spectrum_simulated_reads_90:
    input:
        k13 = "reads/simulated_reads_90.k13.n8.pcon",
        k15 = "reads/simulated_reads_90.k15.n8.pcon",
        k17 = "reads/simulated_reads_90.k17.n8.pcon",
        k19 = "reads/simulated_reads_90.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads_90.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_count(input, output[0])

rule kmer_spectrum_simulated_reads_95:
    input:
        k13 = "reads/simulated_reads_95.k13.n8.pcon",
        k15 = "reads/simulated_reads_95.k15.n8.pcon",
        k17 = "reads/simulated_reads_95.k17.n8.pcon",
        k19 = "reads/simulated_reads_95.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads_95.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_count(input, output[0])
        
rule kmer_spectrum_e_coli_pb:
    input:
        k13 = "reads/SRR8494911.k13.n8.pcon",
        k15 = "reads/SRR8494911.k15.n8.pcon",
        k17 = "reads/SRR8494911.k17.n8.pcon",
        k19 = "reads/SRR8494911.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/e_coli_pb.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_count(input, output[0])

rule kmer_spectrum_e_coli_ont:
    input:
        k13 = "reads/SRR8494940.k13.n8.pcon",
        k15 = "reads/SRR8494940.k15.n8.pcon",
        k17 = "reads/SRR8494940.k17.n8.pcon",
        k19 = "reads/SRR8494940.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/e_coli_ont.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_count(input, output[0])
        
rule kmer_spectrum_s_pneumoniae:
    input:
        k13 = "reads/SRR8556426.k13.n8.pcon",
        k15 = "reads/SRR8556426.k15.n8.pcon",
        k17 = "reads/SRR8556426.k17.n8.pcon",
        k19 = "reads/SRR8556426.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/s_pneumoniae.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_count(input, output[0])

rule kmer_spectrum:
    input:
        "stats/kmer_spectrum/simulated_reads_90.csv",
        "stats/kmer_spectrum/simulated_reads_95.csv",
        "stats/kmer_spectrum/e_coli_ont.csv",
        "stats/kmer_spectrum/e_coli_pb.csv",
        "stats/kmer_spectrum/s_pneumoniae.csv",

rule kmer_spectrum_true_false_simulated_reads_85:
    input:
        true_k13 = "references/CP028309.k13.n8.pcon",
        k13 = "reads/simulated_reads_85.k13.n8.pcon",
        true_k15 = "references/CP028309.k15.n8.pcon",
        k15 = "reads/simulated_reads_85.k15.n8.pcon",
        true_k17 = "references/CP028309.k17.n8.pcon",
        k17 = "reads/simulated_reads_85.k17.n8.pcon",
        true_k19 = "references/CP028309.k19.n8.pcon",
        k19 = "reads/simulated_reads_85.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads_85_true_false.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_true_false(input, output[0])
        
rule kmer_spectrum_true_false_simulated_reads_90:
    input:
        true_k13 = "references/CP028309.k13.n8.pcon",
        k13 = "reads/simulated_reads_90.k13.n8.pcon",
        true_k15 = "references/CP028309.k15.n8.pcon",
        k15 = "reads/simulated_reads_90.k15.n8.pcon",
        true_k17 = "references/CP028309.k17.n8.pcon",
        k17 = "reads/simulated_reads_90.k17.n8.pcon",
        true_k19 = "references/CP028309.k19.n8.pcon",
        k19 = "reads/simulated_reads_90.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads_90_true_false.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_true_false(input, output[0])

rule kmer_spectrum_true_false_simulated_reads_95:
    input:
        true_k13 = "references/CP028309.k13.n8.pcon",
        k13 = "reads/simulated_reads_95.k13.n8.pcon",
        true_k15 = "references/CP028309.k15.n8.pcon",
        k15 = "reads/simulated_reads_95.k15.n8.pcon",
        true_k17 = "references/CP028309.k17.n8.pcon",
        k17 = "reads/simulated_reads_95.k17.n8.pcon",
        true_k19 = "references/CP028309.k19.n8.pcon",
        k19 = "reads/simulated_reads_95.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/simulated_reads_95_true_false.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_true_false(input, output[0])

rule kmer_spectrum_true_false_e_coli_pb:
    input:
        true_k13 = "references/CP028309.k13.n8.pcon",
        k13 = "reads/SRR8494911.k13.n8.pcon",
        true_k15 = "references/CP028309.k15.n8.pcon",
        k15 = "reads/SRR8494911.k15.n8.pcon",
        true_k17 = "references/CP028309.k17.n8.pcon",
        k17 = "reads/SRR8494911.k17.n8.pcon",
        true_k19 = "references/CP028309.k19.n8.pcon",
        k19 = "reads/SRR8494911.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/e_coli_pb_true_false.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_true_false(input, output[0])
        
        
rule kmer_spectrum_true_false_e_coli_ont:
    input:
        true_k13 = "references/CP028309.k13.n8.pcon",
        k13 = "reads/SRR8494940.k13.n8.pcon",
        true_k15 = "references/CP028309.k15.n8.pcon",
        k15 = "reads/SRR8494940.k15.n8.pcon",
        true_k17 = "references/CP028309.k17.n8.pcon",
        k17 = "reads/SRR8494940.k17.n8.pcon",
        true_k19 = "references/CP028309.k19.n8.pcon",
        k19 = "reads/SRR8494940.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/e_coli_ont_true_false.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_true_false(input, output[0])

rule kmer_spectrum_true_false_s_pneumoniae:
    input:
        true_k13 = "references/GCA_002163515.k13.n8.pcon",
        k13 = "reads/SRR8556426.k13.n8.pcon",
        true_k15 = "references/GCA_002163515.k15.n8.pcon",
        k15 = "reads/SRR8556426.k15.n8.pcon",
        true_k17 = "references/GCA_002163515.k17.n8.pcon",
        k17 = "reads/SRR8556426.k17.n8.pcon",
        true_k19 = "references/GCA_002163515.k19.n8.pcon",
        k19 = "reads/SRR8556426.k19.n8.pcon",
    output:
        "stats/kmer_spectrum/s_pneumoniae_true_false.csv"
    resources:
        mem_mb = lambda wcd: round((pow(2, 2 * 19 - 1)/2)/1000000)+10
    run:
        curve.generate_csv_true_false(input, output[0])

        
rule kmer_spectrum_true_false:
    input:
        "stats/kmer_spectrum/simulated_reads_85_true_false.csv",
        "stats/kmer_spectrum/simulated_reads_90_true_false.csv",
        "stats/kmer_spectrum/simulated_reads_95_true_false.csv",
        "stats/kmer_spectrum/e_coli_pb_true_false.csv",
        "stats/kmer_spectrum/e_coli_ont_true_false.csv",        
        "stats/kmer_spectrum/s_pneumoniae_true_false.csv",        
