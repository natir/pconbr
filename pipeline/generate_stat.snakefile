import os 
import sys
sys.path.append(os.path.abspath(os.getcwd()))

from pconbr.kmer_count import curve

include: "pconbr_evaluation.snakefile"

def kmer_spectrum_input(error):
    for k in range(config["kmer_begin"], config["kmer_end"], 2):
        yield f"reads/simulated_reads_{error}.k{k}.spectrum"

def kmer_spectrum_input_csv(error):
    for k in range(config["kmer_begin"], config["kmer_end"], 2):
        yield f"reads/simulated_reads_{error}.k{k}.a0.csv"
        
def true_kmer_input():
    for k in range(config["kmer_begin"], config["kmer_end"], 2):
        yield f"references/CP028309.k{k}.a0.csv"

        
rule kmer_spectrum_simulated_reads:
    input:
        spectrums = lambda wlc: kmer_spectrum_input(wlc.error)
        
    output:
        "stats/kmer_spectrum/simulated_reads_{error}.csv"

    wildcard_constraints:
        error = "\d+"
        
    run:
        curve.generate_csv_count(input, output[0])

rule kmer_spectrum:
    input:
        [f"stats/kmer_spectrum/simulated_reads_{error}.csv" for error in range(config["error_begin"], config["error_end"])]

        
rule kmer_spectrum_true_false_simulated_reads:
    input:
        true = lambda wlc: true_kmer_input(),
        reads = lambda wlc: kmer_spectrum_input_csv(wlc.error),
        
    output:
        "stats/kmer_spectrum/simulated_reads_{error}_true_false.csv"

    wildcard_constraints:
        error = "\d+"
        
    run:
        curve.generate_csv_true_false(input, output[0])
        
rule kmer_spectrum_true_false:
    input:
         [f"stats/kmer_spectrum/simulated_reads_{error}_true_false.csv" for error in range(config["error_begin"], config["error_end"])]
