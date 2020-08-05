include: "count.snakefile"
include: "pconbr_evaluation.snakefile"
include: "generate_stat.snakefile"


rule pconbr_eval:
    input:
        rules.genomic_kmer.input,
        rules.read_kmer.input,

rule pconbr_spectrum:
    input:
        rules.kmer_spectrum.input,
        rules.kmer_spectrum_true_false.input,
        
        
rule all:
    input:
        rules.pconbr_eval.input,
        rules.pconbr_spectrum.input,
        #rules.count_all.input,
        #rules.bacteria.input,
