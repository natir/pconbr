include: "count.snakefile"
include: "parameter_exploration.snakefile"
include: "generate_stat.snakefile"

rule all:
    input:
        #rules.count_all.input,
        rules.genomic_kmer.input,
        rules.read_kmer.input,
        #rules.bacteria.input,
        #rules.kmer_spectrum.input,
        #rules.kmer_spectrum_true_false.input,
