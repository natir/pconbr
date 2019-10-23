
###############################################################################
# Section genomic kmer                                                        #
###############################################################################
rule pcon_genetic:
    input:
        "{path}/{filename}.fasta"
    output:
        "{path}/{filename}.k{kmer_size}.exist"
    shell:
        " && ".join([
            "pcon count -i {input} -o {wildcards.path}/{wildcards.filename}.k{wildcards.kmer_size}.pcon -k {wildcards.kmer_size} -m 1",
           "pcon dump -i {wildcards.path}/{wildcards.filename}.k{wildcards.kmer_size}.pcon -o {output} -a 1 -m exist" 
        ])
        
rule br_genetic:
    input:
        filename="reads/{filename}.fasta",
        exist="references/CP028309.k{kmer_size}.exist"
    output:
        "genetic_kmer/{filename}.k{kmer_size}.s{solidity}.fasta"
    shell:
        "br -i {input.filename} -e {input.exist} -o {output} -s {wildcards.solidity}"

rule genomic_kmer:
    input:
        ["genetic_kmer/simulated_reads.k{}.s{}.fasta".format(k, s) for k in range(9, 21, 2) for s in range(1, 10)]

rule generate_bad_read:
    input:
        "references/CP028309.fasta"
    output:
        "reads/simulated_reads.fasta"
    shell:
        "badread simulate --reference {input} --quantity 100x --identity 95,100,2 --error_model nanopore --seed 42 --start_adapter 0,0 --end_adapter 0,0 | seqtk seq -A - > reads/simulated_reads.fasta"
    
