rule pcon_count:
    input:
        "{path}/{filename}.fasta"
    output:
        "{path}/{filename}.k{kmer_size}.n{nb_bit}.pcon"
    shell:
        "pcon count -i {input} -o {output} -k {wildcards.kmer_size} -m 1 -n {wildcards.nb_bit}",
        
        
rule pcon_dump:
    input:
        "{path}/{filename}.k{kmer_size}.n{nb_bit}.pcon"
    output:
        "{path}/{filename}.k{kmer_size}.n{nb_bit}.a{abundance}.exist"
    shell:
        "pcon dump -i {input} -o {output} -a {wildcards.abundance} -m exist"


rule generate_bad_read:
    input:
        "references/CP028309.fasta"
    output:
        "reads/simulated_reads.fasta"
    shell:
        "badread simulate --reference {input} --quantity 100x --identity 95,100,2 --error_model nanopore --seed 42 --start_adapter 0,0 --end_adapter 0,0 | seqtk seq -A - > reads/simulated_reads.fasta"

###############################################################################
# Section evaluate read                                                       #
###############################################################################
minimap_parameter = {
    "simulated_reads": "-a -x map-ont",
    "SRR8556426": "-a -x map-ont",
    "c_vartiovaarae": "-a -x map-ont",
    "SRR8494940": "-a -x map-ont",
    "SRR8494911": "-a -x map-pb",
    "s_cerevisiae": "-a -x map-ont",
}

reference_path = {
    "simulated_reads": "references/CP028309.fasta",
    "SRR8556426": "references/CP026549.fasta",
    "c_vartiovaarae": None,
    "SRR8494940": "references/CP028309.fasta",
    "SRR8494911": "references/CP028309.fasta",
    "s_cerevisiae": "references/GCA_002163515.fasta",
}

rule mapping:
    input:
        "{path}/{filename}.{parameter}fasta",
    output:
        "{path}/{filename}.{parameter}bam",
    params:
        reference=lambda wcs: reference_path[wcs.filename],
        options=lambda wcs: minimap_parameter[wcs.filename],
    wildcard_constraints:
        filename="[^.]+",
        parameter=".*",
    shell:
        "minimap2 {params.options} {params.reference} {input} | samtools sort - > {output}"
        
rule generate_stat:
    input:
        "{path}/{filename}.{parameter}bam",
    output:
        "{path}/{filename}.{parameter}stats",
    wildcard_constraints:
        filename="[^.]+",
        parameter=".*",
    shell:
        "samtools stats -F 0x900 {input} | grep ^SN | cut -f 2- > {output}"
        
###############################################################################
# Section genomic kmer                                                        #
###############################################################################
rule br_genetic:
    input:
        filename="reads/{filename}.fasta",
        exist="references/CP028309.k{kmer_size}.n{nb_bit}.a1.exist"
    output:
        "genetic_kmer/{filename}.k{kmer_size}.n{nb_bit}.s{solidity}.fasta"
    shell:
        "br -i {input.filename} -e {input.exist} -o {output} -s {wildcards.solidity}"

rule genomic_kmer:
    input:
        ["genetic_kmer/simulated_reads.k{}.n8.pcon".format(k) for k en range(9, 19, 2)],
        ["genetic_kmer/simulated_reads.k{}.n4.s{}.stats".format(k, s) for k in range(9, 19, 2) for s in range(1, 10)],

###############################################################################
# Section genomic kmer                                                        #
###############################################################################
rule br_read:
    input:
        filename="reads/{filename}.fasta",
        exist="reads/{filename}.k{kmer_size}.a{abundance}.exist"
    output:
        "read_kmer/{filename}.k{kmer_size}.a{abundance}.s{solidity}.fasta"
    shell:
        "br -i {input.filename} -e {input.exist} -o {output} -s {wildcards.solidity}"

rule read_kmer:
    input:
        ["genetic_kmer/simulated_reads.k{}.n8.pcon".format(k) for k en range(9, 19, 2)],
        ["read_kmer/simulated_reads.k{}.n4.a{}.s{}.stats".format(k, a, s) for k in range(13, 19, 2) for a in range(1, 15) for s in range(1, 10)]

rule bacteria:
    input:
        "reads/SRR8494940.stats",
        "reads/SRR8494911.stats",
        ["read_kmer/SRR8494940.k{}.n8.pcon".format(k) for k in range(13, 19, 2)],
        ["read_kmer/SRR8494911.k{}.n8.pcon".format(k) for k in range(13, 19, 2)],
        ["read_kmer/SRR8494940.k{}.n4.a{}.s{}.stats".format(k, a, s) for k in range(13, 19, 2) for a in range(1, 15) for s in range(1, 10)],
        ["read_kmer/SRR8494911.k{}.n4.a{}.s{}.stats".format(k, a, s) for k in range(13, 19, 2) for a in range(1, 15) for s in range(1, 10)],
   
