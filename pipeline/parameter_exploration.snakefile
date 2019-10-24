rule pcon:
    input:
        "{path}/{filename}.fasta"
    output:
        "{path}/{filename}.k{kmer_size}.a{abundance}.exist"
    shell:
        " && ".join([
            "pcon count -i {input} -o {wildcards.path}/{wildcards.filename}.k{wildcards.kmer_size}.a{wildcards.abundance}.pcon -k {wildcards.kmer_size} -m 1",
           "pcon dump -i {wildcards.path}/{wildcards.filename}.k{wildcards.kmer_size}.a{wildcards.abundance}.pcon -o {output} -a {wildcards.abundance} -m exist" 
        ])


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
    "s_pneumoniae": "-a -x map-ont",
    "c_vartiovaarae": "-a -x map-ont",
    "e_coli_ont": "-a -x map-ont",
    "e_coli_pb": "-a -x map-pb",
    "s_cerevisiae": "-a -x map-ont",
}

reference_path = {
    "simulated_reads": "references/CP028309.fasta",
    "s_pneumoniae": "references/CP026549.fasta",
    "c_vartiovaarae": None,
    "e_coli_ont": "references/CP028309.fasta",
    "e_coli_pb": "references/CP028309.fasta",
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
        exist="references/CP028309.k{kmer_size}.a1.exist"
    output:
        "genetic_kmer/{filename}.k{kmer_size}.s{solidity}.fasta"
    shell:
        "br -i {input.filename} -e {input.exist} -o {output} -s {wildcards.solidity}"

rule genomic_kmer:
    input:
        ["genetic_kmer/simulated_reads.k{}.s{}.stats".format(k, s) for k in range(9, 19, 2) for s in range(1, 10)]

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
        ["read_kmer/simulated_reads.k{}.a{}.s{}.stats".format(k, a, s) for k in range(9, 19, 2) for a in range(1, 10) for s in range(1, 10)]

