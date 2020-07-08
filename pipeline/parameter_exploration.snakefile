###############################################################################
# Section simulate read                                                       #
###############################################################################
rule generate_bad_read:
    input:
        "references/CP028309.fasta"
        
    output:
        "reads/simulated_reads_{mean_id}.fasta"
        
    shell:
        "badread simulate --reference {input} --quantity 100x --identity {wildcards.mean_id},100,2 --error_model nanopore --seed 42 --start_adapter 0,0 --end_adapter 0,0 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 | seqtk seq -A - > {output}"

###############################################################################
# Section count                                                               #
###############################################################################
rule pcon_count:
    input:
        "{path}/{filename}.fasta"
        
    output:
        "{path}/{filename}.k{kmer_size}.pcon"
        
    resources:
        mem_mb = lambda wcd: (1 << (int(wcd.kmer_size) * 2 - 20)) * 8
        
    shell:
        "pcon count -i {input} -o {output} -k {wildcards.kmer_size}"

rule pcon_dump:
    input:
        "{path}/{filename}.k{kmer_size}.pcon"
        
    output:
        "{path}/{filename}.k{kmer_size}.a{abundance}.solid"

    resources:
        mem_mb = lambda wcd: ((1 << (int(wcd.kmer_size) * 2 - 20)) * 8) + (1 << (int(wcd.kmer_size) * 2 - 20))

    shell:
        "pcon dump -i {input} -a {wildcards.abundance} -s {output}"

###############################################################################
# Section evaluate error rate                                                 #
###############################################################################
minimap_parameter = {
    "ont": "-a -x map-ont",
    "pb": "-a -x map-pb",
}

reads_type = {
    "simulated_reads_95": "ont",
    "simulated_reads_98": "ont",
    "simulated_reads_99": "ont",
    #"SRR8556426": "ont",
    #"c_vartiovaarae": "ont",
    #"SRR8494940": "ont",
    #"SRR8494911": "pb",
    #"s_cerevisiae": "ont",
}

reference_path = {
    "simulated_reads_95": "references/CP028309.fasta",
    "simulated_reads_98": "references/CP028309.fasta",
    "simulated_reads_99": "references/CP028309.fasta",
    #"SRR8556426": "references/CP026549.fasta",
    #"c_vartiovaarae": None,
    #"SRR8494940": "references/CP028309.fasta",
    #"SRR8494911": "references/CP028309.fasta",
    #"s_cerevisiae": "references/GCA_002163515.fasta",
}

def build_index_name(filename):
    return "{}.{}.mmi".format(reference_path[filename].split(".")[0], reads_type[filename])

rule minimap_indexing:
    input:
        "{path}/{filename}.fasta"
        
    output:
        "{path}/{filename}.{reads_type}.mmi"
        
    params:
        options = lambda wcs: minimap_parameter[wcs.reads_type]
        
    shell:
        "minimap2 -t 1 {params.options} -d {output} {input}"

rule mapping:
    input:
        ref = lambda wcs: build_index_name(wcs.filename),
        reads = "{path}/{filename}.{parameter}fasta",
        
    output:
        "{path}/{filename}.{parameter}bam",
        
    params:
        reference = lambda wcs: reference_path[wcs.filename],
        options = lambda wcs: minimap_parameter[reads_type[wcs.filename]],

    threads:
        16
        
    wildcard_constraints:
        filename = "[^.]+",
        parameter = ".*",
        
    shell:
        "minimap2 -t {threads} {params.options} {input.ref} {input.reads} | samtools sort - > {output}"

rule generate_stat:
    input:
        "{path}/{filename}.{parameter}bam",
        
    output:
        "{path}/{filename}.{parameter}stats",
        
    wildcard_constraints:
        filename = "[^.]+",
        parameter = ".*",
        
    shell:
        "samtools stats -F 0x900 {input} | grep ^SN | cut -f 2- > {output}"

###############################################################################
# Section br call                                                             #
###############################################################################
rule br_genetic:
    input:
        filename = "reads/{filename}.fasta",
        solid = "references/CP028309.k{kmer_size}.a0.solid"
        
    output:
        "genetic_kmer/{filename}.k{kmer_size}.s{solidity}.fasta"

    resources:
        mem_mb = lambda wcd: 1 << (int(wcd.kmer_size) * 2 - 20)
        
    shell:
        "br -i {input.filename} -s {input.solid} -c 2 -o {output}"

rule br_read:
    input:
        filename = "reads/{filename}.fasta",
        solid = "reads/{filename}.k{kmer_size}.a{abundance}.solid"
        
    output:
        "read_kmer/{filename}.k{kmer_size}.a{abundance}.s{solidity}.fasta"

    resources:
        mem_mb = lambda wcd: 1 << (int(wcd.kmer_size) * 2 - 20)
        
    shell:
        "br -i {input.filename} -s {input.solid} -c 2 -o {output}"

###############################################################################
# Section global rules                                                        #
###############################################################################
rule genomic_kmer:
    input:
        "reads/simulated_reads_95.stats",
        ["reads/simulated_reads_95.k{}.pcon".format(k) for k in range(9, 21, 2)],
        ["genetic_kmer/simulated_reads_95.k{}.s{}.stats".format(k, s) for k in range(9, 21, 2) for s in range(1, 10)],
        "reads/simulated_reads_98.stats",
        ["reads/simulated_reads_98.k{}.pcon".format(k) for k in range(9, 21, 2)],
        ["genetic_kmer/simulated_reads_98.k{}.s{}.stats".format(k, s) for k in range(9, 21, 2) for s in range(1, 10)],
        "reads/simulated_reads_99.stats",
        ["reads/simulated_reads_99.k{}.pcon".format(k) for k in range(9, 21, 2)],
        ["genetic_kmer/simulated_reads_99.k{}.s{}.stats".format(k, s) for k in range(9, 21, 2) for s in range(1, 10)],

        
rule read_kmer:
    input:
        "reads/simulated_reads_95.stats",
        ["reads/simulated_reads_95.k{}.pcon".format(k) for k in range(9, 19, 2)],
        ["read_kmer/simulated_reads_95.k{}.a{}.s{}.stats".format(k, a, s) for k in range(13, 21, 2) for a in range(1, 16) for s in range(1, 10)],
        "reads/simulated_reads_98.stats",
        ["reads/simulated_reads_98.k{}.pcon".format(k) for k in range(9, 19, 2)],
        ["read_kmer/simulated_reads_98.k{}.a{}.s{}.stats".format(k, a, s) for k in range(13, 21, 2) for a in range(1, 16) for s in range(1, 10)],
        "reads/simulated_reads_99.stats",
        ["reads/simulated_reads_99.k{}.pcon".format(k) for k in range(9, 19, 2)],
        ["read_kmer/simulated_reads_99.k{}.a{}.s{}.stats".format(k, a, s) for k in range(13, 21, 2) for a in range(1, 16) for s in range(1, 10)],

              
# rule bacteria:
#     input:
#         "reads/SRR8494940.stats",
#         "reads/SRR8494911.stats",
#         "reads/SRR8556426.stats",
#         ["reads/SRR8494940.k{}.n8.pcon".format(k) for k in range(13, 19, 2)],
#         ["reads/SRR8494911.k{}.n8.pcon".format(k) for k in range(13, 19, 2)],
#         ["reads/SRR8556426.k{}.n8.pcon".format(k) for k in range(13, 19, 2)],
#         ["read_kmer/SRR8494940.k{}.n4.a{}.s{}.stats".format(k, a, s) for k in range(13, 21, 2) for a in range(1, 16) for s in range(1, 10)],
#         ["read_kmer/SRR8494911.k{}.n4.a{}.s{}.stats".format(k, a, s) for k in range(13, 21, 2) for a in range(1, 16) for s in range(1, 10)],
#         ["read_kmer/SRR8556426.k{}.n4.a{}.s{}.stats".format(k, a, s) for k in range(13, 21, 2) for a in range(1, 16) for s in range(1, 10)],
