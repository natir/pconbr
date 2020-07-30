
import math

def pcon_memory_usage(kmer_size, dump):
    if dump:
        if kmer_size >= 13:
            return math.ceil(
                1.1 *
                (
                    (1 << (kmer_size * 2 - 21)) +
                    (1 << (kmer_size * 2 - 24))
                )
            )
        else:
            return 10
    else:
        if kmer_size >= 11:
            return math.ceil(
                1.1 *
                (
                    1 << (kmer_size * 2 - 21)
                )
            )
        else:
            return 3


def br_memory_usage(kmer_size):
    if kmer_size >= 13:
        return math.ceil(
            1.1 *
            (
                1 << (kmer_size * 2 - 24)
            )
        )
    else:
        return 3


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
        mem_mb = lambda wcd: pcon_memory_usage(int(wcd.kmer_size), False)
        
    shell:
        "pcon count -i {input} -o {output} -k {wildcards.kmer_size}"

rule pcon_dump:
    input:
        "{path}/{filename}.k{kmer_size}.pcon"
        
    output:
        csv = "{path}/{filename}.k{kmer_size}.a{abundance}.csv",
        solid = "{path}/{filename}.k{kmer_size}.a{abundance}.solid",
        spectrum = "{path}/{filename}.k{kmer_size}.spectrum",

    resources:
        mem_mb = lambda wcd: pcon_memory_usage(int(wcd.kmer_size), True)

    shell:
        "pcon dump -i {input} -a {wildcards.abundance} -s {output.solid}, -S {output.spectrum}"


        
###############################################################################
# Section evaluate error rate                                                 #
###############################################################################
minimap_parameter = {
    "ont": "-a -x map-ont",
    "pb": "-a -x map-pb",
}

reads_type = {
    "simulated_reads_95": "ont",
    "simulated_reads_96": "ont",
    "simulated_reads_97": "ont",
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
    "simulated_reads_96": "references/CP028309.fasta",
    "simulated_reads_97": "references/CP028309.fasta",
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
        mem_mb = lambda wcd: br_memory_usage(int(wcd.kmer_size))
        
    shell:
        "br -i {input.filename} -s {input.solid} -c 2 -o {output}"

rule br_read:
    input:
        filename = "reads/{filename}.fasta",
        solid = "reads/{filename}.k{kmer_size}.a{abundance}.solid"
        
    output:
        "read_kmer/{filename}.k{kmer_size}.a{abundance}.s{solidity}.fasta"

    resources:
        mem_mb = lambda wcd: br_memory_usage(int(wcd.kmer_size))
        
    shell:
        "br -i {input.filename} -s {input.solid} -c 2 -o {output}"

###############################################################################
# Section global rules                                                        #
###############################################################################
def sim_genomic_input(error):
    yield f"reads/simulated_reads_{error}.stats"

    for k in range(config["kmer_begin"], config["kmer_end"], 2):
        yield f"reads/simulated_reads_{error}.k{k}.pcon"
        
    for k in range(config["kmer_begin"], config["kmer_end"], 2):
        for s in range(config["solidity_begin"], config["solidity_end"]):
            yield f"genetic_kmer/simulated_reads_{error}.k{k}.s{s}.stats"

def sim_reads_input(error):
    yield f"reads/simulated_reads_{error}.stats"

    for k in range(config["kmer_begin"], config["kmer_end"], 2):
        yield f"reads/simulated_reads_{error}.k{k}.pcon"
        
    for k in range(config["kmer_begin"], config["kmer_end"], 2):
        for a in range(config["abundance_begin"], config["abundance_end"]):
            for s in range(config["solidity_begin"], config["solidity_end"]):
                yield f"genetic_kmer/simulated_reads_{error}.k{k}.a{a}.s{s}.stats"

         
rule genomic_kmer:
    input:
        [sim_genomic_input(error) for error in range(95, 100)]
        
rule read_kmer:
    input:
        [sim_reads_input(error) for error in range(95, 100)]
              
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
