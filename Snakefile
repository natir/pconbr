
configfile: "config.yaml"

include: "rules/utils.snk"
include: "rules/data.snk"
include: "rules/count.snk"
include: "rules/correct.snk"

wildcard_constraints:
    mean_id = '\d+',
    coverage = '\d+',
    kmer_size = '\d+',
    params = '.*',

rule all:
    input:
        rules.data_all.input,
        rules.count_all.input,
        rules.correct_all.input,
