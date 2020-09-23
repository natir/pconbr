
configfile: "etc/config.yaml"

include: "rules/utils.snk"
include: "rules/data.snk"
include: "rules/count.snk"
include: "rules/correct.snk"
include: "rules/filter.snk"
include: "rules/polish.snk"

wildcard_constraints:
    mean_id = '\d+',
    coverage = '\d+',
    kmer_size = '\d+',
    abundance = '\d+',
    params = '.*',
    ratio = "\d+",
    type = '[^\.]+'
    
rule all:
    input:
        rules.data_all.input,
        rules.count_all.input,
        rules.count_ref.input,
        rules.correct_all.input,
        rules.filter_all.input,
        rules.polish_all.input,
