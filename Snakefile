
configfile: "config.yaml"

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
    ratio = "0\.\d+",
    type = "([^\d\.]+(\.e\d\d)?(\.c\d\d)?)",
    
rule all:
    input:
        rules.data_all.input,
        rules.count_all.input,
        rules.correct_all.input,
        rules.filter_all.input,
        rules.polish_all.input,
