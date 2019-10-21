#!/bin/bash

mkdir -p references

# Streptococcus pneumoniae
curl https://www.ebi.ac.uk/ena/browser/api/fasta/CP026549?download=true > references/CP026549.fasta

# Candida vartiovaarae


# E. coli
curl https://www.ebi.ac.uk/ena/browser/api/fasta/CP028309?download=true > references/CP028309.fasta

# Saccharomyces cerevisiae W303
curl https://www.ebi.ac.uk/ena/browser/api/fasta/CM007964.1,CM007965.1,CM007966.1,CM007967.1,CM007968.1,CM007969.1,CM007970.1,CM007971.1,CM007972.1,CM007973.1,CM007974.1,CM007975.1,CM007976.1,CM007977.1,CM007978.1,CM007979.1,CM007980.1,CM007981.1?download=true > references/GCA_002163515.fasta 
