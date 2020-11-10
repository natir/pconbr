#!/bin/bash 

set -e
set -u
set -o pipefail
set -x

function main() {
    if [ "$#" -ne 5 ]; then
	usage
	exit
    fi

    coverage=$1
    thread=$2
    ref=$3
    identity=$4
    output=$5
    pid=$BASHPID
    
    cov_per_thread=$(($coverage / $thread))
    for i in $(seq 1 $thread); do
	run ${cov_per_thread} ${ref} ${identity} ${pid} $i &
    done

    wait

    cat temp_${pid}_*.fastq | seqtk seq -A - | sed -e 's/^\(>[^[:space:]]*\).*/\1/' > ${output}
    rm temp_${pid}_*.fastq

    cat temp_${pid}_*.log
    rm temp_${pid}_*.log
}

function run() {
    coverage=$1
    ref=$2
    identity=$3
    output=temp_${4}_${5}.fastq
    log=temp_${4}_${5}.log
    seed=$((42 + ${5}))
    
    badread simulate --reference ${ref} --quantity ${coverage}x --identity ${identity},100,5 --error_model nanopore --seed ${seed} 2> ${log} > ${output}
}

function usage() {
    echo "multi-thead-badread.sh usage:"
    echo "multi-thead-badread.sh {coverage} {thread} {reference} {identity} {output}"
}

main $@

