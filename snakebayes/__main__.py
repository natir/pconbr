#!/usr/bin/env python3

# standard import
import sys

# pip import
import snakemake
from bayes_opt import BayesianOptimization, UtilityFunction

# project import
from pconbr import identity

def evaluation(k, s):        
    return 1 - float(identity.get_error_rate(generate_filename(k, s)))

def normalize(k, s):
    return {"k": int(k) if int(k) % 2 else int(k) - 1, "s": int(s)}


def generate_filename(k, s):
    return "genetic_kmer/simulated_reads_90.k{}.n4.s{}.stats".format(k, int(s))


def main(args):
    if args is None:
        args = sys.argv[1:]

    optimizer = BayesianOptimization(
        f = None,
        pbounds = {"k": (7, 19), "s": (1, 15)},
        verbose=2,
        random_state=1,
    )
        
    utility = UtilityFunction(kind="ucb", kappa=2.5, xi=0.0)
    
    for _ in range(0, int(args[0])):
        parameters = normalize(**optimizer.suggest(utility))
        while  parameters.values() in optimizer.space:
            parameters = normalize(**optimizer.suggest(utility))
            
        print(parameters)
        print(generate_filename(**parameters))
        snakemake.snakemake("pipeline/parameter_exploration.snakefile", targets=[generate_filename(**parameters)])
        
        optimizer.register(params=parameters, target=evaluation(**parameters))


    return 0
        
if __name__ == "__main__":
    exit(main(sys.argv[1:]))
