#!/usr/bin/env python3

import sys
import argparse

def main(args = None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    parser.add_argument("gfa")
    parser.add_argument("fasta")

    args = vars(parser.parse_args(args))
    
    gfa = args["gfa"]
    fasta = args["fasta"]
    
    with open(fasta, "w") as fh_out:
        with open(gfa) as fh:
            for line in fh:
                if line.startswith("S"):
                    _, tig_id, tig_seq, _ = line.split("\t")
                    fh_out.write(">{}\n{}\n".format(tig_id, tig_seq))


if __name__ == "__main__":
    main(sys.argv[1:])
