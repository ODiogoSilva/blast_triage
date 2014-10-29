#!/usr/bin/python2

import argparse
from Bio.Blast.Applications import NcbiblastxCommandline

parser = argparse.ArgumentParser(description="Application that removes sequences "
											"potentially contaminant sequences from a "
											"file using BLAST on a custom set of genomes")

arg = parser.parse_args()

