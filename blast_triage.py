#!/usr/bin/python2

import argparse
from Bio.Blast.Applications import NcbiblastxCommandline

parser = argparse.ArgumentParser(description="Application that removes sequences "
											"potentially contaminant sequences from a "
											"file using BLAST on a custom set of genomes")

arg = parser.parse_args()


def blast_wrapper(input_file, database_list):
	"""
	:param input_file: String with the input file name
	:param database_list: List containing the name of the databases
	Main function of the BLAST phase of the script. It iterates over a sequence file in
	Fasta format and for each sequence it will BLAST it on the custom blast databases,
	creating an individual file containing the best hit for each one. Then, it will
	parse the BLAST output files and create an instance of the Triage class containing
	the attributes of each best hit.
	:return:
	"""

	file_handle = open(input_file)

	# Start iterating over input file
	for line in file_handle:

		if line.startswith("@"):
			# Retrieving sequence and sequence name
			sequence_code = line[1:].strip()
			sequence = file_handle.next().strip()

			# Creating temporary input file or blast
			temp_handle = open("temp.fas", "w")
			temp_handle.write(">%s\n%s\n" % (sequence_code, sequence))

			# Executing BLAST for each specified database

	file_handle.close()