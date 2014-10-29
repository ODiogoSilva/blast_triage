#!/usr/bin/python2

import argparse
from Bio.Blast.Applications import NcbiblastxCommandline

parser = argparse.ArgumentParser(description="Application that removes sequences "
											"potentially contaminant sequences from a "
											"file using BLAST on a custom set of genomes")

arg = parser.parse_args()


def blast_wrapper(input_file, reference_database, contaminant_database):
	"""
	:param input_file: String with the input file name
	:param reference_database: List containing the path of the reference databases
	:param contaminant_database: List containing the path of the contaminant databases
	Main function of the BLAST phase of the script. It iterates over a sequence file in
	Fasta format and for each sequence it will BLAST it on the custom blast databases.
	The results of the blasts for the different databases will be parsed and compared
	to sort the sequences into contaminants or normal.
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
			temp_query = "temp.fas"
			temp_handle = open(temp_query, "w")
			temp_handle.write(">%s\n%s\n" % (sequence_code, sequence))

			# Creating temporary storage variables
			reference_results, contaminant_results = [], []

			# Executing BLAST for each specified database
			for db in reference_database:
				result = blast_worker(temp_query, db)
				reference_results.append(result)

			for db in contaminant_database:
				result = blast_worker(temp_query, db)
				contaminant_results.append(result)

	file_handle.close()


def blast_worker(query_file, blast_database):
	"""
	Executes biopython's BLAST on a custom database
	:param query_file: String with the input file name
	:param blast_database: String with database path
	:return:
	"""

	# Defining the blast command without specifying the output file name, redirects the
	# output to stdout, which can be parsed within the script
	blastx_cline = NcbiblastxCommandline(cmd="blastn", query=query_file,
										db=blast_database, evalue=0.001, outfmt='6 qseqid'
										' sseqid evalue pident', num_descriptions=1)

	stdout, sterr = blastx_cline()

	# Parse BLAST output
	try:
		vals = stdout.split()
		query_name = vals[0]
		subject_name = vals[1]
		evalue = vals[2]
		ident = vals[3]

		return query_name, subject_name, evalue, ident

	except IndexError:

		return None

