#!/usr/bin/python2

import argparse
import subprocess
from Bio.Blast.Applications import NcbiblastxCommandline

parser = argparse.ArgumentParser(description="Application that removes sequences "
					"potentially contaminant sequences from a "
					"file using BLAST on a custom set of genomes")

parser.add_argument("-in", dest="infile", required=True, help="Provide the FastQ input "
					"file")
parser.add_argument("-r", dest="reference_database", nargs="*", help="Provide the "
					"path to the reference BLAST databases")
parser.add_argument("-c", dest="contaminant_database", nargs="*", help="Provide the "
					"path to the contaminant BLAST databases")
parser.add_argument("-o", dest="output_file", help="Provide the name of the output file")

arg = parser.parse_args()


def blast_wrapper(input_file, reference_database, contaminant_database, output_file):
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
	output_handle = open(output_file, "w")
	log_handle = open(input_file + ".log", "w")

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
				if result is not None:
					reference_results.append(result)

			for db in contaminant_database:
				result = blast_worker(temp_query, db)
				if result is not None:
					contaminant_results.append(result)

			# Triage
			result = triage(reference_results, contaminant_results)
			if result is True:
				output_handle.write(">%s\n%s\n" % (sequence_code, sequence))

			else:
				# The try statement accommodates the scenario where there are only blast
				#  hits on the contaminant database
				try:
					reference_ident = max([float(x[3]) for x in reference_results])
				except IndexError:
					reference_ident = "None"

				contaminant_ident = max([float(x[3]) for x in contaminant_results])
				log_handle.write("%s; %s; %s\n" % (sequence_code, reference_ident, contaminant_ident))

			temp_handle.close()
			# Removing temporary file
			subprocess.Popen(["rm temp.fas"], shell=True).wait()

	file_handle.close()
	output_handle.close()
	log_handle.close()


def triage(reference_result, contaminant_result):
	"""
	From two lists containing the BLAST results for the reference and contaminant
	databases, this function compares the evalue and identity of the query sequence
	:param reference_result: List containing the blast results of the reference databases
	:param contaminant_result: List containing the blast results of the contaminant
	databases
	:return: Either True, if there is no contamination or False, if there is contamination
	"""

	# In case the query sequence has blast hits on the contaminant but not in the
	# reference
	if reference_result is [] and contaminant_result is not []:
		return False

	# In case the query sequence has not hits in neither database, the benefit of doubt
	#  is given and the sequence is saved
	if reference_result is [] and contaminant_result is []:
		return True

	# If there are only blast hits on the reference database
	if reference_result is not [] and contaminant_result is []:
		return True

	#reference_evalues_max = max([float(x[2]) for x in reference_result])
	reference_ident_max = max([float(x[3]) for x in reference_result])

	#contaminant_evalues_max = max([float(x[2]) for x in contaminant_result])
	contaminant_ident_max = max([float(x[3]) for x in contaminant_result])

	if reference_ident_max > contaminant_ident_max:
		return True

	else:
		return False


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


def main():
	# Arguments
	input_file = arg.infie
	reference_databases = arg.reference_database
	contaminant_databases = arg.contaminant_database
	output_file = arg.output_file

	# Execution
	blast_wrapper(input_file, reference_databases, contaminant_databases, output_file)


main()