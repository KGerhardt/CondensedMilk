import sys
import os
import argparse

import multiprocessing

#Build a database
from load_genome import genome_loader
from dna_encoder import dna_encoder
from predict_genes import pyrodigal_manager
from insertion_formatter import data_formatter
from build_database import genes_db

#Merge databases
from combine_databases import db_combiner

#Predict genes for any genomes that don't have them
from sweeten_database import database_sweetener

#Retrieve genome, gene info and output to standard file formats
from extract_db_info import db_extractor

#Functions for encoding genomes and producing/adding to a database
def parallel_process_one_genome(args):
	genome_file, run_gene_prediction = args[0], args[1]
	
	enc = dna_encoder()
	gl = genome_loader(genome_file)
	
	#Load genome seqeunce, find loci of non-calls and replace them with "A" characters
	genome_length, sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices = gl.load_genome()
	
	if run_gene_prediction:
		mn = pyrodigal_manager()
		#Run gene prediction
		genes_by_genome, translation_table, whole_genome_coding_bases, whole_genome_coding_dens = mn.predict_genes(sequences)
	else:
		genes_by_genome, translation_table, whole_genome_coding_bases, whole_genome_coding_dens = None, None, None, None
	
	#Encode text seq to twobit sequences using genomicsqlite

	sequences = enc.compose_twobit(sequences, sequence_order)
	
	return [genome_file, genome_length, sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome, translation_table, whole_genome_coding_bases, whole_genome_coding_dens]
	
def genomes_to_db(genome_directory, output_database, num_threads = 1, in_memory = False, predict_genes = False):
	files = os.listdir(genome_directory)
	full_paths = [os.path.normpath(genome_directory + "/" + f) for f in files]
	
	formatter = data_formatter()
	
	db = genes_db(output_database, in_memory)
	db.open()
	db.initialize()
	db.get_metadata()

	#Don't bother with opening subprocesses if there's only one thread
	if num_threads == 1:	
		#Format and insert should happen in the main process, genome reading and gene prediction in child processes 
		for genome_file in full_paths:
			print("Processing", genome_file)
			
			file_name, genome_length, sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome, \
			translation_table, whole_genome_coding_bases, whole_genome_coding_dens = parallel_process_one_genome([genome_file, predict_genes])
			
			#Collect database indices
			current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = db.give_current_indices()
			#Update formatter indices
			formatter.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
			
			genome_insert, sequence_insert, contig_insert, gene_insert, genomes_to_drop = \
			formatter.add_next_data(file_name, genome_length, translation_table, whole_genome_coding_bases, whole_genome_coding_dens, \
			sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome)
			
			#May be worth allowing gene prediction to be skipped here, done later.
			#Retrieve new indices from the formatter
			current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = formatter.give_current_indices()			
			#Update database indices with new values
			db.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
			
			#Add the new data to database
			db.insert(genome_insert, sequence_insert, contig_insert, gene_insert, genomes_to_drop)
			
			print(genome_file, "complete!")
		
	else:
		print("Running genomes in parallel...")
	
		parallel_args = [(g, predict_genes) for g in full_paths]
	
		pool = multiprocessing.Pool(num_threads)
		
		for results in pool.imap_unordered(parallel_process_one_genome, parallel_args):
			file_name, genome_length, sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome, \
			translation_table, whole_genome_coding_bases, whole_genome_coding_dens = results
			
			#Collect database indices
			current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = db.give_current_indices()
			#Update formatter indices
			formatter.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
			
			genome_insert, sequence_insert, contig_insert, gene_insert, genomes_to_drop = \
			formatter.add_next_data(file_name, genome_length, translation_table, whole_genome_coding_bases, whole_genome_coding_dens, \
			sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome)
			
			#Retrieve new indices from the formatter
			current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = formatter.give_current_indices()			
			#Update database indices with new values
			db.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
			
			#Add the new data to database
			db.insert(genome_insert, sequence_insert, contig_insert, gene_insert, genomes_to_drop)
			
			print(file_name, "complete!")
		
		pool.close()
		pool.join()
	
	db.write_changes()
	db.close()

	print("Condensed Milk complete!")
	
	return 1

def condense_options():
	parser = argparse.ArgumentParser(description='Add genomes to a condensed milk database')
	parser.add_argument('-g', '--genomes', dest='genomes', default = None,
						help='Directory of prokaryotic genome files in FASTA format. These genomes will be compressed, have genes predicted for them using Prodigal, and be added to a database.')
	parser.add_argument('-db', '--database', dest='database', default = None,
						help='Path to a database for condensed milk. If it exists, genomes will be added to the database. If not, the database will be created at this location.')
	parser.add_argument('--memory', dest = 'in_memory', action = 'store_true',
						help='Create a database in memory for processing, then write any changes to disk only at the end of this program. This is faster, has fewer I/O ops than writing directly to disk when processing many genomes, but takes more memory. Use this option if processing many genomes or adding to a very large database on a system without much RAM.')
	parser.add_argument('-p', '--threads', dest = 'threads', default = 1, type = int,
						help='Number of parallel processes to use for predicting genes. Default 1. Parallelized per-file, so using more threads than you have genomes is pointless.')
	
	parser.add_argument('--sweeten', dest = 'predict_genes', action = 'store_true',
						help='Run gene prediction for each genome in the --genomes directory.')


	args = parser.parse_known_args()
	return parser, vars(args[0])

#function for adding genes to a database which is missing gene calls for some/all genomes.
def sweeten_options():
	parser = argparse.ArgumentParser(description='Predict genes for any genomes in a condensed milk database without such predictions and add them to the database. You can also use this module to check if all genomes in the database have genes.')
	parser.add_argument('-db', '--database', dest='database', default = None,
						help='Path to a database for condensed milk. If it exists, genomes will be added to the database. If not, the database will be created at this location.')
	parser.add_argument('-p', '--threads', dest = 'threads', default = 1, type = int,
						help='Number of parallel processes to use for predicting genes. Default 1. Parallelized per-file, so using more threads than you have genomes is pointless.')

	args = parser.parse_known_args()
	return parser, vars(args[0])
	
#Functions for merging two databases - a TODO
def simmer_options():
	parser = argparse.ArgumentParser(description='Combine two Condensed Milk databases.')
	parser.add_argument('--donor', dest='donor', default = None,
						help='Path to a Condensed Milk database. The genomes in this database will be added to the database specified nby --recipient as new genomes in that database. The donor database will not be modified.')
	parser.add_argument('--recipient', dest = 'recipient', default = None,
						help='')

	args = parser.parse_known_args()
	return parser, vars(args[0])

#Functions for loading, writing data from a database

def retrieve_options():
	parser = argparse.ArgumentParser(description='Retrieve data from a condensed milk database')
	
	#Basic info
	parser.add_argument('-db', '--database', dest='database', default = None,
						help='Path to an EXISTING database created by Condensed Milk.')
						
	parser.add_argument('-o', '--output_dir', dest='out', default = 'baked_goods',
						help='An output directory in which genomes, genes, and metadata from this database will be placed.')
			

	#gzip compression of outputs
	parser.add_argument('--compress', dest='compress', action = 'store_true',
						help='gzip any output files from this call. Files are not compressed by default.')
	
	#Query for only one, a few, or many genomes instead of outputting the whole database
	parser.add_argument('--these_genomes', dest = 'genome_list', default = None,
						help = 'Comma-sep list of genomes to retrieve. Use names exactly as written in the metadata files. If this and --these_genomes_file are both absent, all genomes will be retrieved.')
						
	parser.add_argument('--these_genomes_file', dest = 'genome_list_file', default = None,
						help = 'File containing one genome name per row to retrieve. Use names exactly as written in the metadata files. If this and --these_genomes are both absent, all genomes will be retrieved.')
						
	#Extract one or several kinds of data for the genomes indicated. Default if all are false is to do everything.
	parser.add_argument('--describe_JSON', dest = 'do_json', action = 'store_true',
						help = 'Retrieve genome/MAG, contig, and gene metadata from this database and output as a single JSON-file.')
						
	parser.add_argument('--describe_text', dest = 'do_flat_files', action = 'store_true',
						help = 'Retrieve metadata from this database and output as three tab-separated files, one for genome/MAG level data, one for contig level data, and one for gene level data.')
						
	parser.add_argument('--describe_text_single', dest = 'do_single_flat', action = 'store_true',
						help = 'Retrieve genome/MAG, contig, and gene metadata from this database and output as one tab-separated file with as much redundancy as needed to make this happen.')
						
	parser.add_argument('--genome_to_fasta', dest = 'genome_seqs', action = 'store_true',
						help = 'Retrieve genome sequences and recreate the original genome sequence files used as inputs to this database.')		
						
	parser.add_argument('--protein_nt', dest = 'protein_nt', action = 'store_true',
						help = 'Output proteins for each genome into per-genome FASTA files of the nucleotide sequences for each protein')		
						
	parser.add_argument('--protein_aa', dest = 'protein_aa', action = 'store_true',
						help = 'Output proteins for each genome into per-genome FASTA files of the amino acid sequences for each protein')
						
	
	args = parser.parse_known_args()
	return parser, vars(args[0])

def run_retrieve(database, outdir, compress = False, genome_list = None, genome_file = None, json_meta = True, text_meta = True, single_text_meta = True, genome_fasta = True, protein_nt = True, protein_aa = True):
	outdir = os.path.normpath(outdir)
		
	mn = db_extractor(dbpath = database, output_dir = outdir, compress = compress)
	mn.open()
	
	mn.prep_files()
	mn.extract_database_description()
	mn.write_json()
	mn.write_text_triplet()
	mn.write_merged_text()
	mn.iterate_genomes(True, True, True)
	
	mn.close()

def main():
	ok_modules = ["condense", "sweeten", "simmer", "retrieve"]
	if len(sys.argv) < 2:
		writeout = '/'.join(ok_modules)
		print("You need to supply an action to condensed milk. The action should be in the second position of the call, e.g.:")
		print('\t"python3 condensed_milk [{mods}]"'.format(mods=writeout))
		sys.exit()
	else:
		
		module = sys.argv[1]
		if module not in ok_modules:
			print("Module selection '" + str(module) + "' not recognized!")
			print("Module must be one of:")
			for m in ok_modules:
				print("\t", m)
			
			sys.exit()
			
		else:
			if module == "condense":
				parser, args = condense_options()
				
				input_directory = args['genomes']
				output_db = args['database']
				num_threads = args['threads']
				in_mem = args['in_memory']
				do_genes = args['predict_genes']
				
				is_ok = True
				if input_directory is None:
					print("Condensed Milk needs a directory of input genomes in FASTA format!")
					is_ok = False
				if output_db is None:
					print("Condensed Milk needs a path to an output database!")
					
				if is_ok:
					genomes_to_db(genome_directory = input_directory, output_database = output_db, num_threads = num_threads, in_memory = in_mem, predict_genes = do_genes)
				else:
					print("")
					parser.print_help()
					print("")
			
			if module == "sweeten":
				parser, args = sweeten_options()
				database = args['database']
				num_threads = args['threads']
				
				is_ok = True
				if database is None:
					print("Condensed Milk needs a database to add genes to!")
					print("")
					is_ok = False
					parser.print_help()
				else:
					if not os.path.exists(database):
						print("Cannot find a database at", database)
						is_ok = False
						
					if is_ok:
						mn = database_sweetener(dbpath = database, threads = num_threads)
						mn.sweeten()
					else:
						print("")
						parser.print_help()
						print("")
			
			if module == "simmer":
				parser, args = simmer_options()
				donor = args['donor']
				recip = args['recipient']
				
				is_ok = True
				if donor is None or recip is None:
					print("Condensed Milk needs both a donor and a recipient database to combine!")
					is_ok = False
				if is_ok:
					print("")
					print("Adding the contents of", donor, "to", recip)
					mn = db_combiner(donor = donor, recipient = recip)
					mn.run_combine()
					print("Done!")
			
			if module == "retrieve":
				parser, args = retrieve_options()
				
				db = args['database']
				out = args['out']
				compress = args['compress']

				is_ok = True
				if db is None:
					print("You must supply a database to Condensed Milk!")
					is_ok = False
				else:
					if not os.path.exists(db):
						print("Database", db, "not found!")
						is_ok = False
					
					if is_ok:
						run_retrieve(db, out, compress)
		
		


main()