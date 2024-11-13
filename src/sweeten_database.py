import genomicsqlite
import multiprocessing
from predict_genes import pyrodigal_manager
import numpy as np
from insertion_formatter import data_formatter
from build_database import genes_db


#Select genes from a database, 
class database_sweetener:
	def __init__(self, dbpath = None, threads = 1):
		self.dbpath = dbpath
		self.db = genes_db(self.dbpath, in_memory = False)
		
		self.genomes_to_process = None
		self.threads = threads
		
		self.indices = None
	
	def open(self):
		self.db.open()
		
	def close(self):
		self.db.close()
		
	
	def which_genomes(self):
		sql = "SELECT genome_id FROM genome_metadata WHERE has_predicted_genes = 0"
		self.genomes_to_process = []
		for gid in self.db.gsq_conn.execute(sql).fetchall():
			self.genomes_to_process.append(gid[0])
		
	def sweeten(self):
		#Load up the on-disk data
		self.open()
		self.db.get_metadata()
		self.which_genomes()
		
		if len(self.genomes_to_process) == 0:
			print("All genomes in the database already have genes predicted. There isn't anything to do.")
		else:
			print("Predicting genes for:", (len(self.genomes_to_process)), "genomes.")
		
			current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = self.db.give_current_indices()
			
			args = [(g, self.dbpath,) for g in self.genomes_to_process]
			
			formatter = data_formatter()
			
			#We're gonna make all changes in memory
			temp_conn = genes_db(dbpath = ":memory:", in_memory = True)
			temp_conn.open()
			temp_conn.initialize()
			temp_conn.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
			
			if self.threads == 1:
				print("Processing genomes...")
				for arg in args:
					results = parallel_process_genome(arg)
					file_name, genome_length, sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome, \
						translation_table, whole_genome_coding_bases, whole_genome_coding_dens = results
						
					current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = temp_conn.give_current_indices()
					#Update formatter indices
					formatter.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
					
					genome_insert, sequence_insert, contig_insert, gene_insert, genomes_to_drop = \
					formatter.add_next_data(file_name, genome_length, translation_table, whole_genome_coding_bases, whole_genome_coding_dens, \
					sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome)
					
					#May be worth allowing gene prediction to be skipped here, done later.
					#Retrieve new indices from the formatter
					current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = formatter.give_current_indices()			
					#Update database indices with new values
					temp_conn.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
					
					#Add the new data to database
					temp_conn.insert(genome_insert, None, contig_insert, gene_insert, drop_these_genomes = None)
					
					print(file_name, "complete!")
					
			else:
				print("Processing genomes in parallel...")
				pool = multiprocessing.Pool(self.threads)
				for results in pool.imap_unordered(parallel_process_genome, args):
					file_name, genome_length, sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome, \
					translation_table, whole_genome_coding_bases, whole_genome_coding_dens = results
					
					current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = temp_conn.give_current_indices()
					#Update formatter indices
					formatter.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
					
					genome_insert, sequence_insert, contig_insert, gene_insert, genomes_to_drop = \
					formatter.add_next_data(file_name, genome_length, translation_table, whole_genome_coding_bases, whole_genome_coding_dens, \
					sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome)
					
					#May be worth allowing gene prediction to be skipped here, done later.
					#Retrieve new indices from the formatter
					current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id = formatter.give_current_indices()			
					#Update database indices with new values
					temp_conn.set_indices(current_genome_index, current_contig_index, current_gene_index, current_genome_id, current_contig_id, current_gene_id)
					
					#Add the new data to database
					temp_conn.insert(genome_insert, None, contig_insert, gene_insert, drop_these_genomes = None)
					
					print(file_name, "complete!")
					
				pool.close()
				pool.join()
				
			
			genome_md = temp_conn.gsq_conn.execute("SELECT * FROM genome_metadata").fetchall()
			contig_md = temp_conn.gsq_conn.execute("SELECT * FROM contig_metadata").fetchall()
			gene_md = temp_conn.gsq_conn.execute("SELECT * FROM gene_metadata").fetchall()
			
			temp_conn.close()
			
			genomes_to_drop = [(g,) for g in self.genomes_to_process]
			
			self.db.insert(genome_md, None, contig_md, gene_md, drop_these_genomes = genomes_to_drop)
			self.db.write_changes()
			self.close()
		
		
#The version of this function for the sweeten script as opposed to main.
def parallel_process_genome(args):
	genome_id, dbpath = (args[0],), args[1]
	
	genome_sql = "SELECT * FROM genome_metadata WHERE genome_id = ?"
	contig_sql = "SELECT * FROM contig_metadata WHERE genome_id = ? ORDER BY contig_id"
	sequence_sql = "SELECT contig_id, twobit_dna(genome_sequence), loci_of_non_calls_replaced_with_A FROM genome_sequences WHERE genome_id = ? ORDER BY contig_id"
	
	conn = genomicsqlite.connect(dbpath, read_only = True)
	
	gen_tab =  conn.execute(genome_sql, genome_id).fetchall()
	cont_tab =  conn.execute(contig_sql, genome_id).fetchall()
	seq_tab =  conn.execute(sequence_sql, genome_id).fetchall()
	
	conn.close()
	
	gen_tab = gen_tab[0] #First row of what is always a 1-row selection
	
	contig_order = []
	contig_lengths = {}
	descriptions = {}
	
	rev_contig_idx = {}
	for row in cont_tab:
		contig_name = row[1]
		contig_order.append(contig_name)
		contig_id = row[3]
		rev_contig_idx[contig_id] = contig_name
		contig_lengths[contig_name] = row[4]
		descriptions[contig_name] = row[2]
		
	contig_name = None
	contig_id = None
	
	genome_dict = {}
	for sequence_tuple in seq_tab:
		contig_id = sequence_tuple[0]
		contig_name = rev_contig_idx[contig_id]
		
		#Re-insert non-call characters if needed.
		if sequence_tuple[2] is not None:
			a_loci = np.frombuffer(sequence_tuple[2], dtype = np.int32)
			#Positions of the non-call loci on the reverse strand
			seq = list(sequence_tuple[1])
			for locus in a_loci:
				seq[locus] = "N"
			seq = ''.join(seq)
		else:
			seq = sequence_tuple[1]
			
		genome_dict[contig_name] = seq
	
	predictor = pyrodigal_manager()
	
	genes_by_genome, translation_table, whole_genome_coding_bases, whole_genome_coding_dens = predictor.predict_genes(genome_dict)
	
	genome_file = gen_tab[0]
	genome_length = gen_tab[3]
	
	return [genome_file, genome_length, contig_order, None, contig_lengths, descriptions, None, genes_by_genome, translation_table, whole_genome_coding_bases, whole_genome_coding_dens]
	
	

