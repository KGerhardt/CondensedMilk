import genomicsqlite
import os
	
class genes_db:
	def __init__(self, dbpath, in_memory = False):
		self.dbpath = dbpath
		self.gsq_conn = None
		self.in_mem = in_memory
		
		self.genome_index = None
		self.reverse_genome_index = None
		self.contig_index = None
		self.reverse_contig_index = None
		self.gene_index = None
		self.reverse_gene_index = None
		self.gene_prediction_record = None
		
		self.current_genome_index = 0
		self.current_contig_index = 0
		self.current_gene_index = 0
		
	def open(self):
		if self.gsq_conn is None:
			if self.in_mem:
				self.gsq_conn = genomicsqlite.connect(":memory:")
			else:
				self.gsq_conn = genomicsqlite.connect(self.dbpath)
		
	def write_changes(self, path = None):
		#Path exists as an option here for use as part of a python module
		if self.in_mem:	
			if path is None:
				usepath = self.dbpath
			else:
				usepath = path
			
			if os.path.exists(usepath):
				print("Writing changes to existing database...")
				genome_md = self.gsq_conn.execute("SELECT * FROM genome_metadata").fetchall()
				contig_md = self.gsq_conn.execute("SELECT * FROM contig_metadata").fetchall()
				gene_md = self.gsq_conn.execute("SELECT * FROM gene_metadata").fetchall()
				seqs = self.gsq_conn.execute("SELECT * FROM genome_sequences").fetchall()
				dest = genes_db(usepath, in_memory = False)
				dest.open()
				insert_status = dest.insert(genome_md, seqs, contig_md, gene_md, preds)
				dest.gsq_conn.commit()
				dest.gsq_conn.execute("VACUUM")
				dest.gsq_conn.commit()
				dest.close()
			else:
				print("Writing new database...")
				self.gsq_conn.commit()				
				dest = genomicsqlite.connect(usepath)
				self.gsq_conn.backup(dest)
				dest.close()
		else:
			print("Cleaning up additions...")
			self.gsq_conn.commit()
			self.gsq_conn.execute("VACUUM")
			self.gsq_conn.commit()
			
	def close(self):
		if self.gsq_conn is not None:
			self.gsq_conn.close()
					
	def give_current_indices(self):
		return self.genome_index, self.contig_index, self.gene_index, self.current_genome_index, self.current_contig_index, self.current_gene_index
		
	def set_indices(self, genome_index, contig_index, gene_index, next_genome_id, next_contig_id, next_gene_id):
		self.genome_index = genome_index
		self.contig_index = contig_index
		self.gene_index = gene_index
		
		self.current_genome_index = next_genome_id
		self.current_contig_index = next_contig_id
		self.current_gene_index = next_gene_id
		
		self.reverse_genome_index = dict(zip(self.genome_index.values(), self.genome_index.keys()))
		self.reverse_contig_index = dict(zip(self.contig_index.values(), self.contig_index.keys()))
		self.reverse_gene_index = dict(zip(self.gene_index.values(), self.gene_index.keys()))
		
		
	#Prepare the databases' schema
	def initialize(self):
		#SQLite automatically indexes on primary keys, so we mostly avoid indices here
		
		#It would be more efficient to just primary key on genome ID, BUT -
		#queries to this table will almost certainly be of the form "give me this genome" by name, and making something a key
		#causes SQLite to automatically index it for more efficient lookups of this kind.
		genome_table = ["source_file TEXT",
						"genome_name TEXT",
						"genome_id INTEGER",
						"genome_length INTEGER",
						"num_contigs INTEGER",
						"n50 INTEGER",
						"has_predicted_genes BOOLEAN",
						"num_genes INTEGER",
						"translation_table INTEGER",
						"coding_bases INTEGER",
						"coding_density FLOAT",
						"PRIMARY KEY (genome_name, genome_id)"]
						
		genome_table = ', '.join(genome_table)
						
		#Sequences are put into a separate table because the sequence data is large and makes querying for information other than the sequence slow
		#We want the genome metadata table to be fast", so we keep it fast
		genome_sequence_table = ["genome_id INTEGER",
							"contig_id INTEGER PRIMARY KEY",
							"genome_sequence BLOB", #genomicsqlite twobit encoding uses blob
							#To ensure that the sequences can be compressed, we replace all non-call characters with an "A" when inserting
							#We replace these "A" chars with "N" when retrieving the seq
							"loci_of_non_calls_replaced_with_A BLOB", 
							"FOREIGN KEY(genome_id) REFERENCES genome_metadata(genome_id)"]

		genome_sequence_table = ', '.join(genome_sequence_table)
		

		
		#here we have contigs by genome in the case of a multifasta.
		contig_table = ["genome_id INTEGER", 
						"contig_name TEXT",
						"contig_description TEXT",
						"contig_id INTEGER PRIMARY KEY",
						"contig_length INTEGER", 
						"num_genes INTEGER", 
						"coding_bases INTEGER", 
						"coding_density INTEGER",
						"FOREIGN KEY(genome_id) REFERENCES genome_metadata(genome_id)"]
						
		contig_table = ', '.join(contig_table)
		
		#Note to self
		#Eventually it's probably worth updating this to be a generic "features" table genes as CDS feature type
		#and using genomicsqlite's gene range indices to allow for selection of any other types genomic features.
		#Probably just include start/stop index and feature_type column with "feature info" stored as either a JSON or as a semicolon-sep string 
		#Each feature type with corresponding "digest_[CDS/SNP/deletion/insertion" functions
		genes_table = ["genome_id INTEGER",
					"contig_id INTEGER",
					"gene_name TEXT",
					"gene_id INTEGER PRIMARY KEY",
					"strand BOOLEAN",
					#Min and max because start/end are ambiguous wrt. strand - some programs mean "start" as the beginning of the gene even if the start pos is higher than the end.
					#Prodigal means "start" as the minimum position in the genome along the primary strand", even if thats the end of the gene - which it is for every complement strand gene
					"gene_minimum_position INTEGER",
					"gene_maximum_position INTEGER",
					"gene_annotation TEXT",
					"FOREIGN KEY(genome_id) REFERENCES genome_metadata(genome_id)",
					"FOREIGN KEY(contig_id) REFERENCES contig_metadata(contig_id)"]
				
		genes_table = ', '.join(genes_table)
		
		genome_table_sql = "CREATE TABLE IF NOT EXISTS genome_metadata ({table})".format(table = genome_table)
		genome_seq_table_sql = "CREATE TABLE IF NOT EXISTS genome_sequences ({table})".format(table = genome_sequence_table)
		contig_table_sql = "CREATE TABLE IF NOT EXISTS contig_metadata ({table})".format(table = contig_table)
		genes_table_sql = "CREATE TABLE IF NOT EXISTS gene_metadata ({table})".format(table = genes_table)
				
		self.gsq_conn.execute(genome_table_sql)
		self.gsq_conn.execute(genome_seq_table_sql)
		self.gsq_conn.execute(contig_table_sql)
		self.gsq_conn.execute(genes_table_sql)
		self.gsq_conn.commit()
		
	def get_metadata(self):
		if self.gsq_conn is not None:
			self.genome_index = {}
			self.reverse_genome_index = {}
			self.contig_index = {}
			self.reverse_contig_index = {}
			self.gene_index = {}
			self.reverse_gene_index = {}
			self.gene_prediction_record = {}
			
			if self.in_mem:
				#If the database is in memory and an on-disk database already exists, briefly connect and load the on-disk metadata and use it to set the in-memory db metadata
				if os.path.exists(self.dbpath):
					#Load an already-existing database's metadata so that genomes can be added instead of overwriting without loading the db's sequences
					existing = genes_db(self.dbpath, on_disk = True)
					existing.open()
					existing.get_metadata()
					self.genome_index, self.contig_index, self.gene_index, self.current_genome_index, self.current_contig_index, self.current_gene_index = existing.give_current_indices()
					existing.close()
					existing = None			

			#Database is a disk connection	
			else:
				#Load the db's metadata
				#genome_meta = self.gsq_conn.execute("SELECT genome_name, genome_id, has_predicted_genes FROM genome_metadata").fetchall()
				genome_meta = self.gsq_conn.execute("SELECT source_file, genome_id, has_predicted_genes FROM genome_metadata").fetchall()
				contig_meta = self.gsq_conn.execute("SELECT contig_name, contig_id FROM contig_metadata").fetchall()
				gene_meta = self.gsq_conn.execute("SELECT gene_name, gene_id FROM gene_metadata").fetchall()
				
				for row in genome_meta:
					genome, genome_id, has_preds = row[0], row[1], row[2]
					self.genome_index[genome] = genome_id
					self.reverse_genome_index[genome_id] = genome
					self.gene_prediction_record[genome_id] = has_preds
					
				if len(self.genome_index) > 0:
					self.current_genome_index = max(self.genome_index.values())+1
					
				for row in contig_meta:
					contig, contig_id = row[0], row[1]
					self.contig_index[contig] = contig_id
					self.reverse_contig_index[contig_id] = contig
					
				if len(self.contig_index) > 0:
					self.current_contig_index = max(self.contig_index.values())+1
					
				for row in gene_meta:
					gene, gene_id = row[0], row[1]
					self.gene_index[gene] = gene_id
					self.reverse_gene_index[gene_id] = gene
					
				if len(self.gene_index) > 0:
					self.current_gene_index = max(self.gene_index.values())+1
				
	#Set/update indices and insert
	def insert(self, genome_insert, sequence_insert, contig_insert, gene_insert = None, drop_these_genomes = None):
		#If a genome is being added to a database for a second time, first delete its records everywhere
		if drop_these_genomes is not None:
			if genome_insert is not None:
				self.gsq_conn.executemany('DELETE FROM genome_metadata WHERE genome_id = ?', drop_these_genomes)
			if sequence_insert is not None:
				self.gsq_conn.executemany('DELETE FROM genome_sequences WHERE genome_id = ?', drop_these_genomes)
			if contig_insert is not None:
				self.gsq_conn.executemany('DELETE FROM contig_metadata WHERE genome_id = ?', drop_these_genomes)
			if gene_insert is not None:
				self.gsq_conn.executemany('DELETE FROM gene_metadata WHERE genome_id = ?', drop_these_genomes)
			
		if genome_insert is not None:
			self.gsq_conn.executemany('INSERT INTO genome_metadata VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', genome_insert)
		
		if sequence_insert is not None:
			self.gsq_conn.executemany('INSERT INTO genome_sequences VALUES (?, ?, ?, ?)', sequence_insert)
		
		if contig_insert is not None:
			self.gsq_conn.executemany('INSERT INTO contig_metadata VALUES (?, ?, ?, ?, ?, ?, ?, ?)', contig_insert)
		
		if gene_insert is not None:
			self.gsq_conn.executemany('INSERT INTO gene_metadata VALUES (?, ?, ?, ?, ?, ?, ?, ?)', gene_insert)
		
		
		
		