import os
#Class for preparing genome, gene, contig, seq data into the respective lists of tuples expected by the database.
class data_formatter:
	def __init__(self):
		self.genome_index = None
		self.contig_index = None
		self.gene_index = None
		
		self.next_genome_id = None
		self.next_contig_id = None
		self.next_gene_id = None
		
		self.genomes_to_drop = None
				
	def set_indices(self, genome_index, contig_index, gene_index, next_genome_id, next_contig_id, next_gene_id):
		self.genome_index = genome_index
		self.contig_index = contig_index
		self.gene_index = gene_index
		
		self.next_genome_id = next_genome_id
		self.next_contig_id = next_contig_id
		self.next_gene_id = next_gene_id

	def get_n50(self, genome_length, contig_lengths):
		n50_arr = list(contig_lengths.values())
			
		#Bitwise and 1 is true if odd
		if genome_length & 1:
			n50_length = (genome_length + 1) / 2
		else:
			n50_length = genome_length / 2
			
		n50_arr.sort()
		list_idx = 0
		#No need to check if there's only one contig
		if len(n50_arr) > 1:
			while n50_length > 0:
				n50_length -= n50_arr[list_idx]
				list_idx += 1
		
		#Max list idx just means the n50 is the longest contig.
		if list_idx == len(n50_arr):
			list_idx -= 1
		
		n50_value = n50_arr[list_idx]
			
		return n50_value
	
	#Each of the following functions formats data into lists of tuples for insertion into a database

	def update_indices(self, genome_name, sequence_order, genes_dict):
		self.genomes_to_drop = []
		
		if genome_name not in self.genome_index:
			self.genome_index[genome_name] = self.next_genome_id
			self.next_genome_id += 1
		else:
			self.genomes_to_drop.append((self.genome_index[genome_name], ))
		
		if len(self.genomes_to_drop) == 0:
			self.genomes_to_drop = None
		
		for contig in sequence_order:
			if contig not in self.contig_index:
				self.contig_index[contig] = self.next_contig_id
				self.next_contig_id += 1
			
			if genes_dict is not None:
				for gene in genes_dict[contig]['genes']:
					if gene not in self.gene_index:
						self.gene_index[gene] = self.next_gene_id
						self.next_gene_id += 1
	'''
	genome_table = [	"source_file TEXT",
						"genome_name TEXT",
						"genome_description TEXT",
						"genome_id INTEGER",
						"genome_length INTEGER",
						"num_contigs INTEGER",
						"n50 INTEGER",
						"num_genes INTEGER",
						"translation_table INTEGER",
						"coding_bases INTEGER",
						"coding_density FLOAT",
						"PRIMARY KEY (genome_name, genome_id)"]
	'''
	
	def create_genome_metadata_insert(self, genome_name, genome_length, num_contigs, n50, has_genes, num_genes, tt, coding_bases, coding_dens):
		insertion = []
		genome_id = self.genome_index[genome_name]
		genome_name_truncated = os.path.basename(genome_name)
		
		
		next_row = (genome_name, genome_name_truncated, genome_id, genome_length, num_contigs, n50, has_genes, num_genes, tt, coding_bases, coding_dens,)
		insertion.append(next_row)
		
		if len(insertion) == 0:
			insertion = None
		
		return insertion
		

	'''
	genome_sequence_table = ["genome_id INTEGER"
							"contig_id INTEGER PRIMARY KEY",
							"genome_sequence BLOB", #genomicsqlite twobit encoding uses blob
							#To ensure that the sequences can be compressed, we replace all non-call characters with an "A" when inserting
							#We replace these "A" chars with "N" when retrieving the seq
							"loci_of_non_calls_replaced_with_A BLOB", 
							"FOREIGN KEY(genome_id) REFERENCES genome_metadata(genome_id)"]
	'''
	def create_genome_seq_insert(self, genome_name, seq_order, encoded_seqs, loci_bytestrings):
		genome_id = self.genome_index[genome_name]
		insertion = []
		for contig in seq_order:
			next_row = (genome_id, self.contig_index[contig], encoded_seqs[contig], loci_bytestrings[contig],)
			insertion.append(next_row)
			
		if len(insertion) == 0:
			insertion = None
				
		return insertion
	
	#Format data into tuples of:
	'''
	contig_table = ["genome_id INTEGER", 
					"contig_name TEXT", 
					"contig_id INTEGER PRIMARY KEY",
					"contig_length INTEGER", 
					"num_genes INTEGER", 
					"coding_bases INTEGER", 
					"coding_density INTEGER",
					"FOREIGN KEY(genome_id) REFERENCES genome_metadata(genome_id)"]
	'''
	def create_contig_metadata_insert(self, genome_name, seq_order, seq_descriptions, contig_lengths, genes_dict):
		genome_id = self.genome_index[genome_name]
		insertion = []
		for contig in seq_order:
			desc = seq_descriptions[contig]
			contig_id = self.contig_index[contig]
			contig_length = contig_lengths[contig]
			
			if genes_dict is not None:
				num_genes = genes_dict[contig]['num_genes']
				coding_bases = genes_dict[contig]['coding_bases']
				coding_density = genes_dict[contig]['coding_density']
			else:
				num_genes = None
				coding_bases = None
				coding_density = None
				
			next_row = (genome_id, contig, desc, contig_id, contig_length, num_genes, coding_bases, coding_density,)
			insertion.append(next_row)
			
		if len(insertion) == 0:
			insertion = None
			
		return insertion
		
	#Format data into tuples of:
	'''
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
	'''
	def create_gene_metadata_insert(self, genome_name, seq_order, gene_dict):
		insertion = []
		if gene_dict is not None:
			genome_id = self.genome_index[genome_name]
			for contig in seq_order:
				contig_id = self.contig_index[contig]
				for gene in gene_dict[contig]['genes']:
					gene_id = self.gene_index[gene]
					#(gene_strand, gene_start, gene_end, sequence_description)
					strand = gene_dict[contig]['genes'][gene][0]
					min_position = gene_dict[contig]['genes'][gene][1]
					max_position = gene_dict[contig]['genes'][gene][2]
					gene_annot = gene_dict[contig]['genes'][gene][3]
					
					next_row = (genome_id, contig_id, gene, gene_id, strand, min_position, max_position, gene_annot,)
					insertion.append(next_row)
		
		if len(insertion) == 0:
			insertion = None
		
		return insertion
		
	def add_next_data(self, genome_name, genome_length, translation_table, whole_genome_coding_bases, whole_genome_coding_dens, \
					sequence_order, sequences, contig_lengths, descriptions, noncalls_replaced_with_A_indices, genes_by_genome):
		
		#Add new genomes, contigs, genes to indices
		self.update_indices(genome_name, sequence_order, genes_by_genome)
		
		n50 = self.get_n50(genome_length, contig_lengths)
		
		num_contigs = len(sequence_order)

		if genes_by_genome is not None:
			has_genes = True
			num_genes = 0
			for contig in genes_by_genome:
				num_genes += len(genes_by_genome[contig]['genes'])
		else:
			has_genes = False
			num_genes = None
				
		
		genome_md_ins = self.create_genome_metadata_insert(genome_name, genome_length, num_contigs, n50, has_genes, num_genes, translation_table, whole_genome_coding_bases, whole_genome_coding_dens)
		
		if sequences is not None:
			genome_seq_ins = self.create_genome_seq_insert(genome_name, sequence_order, sequences, noncalls_replaced_with_A_indices)
		else:
			genome_seq_ins = None
		
		contig_md_ins = self.create_contig_metadata_insert(genome_name, sequence_order, descriptions, contig_lengths, genes_by_genome, )
		
		if has_genes:
			gene_md_ins = self.create_gene_metadata_insert(genome_name, sequence_order, genes_by_genome)
		else:
			gene_md_ins = None
		
		
		return genome_md_ins, genome_seq_ins, contig_md_ins, gene_md_ins, self.genomes_to_drop
		
	def give_current_indices(self):
		return self.genome_index, self.contig_index, self.gene_index, self.next_genome_id, self.next_contig_id, self.next_gene_id
	
	