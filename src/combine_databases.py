import genomicsqlite
from build_database import genes_db

class db_combiner:
	def __init__(self, recipient = None, donor = None):
		self.rec_path = recipient
		self.don_path = donor
		
		self.donor = genes_db(self.don_path)
		self.recipient = genes_db(self.rec_path)
		
		self.donor_to_rec_genomes = None
		self.donor_to_rec_contigs = None
		self.donor_to_rec_genes = None
		
		self.genomes_to_add = None
		
	def open(self):
		self.donor.open()
		self.recipient.open()
		
	def close(self):
		self.donor.close()
		self.recipient.close()
		
	def reindex(self):
		self.recipient.get_metadata()
		self.donor.get_metadata()
		
		self.donor_to_rec_genomes = {}
		self.donor_to_rec_contigs = {}
		self.donor_to_rec_genes = {}
		
		for donor_genome_id in sorted(self.donor.reverse_genome_index):
			genome_name = self.donor.reverse_genome_index[donor_genome_id]
			if genome_name not in self.recipient.genome_index:
				self.donor_to_rec_genomes[donor_genome_id] = self.recipient.current_genome_index
				self.recipient.current_genome_index += 1
			
		for donor_contig_id in sorted(self.donor.reverse_contig_index):
			contig_name = self.donor.reverse_contig_index[donor_contig_id]
			if contig_name not in self.recipient.contig_index:
				self.donor_to_rec_contigs[donor_contig_id] = self.recipient.current_contig_index
				self.recipient.current_contig_index += 1

			
		for donor_gene_id in sorted(self.donor.reverse_gene_index):
			gene_name = self.donor.reverse_gene_index[donor_gene_id]
			if gene_name not in self.recipient.gene_index:
				self.donor_to_rec_genes[donor_gene_id] = self.recipient.current_gene_index
				self.recipient.current_gene_index += 1
			
		self.genomes_to_add = list(self.donor_to_rec_genomes.keys())
		self.genomes_to_add.sort()
		
	def split_to_groups(self):
		sql_isin_limit = 998
		if len(self.genomes_to_add) > sql_isin_limit:
			grps = []
			for i in range(0, (len(self.genomes_to_add) // sql_isin_limit)+1):
				next_group = self.genomes_to_add[(i*sql_isin_limit):((i+1)*sql_isin_limit)]
				if len(next_group) > 0:
					grps.append(next_group)
			self.genomes_to_add = grps
				
		else:
			self.genomes_to_add = [self.genomes_to_add]
	
	def add_donor_to_recipient(self):
		print(len(self.genomes_to_add), "new genomes will be added to", self.rec_path)
		#Go one by one for sequences; make sure memory is never an issue
		for genome_id in self.genomes_to_add:
			new_genome_id = self.donor_to_rec_genomes[genome_id]
			next_seqs = self.donor.gsq_conn.execute("SELECT * FROM genome_sequences WHERE genome_id = ?", (genome_id,)).fetchall()
			reformatted = []
			for seq in next_seqs:
				seq = list(seq)
				seq[0] = new_genome_id
				new_contig_id = self.donor_to_rec_contigs[seq[1]]
				seq[1] = new_contig_id
				reformatted.append(tuple(seq))
			
			self.recipient.gsq_conn.executemany("INSERT INTO genome_sequences VALUES (?, ?, ?, ?)", reformatted)
			reformatted = None
		
		#Go group by group for metadata - chunkwise makes more sense.
		self.split_to_groups()
		
		reformatted_genomes = []
		reformatted_contigs = []
		reformatted_genes = []
		
		#Select data chunkwise
		for genome_group in self.genomes_to_add:
			binding_format = ', '.join(['?']*len(genome_group))
			sql1 = "SELECT * FROM genome_metadata WHERE genome_id IN ({genome_ids})".format(genome_ids = binding_format)
			sql2 = "SELECT * FROM contig_metadata WHERE genome_id IN ({genome_ids})".format(genome_ids = binding_format)
			sql3 = "SELECT * FROM gene_metadata WHERE genome_id IN ({genome_ids})".format(genome_ids = binding_format)
			donor_genomes = self.donor.gsq_conn.execute(sql1, genome_group).fetchall()

			reformatted_genomes = []
			reformatted_contigs = []
			reformatted_genes = []
			
			for row in donor_genomes:
				row = list(row)
				row[2] = self.donor_to_rec_genomes[row[2]]
				row = tuple(row)
				reformatted_genomes.append(row)
				
			donor_genomes = None
			
			donor_contigs = self.donor.gsq_conn.execute(sql2, genome_group).fetchall()
			
			for row in donor_contigs:
				row = list(row)
				row[0] = self.donor_to_rec_genomes[row[0]]
				row[3] = self.donor_to_rec_contigs[row[3]]
				row = tuple(row)
				reformatted_contigs.append(row)
			
			donor_contigs = None
			
			donor_genes = self.donor.gsq_conn.execute(sql3, genome_group).fetchall()
			for row in donor_genes:
				row = list(row)
				row[0] = self.donor_to_rec_genomes[row[0]]
				row[1] = self.donor_to_rec_contigs[row[1]]
				row[3] = self.donor_to_rec_genes[row[3]]
				
				row = tuple(row)
				reformatted_genes.append(row)
				
		
		#Check for no-inserts
		if len(reformatted_genomes) == 0:
			reformatted_genomes = None
		if len(reformatted_contigs) == 0:
			reformatted_contigs = None
		if len(reformatted_genes) == 0:
			reformatted_genes = None
				
		#Use the database insert method
		self.recipient.insert(reformatted_genomes, None, reformatted_contigs, reformatted_genes)
		self.recipient.write_changes()
		
	def run_combine(self):
		self.open()
		self.reindex()
		self.add_donor_to_recipient()
		self.close()
		
		
