import pyrodigal
import io
import numpy as np

class pyrodigal_manager:
	def __init__(self, trans_tables = [11, 4], meta = False, verbose = False):
		self.genome_dict = None
		self.training_seq = None
		
		self.genome_length = 0
		self.contig_lengths = None
		
		self.meta = meta
		self.verbose = verbose
		
		self.training_seq = None
		self.gene_finder = None
		
		self.coding_bases_by_table = {}
		self.coding_density_by_table = {}
		
		self.trans_tables = trans_tables
		
		self.current_genes = None
		self.winning_coding_dens = 0
		self.winning_table = None
		
		self.protein_seqs = {}
		self.protien_length_nt = 0
		
	#Reset to initial conditions so that the current use of this object doesn't interfere with the next
	def reset(self):
		self.genome_dict = None
		self.training_seq = None
		
		self.genome_length = 0
		self.contig_lengths = {}
		
		self.training_seq = None
		self.gene_finder = None
		
		self.coding_bases_by_table = {}
		self.coding_density_by_table = {}
				
		self.current_genes = None
		self.winning_coding_dens = 0
		self.winning_table = None
		
		self.protein_seqs = {}
		self.protein_length_nt = 0
		
	def prep_genome_dict_for_prediction(self, genome_dict):
		self.reset()
		self.genome_dict = {}
		for g in genome_dict:
			#asbytes = g.encode(encoding = "ascii")
			contig_length = len(genome_dict[g])
			self.contig_lengths[g] = contig_length
			self.genome_dict[g] = genome_dict[g].encode(encoding = "ascii")
		
		self.genome_length = sum(self.contig_lengths.values())
		
	def prep_training_seq(self):
		running_sum = 0
		seqs_added = 0
		if self.training_seq is None:
			self.training_seq = []
			for seq in self.genome_dict:
				running_sum += len(self.genome_dict[seq])
				if seqs_added > 0:
					#Prodigal interleaving logic - add this breaker between sequences, starting at sequence 2
					self.training_seq.append(b'TTAATTAATTAA')
					running_sum += 12
					
				seqs_added += 1
					
				#Handle excessive size
				if running_sum >= 32000000:					
					print("Warning:  Sequence is long (max 32000000 for training).")
					print("Training on the first 32000000 bases.")
				
					to_remove = running_sum - 32000000
					
					#Remove excess characters
					cut_seq = self.genome_dict[seq][:-to_remove]
					#Add the partial seq
					self.training_seq.append(cut_seq)
					
					#Stop the loop and move to training
					break
				
				#add in a full sequence
				self.training_seq.append(self.genome_dict[seq])

			if seqs_added > 1:
				self.training_seq.append(b'TTAATTAATTAA')
				
			self.training_seq = b''.join(self.training_seq)
		
		if len(self.training_seq) < 20000:
			if self.verbose:
				print("Can't train on 20 thousand or fewer characters. Switching to meta mode.")
			self.manager = pd.GeneFinder(meta=True)
			self.meta = Trues
		else:
			if self.verbose:
				#G is 71, C is 67; we're counting G + C and dividing by the total.
				gc = round(((self.training_seq.count(67) + self.training_seq.count(71))/ len(self.training_seq)) * 100, 2)
				print(len(self.training_seq), "bp training seq created,", gc, "pct GC")
				
	def train(self, table = 11):
		if not self.meta:
			self.gene_finder = pyrodigal.GeneFinder(meta = False)
			self.gene_finder.train(self.training_seq, translation_table = table)
		else:
			self.gene_finder = pyrodigal.GeneFinder(meta = True)
			
	def predict(self):
		current_table = self.gene_finder.training_info.translation_table
		
		if current_table is None:
			current_table = "meta"
			
		self.coding_bases_by_table[current_table] = 0
		
		next_gene_set = {}
		
		if self.gene_finder is not None:
			for seq in self.genome_dict:
				genes = self.gene_finder.find_genes(self.genome_dict[seq])
				next_gene_set[seq] = genes
				for g in genes:	
					self.coding_bases_by_table[current_table] += len(g.sequence())
		
		self.coding_density_by_table[current_table] = self.coding_bases_by_table[current_table] / self.genome_length
		
		if self.current_genes is None:
			self.current_genes = next_gene_set
			self.winning_coding_dens = self.coding_density_by_table[current_table]
			self.winning_table = current_table
		else:
			if self.coding_density_by_table[current_table] > (self.winning_coding_dens * 1.1):
				if self.verbose:
					print("Translation table", self.winning_table, "had coding density of", self.winning_coding_dens)
					print("New translation table", current_table, "more than 10% better, coding density of", self.coding_density_by_table[current_table])
					print("New table's predicted genes will be used")
				self.current_genes = next_gene_set
				self.winning_coding_dens = self.coding_density_by_table[current_table]
				self.winning_table = current_table
			
		
	def extract_data(self):
		genes_text = io.StringIO("")
		
		for seq in self.current_genes:
			self.current_genes[seq].write_gff(genes_text, seq)
			
		genes_text = genes_text.getvalue()
		
		digested_gff = self.digest_gff(genes_text)
				
		return digested_gff
		
	def digest_gff(self, gff_text):
		genes_by_genome = {}
		seqdata = None
		moddata = None

		seqid = None
		seqlen = None
		
		lines = gff_text.splitlines()
		for line in lines:
			if line.startswith("#"):
				'''
				GFF3 header looks like:
				##gff-version  3
				# Sequence Data: seqnum=27;seqlen=56291;seqhdr="2012990006.a:JC3AJCVIAssemblies_1106445189678"
				# Model Data: version=pyrodigal.v3.5.2;run_type=Single;model="Ab initio";gc_cont=47.44;transl_table=11;uses_sd=0
				'''
				#GFF version header. Skip it.
				if line.startswith("##"):
					pass
				#Sequence data line: we want the seqlen and the seqhdr values
				if line.startswith("# Sequence Data: "):
					seqdata = line.split("Sequence Data: ")[1]
					seqdata = seqdata.split(";")
					seqlen = seqdata[1]
					seqlen = int(seqlen[7:])
					
					seqid = seqdata[2]
					seqid = seqid[8:-1]
					
					if seqid not in genes_by_genome:
						genes_by_genome[seqid] = {}
						genes_by_genome[seqid]["model_data"] = None
						genes_by_genome[seqid]["genes"] = {}
				
				#We're just gonna leave this as a string in the database
				if line.startswith("# Model Data: "):
					moddata = line.split('Model Data: ')[1]
					genes_by_genome[seqid]["model_data"] = moddata
					
			else:
				#GFF3 content looks like:
				#Assemblies_1106445189264_66;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.567;conf=100.00;score=162.40;cscore=152.84;sscore=9.56;rscore=0.05;uscore=4.36;tscore=3.89;
				
				segs = line.strip().split()
				gene_strand = segs[6]
				
				#Convert GFF3 format +/- to True/False for storage in sqlitedb
				if gene_strand == "+":
					gene_strand = True
				else:
					gene_strand = False
					
				gene_start = int(segs[3])
				gene_end = int(segs[4])

				#These are generated with pyrodigal, so they will always look like this:
				#ID=2012990006.a:JC3AJCVIAssemblies_1106445189678_57;partial=00;start_type=GTG;rbs_motif=None;rbs_spacer=None;gc_cont=0.619;conf=100.00;score=104.60;cscore=105.19;sscore=-0.59;rscore=0.05;uscore=1.61;tscore=-1.60;
				attributes = segs[8]
				attributes = attributes[3:-1]
				attributes = attributes.split(';')
				
				gene_id = attributes[0]
				#Format as prodigal gene header
				sequence_description = ';'.join(attributes)
				
				genes_by_genome[seqid]["genes"][gene_id] = (gene_strand, gene_start, gene_end, sequence_description)
							
				
		return genes_by_genome
	
	def summarize_genes_by_genome(self, genes_by_genome):
		genome_coding_bases = 0
	
		for seq in genes_by_genome:
			num_genes = len(genes_by_genome[seq]['genes'])
			genes_by_genome[seq]['num_genes'] = num_genes

			genes_by_genome[seq]['contig_length'] = self.contig_lengths[seq]
			
			#create an array of essentially bools (numpy uses an int8 as bool under the hood) to represent the number of unique coding bases in this contig
			coding_bases = np.zeros(self.contig_lengths[seq], dtype = np.int8)
			
			for g in genes_by_genome[seq]['genes']:
				gene_start = genes_by_genome[seq]['genes'][g][1]
				gene_end = genes_by_genome[seq]['genes'][g][2]
				#prodigal gene starts are 1-indexed, python is 0-indexed; start -1 to correct, end is left as-is because gene ends are inclusive
				coding_bases[(gene_start-1):(gene_end)] = 1
			
			genes_by_genome[seq]['coding_bases'] = int(np.sum(coding_bases))
			genome_coding_bases += genes_by_genome[seq]['coding_bases']
			genes_by_genome[seq]['coding_density'] = float(genes_by_genome[seq]['coding_bases'] / genes_by_genome[seq]['contig_length'])
			
		genome_coding_bases = int(genome_coding_bases)
			
		return genes_by_genome, genome_coding_bases
			
		
	def predict_genes(self, genome_dict):
		self.prep_genome_dict_for_prediction(genome_dict)
		self.prep_training_seq()
		for table in self.trans_tables:
			self.train(table = table)
			self.predict()
		
		genes_by_genome = self.extract_data()
		genes_by_genome, whole_genome_coding_bases = self.summarize_genes_by_genome(genes_by_genome)
		
		whole_genome_coding_dens = float(whole_genome_coding_bases / self.genome_length)
		
		'''
		Below are the tables which are filled with data from this part
		
		#From this script: num_genes, translation table, coding_bases, coding_density
		genome_table = ["genome_name TEXT",
						"genome_description TEXT",
						"genome_id INTEGER", #This will be a placeholder with value -1 in this script
						"genome_length INTEGER",
						"num_contigs INTEGER",
						"n50 INTEGER",
						"num_genes INTEGER",
						"translation_table INTEGER",
						"coding_bases INTEGER",
						"coding_density FLOAT",
						"PRIMARY KEY (genome_name, genome_id)"]
		
		#From this script:
		#Each seq is a contig
		genes_by_genome[seq]['num_genes']
		genes_by_genome[seq]['contig_length'] = self.contig_lengths[seq]
		genes_by_genome[seq]['coding_bases'] = np.sum(coding_bases)
		genome_coding_bases += genes_by_genome[seq]['coding_bases']
		genes_by_genome[seq]['coding_density'] = genes_by_genome[seq]['coding_bases'] / genes_by_genome[seq]['contig_length']
		
		contig_table = ["genome_id INTEGER", 
						"contig_name TEXT", 
						"contig_id INTEGER PRIMARY KEY", #This will be a placeholder with value -1 in this script
						"contig_length INTEGER", 
						"num_genes INTEGER", 
						"coding_bases INTEGER", 
						"coding_density INTEGER",
						"FOREIGN KEY(genome_id) REFERENCES genome_metadata(genome_id)"]
						
		genes_table = ["contig_id INTEGER", #This will be a placeholder with contig_name in this script
					"gene_name TEXT",
					"gene_id INTEGER PRIMARY KEY", #This will be a placeholder with value -1 in this script
					"strand BOOLEAN",
					#Min and max because start/end are ambiguous wrt. strand - some programs mean "start" as the beginning of the gene even if the start pos is higher than the end.
					#Prodigal means "start" as the minimum position in the genome along the primary strand", even if thats the end of the gene - which it is for every complement strand gene
					"gene_minimum_position INTEGER",
					"gene_maximum_position INTEGER",
					"gene_annotation TEXT",
					"FOREIGN KEY(contig_id) REFERENCES contig_metadata(contig_id)"]
		'''
		
		return genes_by_genome, self.winning_table, whole_genome_coding_bases, whole_genome_coding_dens
