import genomicsqlite
import pyhmmer
import numpy as np

import gzip
import json

import os
import io

class db_extractor:
	def __init__(self, dbpath = None, output_dir = "baked_goods", compress = False, in_mem = False):
		self.gsq_conn = None
		
		self.db_meta = None
		
		self.compress = compress
		
		self.dbpath = dbpath
		self.outdir = output_dir
		self.in_memory = in_mem
		
		self.genome_index = None
		self.reverse_genome_index = None
		self.contig_index = None
		self.reverse_contig_index = None
		self.gene_index = None
		self.reverse_gene_index = None
		
		self.current_genome_index = 0
		self.current_contig_index = 0
		self.current_gene_index = 0
		
		self.metadata_prefix = None
		self.genomes_prefix = None
		self.protein_nt_prefix = None
		self.protein_aa_prefix = None
		
		self.dna_alpha = pyhmmer.easel.Alphabet.dna()
		
		self.rc = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
		
	def open(self):
		if self.gsq_conn is None:
			if os.path.exists(self.dbpath):
				if self.in_memory:
					self.gsq_conn = genomicsqlite.connect(":memory:")
					if os.path.exists(self.dbpath):
						#Load an already-existing database into memory
						existing = genomicsqlite.connect(self.dbpath)
						existing.backup(self.gsq_conn)
						existing.close()				
				else:
					self.gsq_conn = genomicsqlite.connect(self.dbpath, read_only = True)
					
				self.get_metadata()
			else:
				print("No database found!")	
		
	def close(self):
		if self.gsq_conn is not None:
			self.gsq_conn.close()
	
	def ready_converter(self):
		self.pyhmmer_converter = None
	
	def get_metadata(self):
		if self.gsq_conn is not None:
			self.genome_index = {}
			self.reverse_genome_index = {}
			self.contig_index = {}
			self.reverse_contig_index = {}
			self.gene_index = {}
			self.reverse_gene_index = {}
		
			if self.in_memory:
				#If the database is in memory and an on-disk database already exists, briefly connect and load the on-disk metadata and use it to set the in-memory db metadata
				if os.path.exists(self.dbpath):
					#Load an already-existing database's metadata so that genomes can be added instead of overwriting without loading the db's sequences
					existing = genes_db(self.dbpath, on_disk = True)
					existing.open()
					existing.get_metadata()
					self.genome_index, self.contig_index, self.gene_index, self.current_genome_index, self.current_contig_index, self.current_gene_index = existing.give_current_indices()
					existing.close()
					existing = None
			else:
				#Load the db's metadata
				genome_meta = self.gsq_conn.execute("SELECT genome_name, genome_id FROM genome_metadata").fetchall()
				contig_meta = self.gsq_conn.execute("SELECT contig_name, contig_id FROM contig_metadata").fetchall()
				gene_meta = self.gsq_conn.execute("SELECT gene_name, gene_id FROM gene_metadata").fetchall()
				
				for row in genome_meta:
					genome, genome_id = row[0], row[1]
					self.genome_index[genome] = genome_id
					self.reverse_genome_index[genome_id] = genome
					
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
		
	def extract_database_description(self):
		self.db_meta = {
						'origin_database':self.dbpath,
						'num_genomes':len(self.genome_index),
						'genome_index':self.genome_index,
						'genomes':{}
						}
		if self.gsq_conn is not None:
			genomes = self.gsq_conn.execute("SELECT * FROM genome_metadata").fetchall()
			contigs = self.gsq_conn.execute("SELECT * FROM contig_metadata").fetchall()
			genes = self.gsq_conn.execute("SELECT * FROM gene_metadata").fetchall()
			#test_data/genomes_2/2013843002_3.fna.gz 2013843002_3.fna.gz 0 1385230 34 107199 1462 11 1246329 0.8997271211278993
			for row in genomes:
				src = row[0]
				genome_name = row[1]
				genome_id = row[2]
				genome_length = row[3]
				num_contigs = row[4]
				n50 = row[5]
				num_genes_in_genome = row[6]
				translation_table = row[7]
				coding_bases_in_genome = row[8]
				coding_base_fraction = row[9]
				
				
				self.db_meta['genomes'][genome_name] = {
											'source_file':src,
											'genome_id_in_database': genome_id,
											'genome_length': genome_length,
											'number_of_contigs': num_contigs,
											'n50': n50,
											'number_of_genes_in_genome': num_genes_in_genome,
											'translation_table_for_gene_prediction': translation_table,
											'coding_base_count_in_genome': coding_bases_in_genome,
											'fraction_of_coding_bases_in_genome': coding_base_fraction,
											'contigs': {}
											}
			
			#10 2029527003.a:APTF_contig21503 [description] 1401 42755 43 37933 0.8872178692550579
			for row in contigs:
				genome_id = row[0]
				genome_name = self.reverse_genome_index[genome_id]
				contig_name = row[1]
				contig_description = row[2]
				contig_id = row[3]
				contig_length = row[4]
				num_genes_in_contig = row[5]
				coding_bases_in_contig = row[6]
				coding_fraction_in_contig = row[7]
				

				if len(contig_description) == 0:
					contig_description = "SeqID_only"
				
				self.db_meta['genomes'][genome_name]['contigs'][contig_name] = {
																			'contig_id_in_database':contig_id,
																			'contig_description':contig_description,
																			'contig_length':contig_length,
																			'number_of_genes_in_contig':num_genes_in_contig, 
																			'coding_base_count_in_contig':coding_bases_in_contig, 
																			'fraction_of_coding_bases_in_contig':coding_fraction_in_contig, 
																			'genes':{} 
																			}
			
			#10 1401 2029527003.a:APTF_contig21503_43 26051 0 40330 42186 2029527003.a:APTF_contig21503_43 # partial=00 # start_type=ATG # rbs_motif=TAA # rbs_spacer=8bp # gc_cont=0.415 # conf=99.99 # score=215.30 # cscore=209.45 # sscore=5.84 # rscore=2.04 # uscore=-0.20 # tscore=4.00													
			for row in genes:
				genome_id = row[0]
				contig_id = row[1]
				genome_name = self.reverse_genome_index[genome_id]
				contig_name = self.reverse_contig_index[contig_id]
				
				gene_name = row[2]
				gene_id = row[3]
				gene_strand = row[4]
				if gene_strand:
					gene_strand = "+"
					gene_start_position = row[5]
					gene_end_position = row[6]
				else:
					gene_strand = "-"
					gene_start_position = row[6]
					gene_end_position = row[5]
					
				gene_annotation = row[7]
				
				self.db_meta['genomes'][genome_name]['contigs'][contig_name]['genes'][gene_name] = {
																								'gene_id_in_database':gene_id,
																								'strand':gene_strand,
																								'one_indexed_gene_start':gene_start_position,
																								'one_indexed_gene_end':gene_end_position,
																								'annotation':gene_annotation
																								}
				
	
	def check_or_make_output_dir(self, path):
		if not os.path.exists(path):
			os.makedirs(path, exist_ok = True)
	
	#Setup output
	def prep_files(self):
		self.check_or_make_output_dir(self.outdir)
		
		self.metadata_prefix = os.path.normpath(self.outdir + "/metadata")
		self.genomes_prefix = os.path.normpath(self.outdir + "/genome_fasta")
		self.protein_nt_prefix = os.path.normpath(self.outdir + "/gene_nt_fasta")
		self.protein_aa_prefix = os.path.normpath(self.outdir + "/gene_aa_fasta")
		
		self.check_or_make_output_dir(self.metadata_prefix)
		self.check_or_make_output_dir(self.genomes_prefix)
		self.check_or_make_output_dir(self.protein_nt_prefix)
		self.check_or_make_output_dir(self.protein_aa_prefix)
		
	def format_bp_per_line(self, genome_string, num_char_per_line = 70):
		#ceiling funciton without the math module
		ceiling = int(round((len(genome_string)/num_char_per_line)+0.5, 0))
		formatted = '\n'.join([genome_string[(i*num_char_per_line):(i+1)*num_char_per_line] for i in range(0, ceiling)])
		return formatted
	
	#Function for handling the zipping of text if needed and the opening/writing to a file.
	def write_to_file(self, file, text):
		if self.compress:
			text = text.encode(encoding = "ascii")
			text = gzip.compress(text)
			with open(file+".gz", "wb") as outwriter:
				outwriter.write(text)
		else:
			with open(file, "w") as outwriter:
				outwriter.write(text)
	
	#Write all of the metadata in a friendly json
	def write_json(self):
		#Will have to check what appropriate indent is, it's a pretty deep dict
		json_data = json.dumps(self.db_meta, indent = 4)
		self.write_to_file(os.path.normpath(self.metadata_prefix+"/metadata_JSON.txt"), json_data)
	
	#Three metadata files for genomes, contigs, genes. Written with human readability in mind.
	def write_text_triplet(self):
		genome_header = ['source_file', 'genome_name', 'genome_length', 'num_contigs', 'n50', 'num_genes', 'trans_table', 'coding_bases', 'coding_fraction']
		contig_header = ['genome_name', 'contig_name', 'contig_description', 'contig_length', 'num_genes', 'coding_bases', 'coding_fraction']
		gene_header = ['genome_name', 'contig_name', 'gene_name', 'gene_strand', 'gene_start', 'gene_end', 'annotation']
		
		genome_output = [genome_header]
		contig_output = [contig_header]
		gene_output = [gene_header]
		
		for genome in self.db_meta['genomes']:
			#add genome data
			src = self.db_meta['genomes'][genome]['source_file']
			glen = self.db_meta['genomes'][genome]['genome_length']
			nctg = self.db_meta['genomes'][genome]['number_of_contigs']
			n50 = self.db_meta['genomes'][genome]['n50']
			ngenes = self.db_meta['genomes'][genome]['number_of_genes_in_genome']
			tt = self.db_meta['genomes'][genome]['translation_table_for_gene_prediction']
			cbg = self.db_meta['genomes'][genome]['coding_base_count_in_genome']
			cfg = self.db_meta['genomes'][genome]['fraction_of_coding_bases_in_genome']
			
			next_genome_info = [src, genome, str(glen), str(nctg), str(n50), str(ngenes), str(tt), str(cbg), str(round(cfg, 6))]
			genome_output.append(next_genome_info)
			
			for contig in self.db_meta['genomes'][genome]['contigs']:		
				desc = self.db_meta['genomes'][genome]['contigs'][contig]['contig_description']			
				ctl = self.db_meta['genomes'][genome]['contigs'][contig]['contig_length']
				ngenes_ct = self.db_meta['genomes'][genome]['contigs'][contig]['number_of_genes_in_contig']
				cbct = self.db_meta['genomes'][genome]['contigs'][contig]['coding_base_count_in_contig']
				cfct = self.db_meta['genomes'][genome]['contigs'][contig]['fraction_of_coding_bases_in_contig']
				
				next_contig_info = [genome, contig, desc, str(ctl), str(ngenes_ct), str(cbct), str(round(cfct, 6))]
				contig_output.append(next_contig_info)
				
				for gene in self.db_meta['genomes'][genome]['contigs'][contig]['genes']:				
					strand = self.db_meta['genomes'][genome]['contigs'][contig]['genes'][gene]['strand']
					st = self.db_meta['genomes'][genome]['contigs'][contig]['genes'][gene]['one_indexed_gene_start']
					nd = self.db_meta['genomes'][genome]['contigs'][contig]['genes'][gene]['one_indexed_gene_end']
					annot = self.db_meta['genomes'][genome]['contigs'][contig]['genes'][gene]['annotation']
			
					next_gene_info = [genome, contig, gene, strand, str(st), str(nd), annot]
					gene_output.append(next_gene_info)
			
		genome_output = ['\t'.join(l) for l in genome_output]
		contig_output = ['\t'.join(l) for l in contig_output]
		gene_output = ['\t'.join(l) for l in gene_output]
		
		genome_output = '\n'.join(genome_output) + "\n"
		contig_output = '\n'.join(contig_output) + "\n"
		gene_output = '\n'.join(gene_output) + "\n"
		
		self.write_to_file(os.path.normpath(self.metadata_prefix+"/genome_tsv.txt"), genome_output)
		self.write_to_file(os.path.normpath(self.metadata_prefix+"/contig_tsv.txt"), contig_output)
		self.write_to_file(os.path.normpath(self.metadata_prefix+"/gene_tsv.txt"), gene_output)

	#Highly redundant flat file containing all of the metadata, but with row-wise replication of genome and contig data per gene.
	def write_merged_text(self):
		merged_header = ['source_file', 'genome_name', 'genome_length', 'num_contigs', 'n50', 'num_genes_in_genome', 'trans_table', 'coding_bases_in_genome', \
		'coding_fraction_in_genome', 'contig_name', 'contig_description', 'contig_length', 'num_genes_in_contig', 'coding_bases_in_contig', 'coding_fraction_in_contig' \
		'gene_name', 'gene_strand', 'gene_start', 'gene_end', 'annotation']
		
		merged_output = [merged_header]

		for genome in self.db_meta['genomes']:
			#add genome data
			src = self.db_meta['genomes'][genome]['source_file']
			glen = self.db_meta['genomes'][genome]['genome_length']
			nctg = self.db_meta['genomes'][genome]['number_of_contigs']
			n50 = self.db_meta['genomes'][genome]['n50']
			ngenes = self.db_meta['genomes'][genome]['number_of_genes_in_genome']
			tt = self.db_meta['genomes'][genome]['translation_table_for_gene_prediction']
			cbg = self.db_meta['genomes'][genome]['coding_base_count_in_genome']
			cfg = self.db_meta['genomes'][genome]['fraction_of_coding_bases_in_genome']
						
			for contig in self.db_meta['genomes'][genome]['contigs']:		
				desc = self.db_meta['genomes'][genome]['contigs'][contig]['contig_description']			
				ctl = self.db_meta['genomes'][genome]['contigs'][contig]['contig_length']
				ngenes_ct = self.db_meta['genomes'][genome]['contigs'][contig]['number_of_genes_in_contig']
				cbct = self.db_meta['genomes'][genome]['contigs'][contig]['coding_base_count_in_contig']
				cfct = self.db_meta['genomes'][genome]['contigs'][contig]['fraction_of_coding_bases_in_contig']
				
				for gene in self.db_meta['genomes'][genome]['contigs'][contig]['genes']:				
					strand = self.db_meta['genomes'][genome]['contigs'][contig]['genes'][gene]['strand']
					st = self.db_meta['genomes'][genome]['contigs'][contig]['genes'][gene]['one_indexed_gene_start']
					nd = self.db_meta['genomes'][genome]['contigs'][contig]['genes'][gene]['one_indexed_gene_end']
					annot = self.db_meta['genomes'][genome]['contigs'][contig]['genes'][gene]['annotation']
			
					next_merged_info = [src, genome, str(glen), str(nctg), str(n50), str(ngenes), str(tt), str(cbg), \
					str(round(cfg, 6)), contig, desc, str(ctl), str(ngenes_ct), str(cbct), str(round(cfct, 6)),\
					gene, strand, str(st), str(nd), annot]
					
					merged_output.append(next_merged_info)
			
		merged_output = ['\t'.join(l) for l in merged_output]
		
		merged_output = '\n'.join(merged_output) + "\n"
		
		self.write_to_file(os.path.normpath(self.metadata_prefix+"/merged_metadata_tsv.txt"), merged_output)

	#Use pyhmmer to manage reverse complement, nt -> AA translation
	def gene_to_aa(self, seq, strand, trans_table):
		nt_seq = pyhmmer.easel.TextSequence(sequence = seq) #Store the sequence as a 
		nt_seq = nt_seq.digitize(alphabet = self.dna_alpha)
		
		if strand == "-":
			nt_seq.reverse_complement(inplace=True)
			
		gencode = pyhmmer.easel.GeneticCode(translation_table = trans_table)
		aa_seq = nt_seq.translate(genetic_code = gencode)
		
		nt_seq = nt_seq.textize()
		aa_seq = aa_seq.textize()
		
		return nt_seq.sequence, aa_seq.sequence

	
	def iterate_genomes(self, do_genomes = True, do_genes_nt = False, do_genes_aa = False):
		for genome in self.db_meta['genomes']:
			self.load_genome(genome, do_genomes, do_genes_nt, do_genes_aa)
			
	
	def load_genome(self, genome_name, do_genome = True, do_nt = False, do_aa = False):
		this_translation_table = self.db_meta['genomes'][genome_name]['translation_table_for_gene_prediction']
		
		genome_seqs = []
		nt_seqs = []
		aa_seqs = []
		
		genome_id = self.genome_index[genome_name]
		#We carefully ensured that contigs were inserted in the correct sequence order earlier, which also ensures they are loaded in the right order here.
		sql = "SELECT contig_id, twobit_dna(genome_sequence), loci_of_non_calls_replaced_with_A FROM genome_sequences WHERE genome_id = ? ORDER BY contig_id"
		for sequence_tuple in self.gsq_conn.execute(sql, (genome_id,)).fetchall():
			
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
			
			contig_id = sequence_tuple[0]
			contig_name = self.reverse_contig_index[contig_id]
			
			desc = self.db_meta['genomes'][genome_name]['contigs'][contig_name]['contig_description']
			
			if do_genome:
				formatted = self.format_bp_per_line(seq, 70)
				genome_seqs.append(">"+contig_name+" "+desc)
				genome_seqs.append(formatted)
			
			if do_nt or do_aa:
				for gene in self.db_meta['genomes'][genome_name]['contigs'][contig_name]['genes']:
					
					annot = self.db_meta['genomes'][genome_name]['contigs'][contig_name]['genes'][gene]['annotation']
					strand = self.db_meta['genomes'][genome_name]['contigs'][contig_name]['genes'][gene]['strand']
					
					gene_start = self.db_meta['genomes'][genome_name]['contigs'][contig_name]['genes'][gene]['one_indexed_gene_start']
					gene_end = self.db_meta['genomes'][genome_name]['contigs'][contig_name]['genes'][gene]['one_indexed_gene_end']
					
					gene_low, gene_high = min([gene_start, gene_end]), max(gene_start, gene_end)
					
					#prodigal header format
					
					if strand == "+":
						num_strand = "1"
					else:
						num_strand = "-1"
					
					gene_printout = " # ".join([">"+gene, str(gene_low), str(gene_high), num_strand, annot])
					
					gene_low -= 1 #Adjust to zero index
					
					#Collect reference sequence
					gene_seq = seq[gene_low:gene_high]
					
					nt_seq, aa_seq = self.gene_to_aa(gene_seq, strand, this_translation_table)
					
					if do_nt:
						nt_seqs.append(gene_printout)
						nt_formatted = self.format_bp_per_line(nt_seq, 70)
						nt_seqs.append(nt_formatted)
						
					if do_aa:
						aa_seqs.append(gene_printout)
						aa_formatted = self.format_bp_per_line(aa_seq, 60)
						aa_seqs.append(aa_formatted)
		
		output_filename_genomes = os.path.normpath(self.genomes_prefix + "/"+genome_name+"_genomic_fasta.txt")
		output_filename_nt = os.path.normpath(self.protein_nt_prefix + "/"+genome_name+"_gene_nt_fasta.txt")
		output_filename_aa = os.path.normpath(self.protein_aa_prefix + "/"+genome_name+"_gene_aa_fasta.txt")
		
		if do_genome:
			genome_seqs = '\n'.join(genome_seqs) + "\n"
			self.write_to_file(output_filename_genomes, genome_seqs)
			
		if do_nt:
			nt_seqs = '\n'.join(nt_seqs) + "\n"
			self.write_to_file(output_filename_nt, nt_seqs)
		
		if do_aa:
			aa_seqs = '\n'.join(aa_seqs) + "\n"
			self.write_to_file(output_filename_aa, aa_seqs)		

				

	
	