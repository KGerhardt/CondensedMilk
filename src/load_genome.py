from quick_read import agnostic_reader
import os
import numpy as np

class genome_loader:
	def __init__(self, genome_file):
		self.fp = genome_file
		self.exists = os.path.exists(genome_file)
		
		self.file_basename = os.path.basename(genome_file)
		
		self.sequence_order = None
		self.sequences = None
		self.descriptions = None
		
		self.noncalls_replaced_with_A_indices = None
		#ASCII A = 65, C = 67, G = 71, T = 84
		self.acceptable_chars = np.array([65, 67, 71, 84], dtype = np.int32)
		
	def load_genome(self):
		sequence_order = []
		sequences = {}
		descriptions = {}
		noncalls_replaced_with_A_indices = {}
		genome_length = 0
		
		reader = agnostic_reader(self.fp)
		contents = reader.read()
		contents = contents.split(">")
		
		#Properly formatted FASTA files will have a first item '' after split on '>'. 
		#This is a string of length 0 prepending the first '>' character that is the first char in such a file.
		contents = contents[1:]
		for c in contents:
			lines = c.splitlines()
			defline = lines[0]
			seq = ''.join(lines[1:])
			seq = seq.upper() #Consistency
			seq = seq.replace('\x00', '')
				
			seqid = defline.split()[0]
			description = defline[len(seqid):]
			
			#Check for Ns or other non-call characters
			if set(seq) > set("ATCG"):
				#Convert the sequences' characters to a list of ASCII number representations
				as_ascii = np.array(list(seq.encode("ascii")), dtype = np.int32)
				#Find the loci of non-ATCG characters
				bad_char_loci = np.nonzero(np.logical_not(np.isin(as_ascii, self.acceptable_chars)))[0]
				#Replace those characters with "A" and record the locations
				as_ascii[bad_char_loci] = 67
				
				bad_char_loci = bad_char_loci.astype(np.int32)
				
				#Convert to bytestrings for insertion into db. It's easier and faster to load these with for numpy.frombuffer()
				noncalls_replaced_with_A_indices[seqid] = bad_char_loci.tobytes(order='C') #Explicit about the order to avoid endianness issues
				seq = b''.join(as_ascii) #turn it into a bytestring
				seq = seq.decode("ascii") #Convert back to text
				seq = seq.replace('\x00', '')
			else:
				#SQLite encodes these as nulls - clear indication a seq had no no-calls
				noncalls_replaced_with_A_indices[seqid] = None
				
			sequence_order.append(seqid)
			sequences[seqid] = seq
			descriptions[seqid] = description
			genome_length += len(seq)
			
		return genome_length, sequence_order, sequences, descriptions, noncalls_replaced_with_A_indices