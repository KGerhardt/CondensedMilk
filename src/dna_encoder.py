import genomicsqlite

#Compress sequences to bytestrings to save memory. Probably unnecessary but python pickling for memory passing is a pain with multiprocessing
class dna_encoder:
	def __init__(self):
		self.seq = None
		self.twobit = None
		
		
	#Creating a database in memory and encoding then retrieving the sequence that way isn't super 
	#efficient, but it's zero I/O and ensures that the sequence is ready to be inserted into a database 
	#as-is rather than requiring encoding in main
	def compose_twobit(self, seq_dict, seq_order):
		insertable = []
		for seq in seq_order:
			next_row = (seq_dict[seq],)
			insertable.append(next_row)
	
		conn = genomicsqlite.connect(":memory:")
		conn.execute("CREATE TABLE seqs (genome_seq BLOB)")
		if len(insertable) > 1:
			conn.executemany('INSERT INTO seqs VALUES (nucleotides_twobit(?))', insertable)
		else:
			conn.execute('INSERT INTO seqs VALUES (nucleotides_twobit(?))', insertable)
			
		conn.commit()
		
		genome_sequences = conn.execute("SELECT genome_seq FROM seqs").fetchall()
		
		conn.close()
		#List of tuples to list of blobs
		genome_sequences = [g[0] for g in genome_sequences]
		
		results = dict(zip(seq_order, genome_sequences))
		
		return results
		
	def decompose_twobit(self, twobit_seq, seq_order):
		insertable = []
		for seq in seq_order:
			next_row = (twobit_seq[seq],)
			insertable.append(next_row)
	
		conn = genomicsqlite.connect(":memory:")
		conn.execute("CREATE TABLE seqs (genome_seq BLOB)")
		if len(insertable) > 1:
			conn.executemany('INSERT INTO seqs VALUES (?)', insertable)
		else:
			conn.execute('INSERT INTO seqs VALUES (?)', insertable)
			
		conn.commit()
		
		genome_sequences = conn.execute("SELECT twobit_dna(genome_seq) FROM seqs").fetchall()
		
		conn.close()
		#List of tuples to list of blobs
		genome_sequences = [g[0] for g in genome_sequences]
		
		results = dict(zip(seq_order, genome_sequences))
		
		return results
		
		


