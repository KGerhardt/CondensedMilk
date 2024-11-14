# CondensedMilk
Tool for efficiently predicting genes from genomes and storing those genome and gene sequences

CondensedMilk is a tool for taking a set of prokaryotic genome sequences, predicting genes for those sequences, and storing both in an efficient database that allows for their retrieval in part and whole. You give it a set of genomes, and it gives you a database storing those genomes and the coordinates - but not the sequences - of their genes. In reverse, you can rapidly recover the original genome sequences or use the coordinates of the genes to retrieve them in either nucleotide or animo acid formats. Easygenes also supports gzipped inputs and gzipping outputs for the sake of sanity.

Internally, it utilizes GenomicSQLite (https://github.com/mlin/GenomicSQLite) to store DNA sequences in a twobit representation and in a compressed database format to reduce the overall footprint of the database and accelerate the loading time of DNA sequences, Pyrodigal (https://github.com/althonos/pyrodigal/) to manage gene prediction, and PyHMMER (https://github.com/althonos/pyhmmer/) to efficiently manage nucleotide to amino acid translation.

Currently this project supports four actions:

(1) python3 src/main.py condense -g [genome_directory] -db [database_path] -p [threads] [--sweeten]
(2) python3 src/main.py sweeten -db [database_path] -p [threads]
(3) python3 src/main.py simmer --donor [donor_database] --recipient [recipient_database]
(4) python3 src/main.py bake -db [database_path] -o [output_directory] [--compress]

(1) Handles compressing and adding genomes to a database, optionally predicting genes at the same time if --sweeten is used.
(2) Checks the database for any genomes which have not yet had genes predicted for them and runs gene prediction. Translation tables 4 and 11 are compared for performance for each genome and the table which predicts a higher coding density in the genome is selected.
(3) Allows a user to combine two Condensed Milk databases.
(4) Allows a user to retrieve (A) a summary of the contents of the database in JSON or tabular formats, and (B) to recover the original sequences of the genomes in FASTA format and to output genes in either nucleotide FASTA or amino acid FASTA formats, each on a per-genome basis. If --compress is used, the outputs will be gzipped.

The core behavior of this tool is done and it is ready to use on a command line, but I will be adding a Python API for key behaviors.
