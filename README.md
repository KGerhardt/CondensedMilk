# easygenes
Tool for efficiently predicting genes from genomes and storing those genome and gene sequences

Easygenes is a tool for taking a set of prokaryotic genome sequences, predicting genes for those sequences, and storing both in an efficient database that allows for their retrieval in part and whole. You give it a set of genomes, and it gives you a database storing those genomes and the coordinates - but not the sequences - of their genes. In reverse, you can rapidly recover the original genome sequences or use the coordinates of the genes to retrieve them in either nucleotide or animo acid formats. Easygenes also supports gzipped inputs and gzipping outputs for the sake of sanity.

Internally, it utilizes GenomicSQLite (https://github.com/mlin/GenomicSQLite) to store DNA sequences in a twobit representation and in a compressed database format to reduce the overall footprint of the database and accelerate the loading time of DNA sequences, Pyrodigal (https://github.com/althonos/pyrodigal/) to manage gene prediction, and PyHMMER (https://github.com/althonos/pyhmmer/) to efficiently manage nucleotide to amino acid translation.

Currently this project supports two actions:

(1) python3 src/main.py condense -g [genome_directory] -db [database_path] -p [threads]
(2) python3 src/main.py retrieve -db [database_path] -o [output_directory] --compress

(1) handles building a database and predicting genes, (2) handles retrieving those genes.

This is still a work in progress, but this is in a state with some useful functionality.
