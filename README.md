nucleoHMM
=======
### Hidden Markov Model-based nucleossome positioning prediction

# Running
This program has been tested in Python 3 and does not require any additional packages.

Example usage to predict nucleossome positioning for all sequences in input file:

    `python3 nucleoHMM.py tests/seq.fa`

# Training the model
nucleoHMM provides two naïve models based on single- and dinucleotide frequencies.

Currently, training the model is only implemented for the single nucleotide model. You can train using a file with nucleotide positions and bound states in a tab-delimited file (see Output and an example in tests/yeast_chr1.tsv).

Train the model with the 't' option and a training file as argument:

    `python3 nucleoHMM.py -t tests/train.tsv tests/yeast_chr1.tsv`

## Generating random training set
Supply an integer argument to the '-g' option to generate random training data with that length:

    `python3 nucleoHMM.py -g 10000 tests/seq.fa`

# Input
nucleoHMM accepts [FASTA](http://en.wikipedia.org/wiki/FASTA_format) files as input. Currently it supports only the DNA alphabet (A,C,G,T).

# Output
nucleoHMM outputs one tab-delimited file for each sequence with nucleotide positions and bound states.

Output files are named as the fasta identifier. If empty, files will be named based on sequence number with the following prefix: 'NucleoHMM_seq<seqnumber>.tsv'

Output directory can be specified with the -o option.

When training model on data, nucleoHMM saves the model in binary format but also in human-readable JSON.

# Other options
Run `python3 nucleoHMM.py` to see all options.
