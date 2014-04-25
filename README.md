nucleoHMM
=======

Hidden-Markov Model-based nucleossome positioning prediction
-----------

# Running
This program has been tested in Python 3 and does not require any additional packages.

List all options:
    `python3 nucleoHMM.py`

Example usage to predict nucleossome positioning for all sequences in input file:

    `python3 nucleoHMM.py tests/seq.tsv`

# Training the model
You can train the model using with nucleotide positions and bound states in a tab-delimited file.

See an example in tests/train.tsv.

You can train the model by specifying the 't' option and giving a training file as argument:
    `python3 nucleoHMM.py -t tests/train.tsv tests/seq.tsv`

# Input
In the current state, nucleoHMM accepts files with one sequence per line containing only the DNA alphabet (A,C,G,T).

In the future will accept fasta format as input.

# Output
nucleoHMM outputs one tab-delimited file for each sequence with nucleotide positions and bound states.

When training model on data, nucleoHMM saves the model in binary format but also in human-readable JSON.

# Other options
Run `python3 nucleoHMM.py` to see all options.