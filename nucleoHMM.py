#!/usr/bin/env python
"""
nucleoHMM v0.1
Andre Rendeiro <afrendeiro at gmail.com> - 2014

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.
"""

import os, sys, csv, logging, pickle, json, pprint, random
from optparse import OptionParser
from math import *


def main():
    # option parser
    usage = 'python NucleoHMM.py [-m -g -t -o -l] seq.tsv'
    parser = OptionParser(usage=usage)

    parser.add_option("-m", "--model",
    type="str", dest="model", default='single',
    help="Nucleotide model to use. Options: 'single' or 'dinucleotide'; Default: 'single'")

    parser.add_option("-g", "--generateRandom",
    type="int", dest="generate",
    help="Generates random training set with length of argument (int) and exits. File written in 'models/random'")

    parser.add_option("-t", "--train",
    type="str", dest="train",
    help="Data file to train model.")

    parser.add_option("--model-output",
    type="str", dest="modeloutput", default="user",
    help="File to ouput trained model. Default to: models/user")

    parser.add_option("--train-only",
    dest="trainonly", action="store_true",
    help="Flag to only train model and then exit.")

    parser.add_option("-o", "--outputDir",
    type="str", dest="outdir", default='output/',
    help="""Directory for output files. Default: to: output/
    Files will be named as the fasta identifier, if provided.
    If empty, files will be named based on sequence number with the following prefix: NucleoHMM_seq<seqnumber>.tsv""")

    parser.add_option("-l", "--log",
    type="str", dest="log", default='logs/log.txt',
    help="directory/filename to store log file to: default = logs/log.txt")

    # read arguments and options
    (options, args) = parser.parse_args()

    if len(args) > 1 or (len(args) == 0 and options.trainonly == None and type(options.generate) != type(int())):
        # return help mesage if argument number is incorrect
        print(__doc__)
        parser.print_help()
        sys.exit(0)
    
    # log configuration
    logging.basicConfig(filename=options.log,level=logging.DEBUG)

    # prepare directories for logs and output
    if not os.path.exists(os.path.dirname(options.log)):
        os.makedirs(os.path.dirname(options.log))
        logging.info("Created directory to store logs: %s" % options.log)
    if not os.path.exists(os.path.dirname(options.outdir)):
        os.makedirs(os.path.dirname(options.outdir))
        logging.info("Created directory for output: %s" % options.outdir)

    global alphabet
    alphabet = ["A", "C", "T", "G"]
    
    global states
    states = ('Bound', 'Not-bound')
    
    ####### Generate random training data
    if options.generate:
        generateData('tests/random.tsv', options.generate, alphabet=alphabet, states=states)
        logging.info("Program finished. Generated random training data. Output: 'tests/random.tsv'.")
        print("Program finished. Generated random training data. Output: 'tests/random.tsv'.")
        sys.exit(0)

    ####### Train model on data
    if options.train:
        print("Training model with provided data...")
        # make sure training file exist
        if not os.path.isfile(options.train):
            logging.info("Training file %s does not exist. Aborting." % options.train)
            print("Training file %s does not exist. Aborting." % options.train)
            sys.exit(0)

        trainingData = openData(options.train)

        transitionProb, emissionProb = trainModel(trainingData, options.model)
        
        # Add initial parameters to model (assuming start always on 'Not-bound' state!)
        startProb = {'Bound': 0.5, 'Not-bound': 0.5}

        saveModel("models/" + options.modeloutput, states, startProb, transitionProb, emissionProb)

        # Print trained parameters
        print("Model trained.")
        print("Start probabilities:")
        pprint.pprint(startProb, indent=4)
        print("Transition probabilities:")
        pprint.pprint(transitionProb, indent=4)
        print("Emission probabilities:")
        pprint.pprint(emissionProb, indent=4)

        if options.trainonly:
            logging.info("Program finished. Trained %s nucleotide model with provided data. Saved in '%s'." % (options.model, "models/" + options.modeloutput))
            print("Program finished. Trained %s nucleotide model with provided data. Saved in '%s'." % (options.model, "models/" + options.modeloutput))
            sys.exit(0)

    ####### Predict using model
    # get infile argument
    infile = args[0]
    # make sure file exists
    if not os.path.isfile(infile):
        logging.info("Input file %s does not exist. Aborting." % infile)
        print("Input file %s does not exist. Aborting." % infile)
        sys.exit(0)
    
    seqs = openSeqs(infile)

    # if not trained, use preovided model, else is already loaded
    if not options.train:
        # load predefined model
        if options.model == "single":
            states, startProb, transitionProb, emissionProb = loadModel("models/" + options.model)
            logging.info('Chose the single nucleotide model')
        elif options.model == "dinucleotide":
            states, startProb, transitionProb, emissionProb = loadModel("models/" + options.model)
            logging.info('Chose the single dinucleotide model')

    i = 1
    
    for ID, SEQ in seqs:
        ## preprocess sequence
        SEQ = SEQ.strip().rstrip().upper()
                
        # make list with string sequence
        SEQ = list(SEQ)

        # process sequence according to nucleotide model
        if options.model == "dinucleotide":
            print("dinuc")
            newseq = []
            for n in range(0,len(SEQ)-1):
                newseq.append(SEQ[n] + SEQ[n + 1])
            SEQ = newseq

        ## run viterbi
        path, probability = viterbi(states, startProb, transitionProb, emissionProb, SEQ)

        logging.info('Found Viterbi path for sequence %s with probability %e.' % (ID, probability))

        ## write output to files
        # add header
        out = [["Position", "State"]]
        # add remaining positions and state
        for pos in range(len(path)):
            # process sequence according to nucleotide model
            if options.model == "single":
                out.append([SEQ[pos], path[pos]])
            elif options.model == "dinucleotide":
                print("dinuc")
                out.append([SEQ[pos][0], path[pos]])
        
        if ID == "":
            writeData(out, options.outdir + "NucleoHMM_seq" + str(i) + ".tsv")
            print("Processed sequence nº %d. Output:" % (i), "'" + options.outdir + "NucleoHMM_seq" + str(i) + ".tsv'")
        else:
            writeData(out, options.outdir + "/" + ID + ".tsv")
            print("Processed sequence nº %d. Output:" % (i), "'" + options.outdir + "/" + ID + ".tsv'")
        
        i += 1

    logging.info('Program finished. Processed %d sequences' % (i - 1))
    print("Program finished. Processed %d sequences" % (i - 1))

def generateData(outfile, length, alphabet, states):
    """ Creates random data that can be used to train the model."""
    with open(outfile, 'w') as tsvfile:
        wr = csv.writer(tsvfile, delimiter='\t')
        for i in range(0, length):
            wr.writerow((random.choice(alphabet), random.choice(states)))
    tsvfile.close()
    logging.info('Succesfully created random training data. Output: %s' % outfile)

def openData(infile):
    """ Opens data file 'infile' with nucleossome data to train. Returns list of tuples from all data lines in file."""
    with open(infile) as data_file:
        i = 0
        lines = []
        for line in data_file:
            # skip header line # fix this to detect presence of header
            if i == 0:
                i += 1
                continue
            if (line == "" or line == "\n"):
                raise Exception("Empty sequences in training data.")
                logging.info('Empty sequences in training data. Terminating')
                sys.exit(0)
            position, state = line.rstrip().split("\t")
            lines.append((position, state))
            i += 1
    return lines

def trainModel(data, model='single'):
    """ Trains the model with provided data. Outputs the model transition and emission probabilities."""
    # Initialize model structure with 0 counts
    transitionProb = {'Bound' : {'Bound' : 0, 'Not-bound' : 0}, 'Not-bound' : {'Bound' : 0, 'Not-bound' : 0}}

    if model == 'single':
        emissionProb = {'Bound' : {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}, 'Not-bound' : {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}}
    elif model == 'dinucleotide':
        emissionProb = {
            'Bound' : {'AA': 0, 'AC': 0, 'AT': 0, 'AG': 0, 'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0, 
                'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0, 'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0 },
            'Not-bound' : {'AA': 0, 'AC': 0, 'AT': 0, 'AG': 0, 'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0, 
                'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0, 'TA': 0, 'TC': 0, 'TG': 0, 'TT': 0 }
        }
        # process sequence according to nucleotide model
        newdata = []
        for t in range(0,len(data) - 1):
            newdata.append((data[t][0] + data[t + 1][0], data[t][1]))
        data = newdata
    
    # keep track of state occurance
    stateOccur = {'Bound' : 0.0, 'Not-bound' : 0.0}

    # force a start 'Not-bound' state
    prevState = 'Not-bound'
    stateOccur[prevState] += 1

    for t in range(0, len(data)):
        curPos = data[t][0]
        curState = data[t][1]
        stateOccur[curState] += 1
        
        # Count state transition occurence
        transitionProb[prevState][curState] += 1

        # Count state emission occurence
        emissionProb[curState][curPos] += 1

        # Update states for next t
        prevPos = curPos
        prevState = curState

    # force a end 'Not-bound' state
    curState = 'Not-bound'
    transitionProb[prevState][curState] += 1

    # Transform into probabilities (conditional)
    for prevState in transitionProb:
        for curState in transitionProb[prevState]:
            transitionProb[prevState][curState] /= stateOccur[prevState]

    for curState in emissionProb:
        for curPos in emissionProb[curState]:
            # compensate for extra end state
            if curState == 'Not-bound':
                emissionProb[curState][curPos] /= stateOccur[curState] - 1
            else:
                emissionProb[curState][curPos] /= stateOccur[curState]

    logging.info('Succesfully trained model with user data')
    return transitionProb, emissionProb

def saveModel(modelOutput, states, startProb, transitionProb, emissionProb):
    """ Saves user trained model as pickle binary object """
    with open(modelOutput, 'wb') as out:
        pickle.dump(states, out)
        pickle.dump(startProb, out)
        pickle.dump(transitionProb, out)
        pickle.dump(emissionProb, out)
    out.close()
    with open(modelOutput + ".json", 'w') as out:
        json.dump(states, out, indent = 4)
        json.dump(startProb, out, indent = 4)
        json.dump(transitionProb, out, indent = 4)
        json.dump(emissionProb, out, indent = 4)
    out.close()
    logging.info('Succesfully saved trained model to %s' % modelOutput)
    logging.info('Succesfully wrote trained model to %s' % modelOutput)

def loadModel(pickle_obj):
    """ Loads model from pickle binary file """
    with open(pickle_obj, 'rb') as model:
        states = pickle.load(model)
        startProb = pickle.load(model)
        transitionProb = pickle.load(model)
        emissionProb = pickle.load(model)        
    model.close()
    logging.info('Succesfully loaded %s as model' % pickle_obj)
    return states, startProb, transitionProb, emissionProb

def openSeqs(infile):
    """ Opens input file 'infile'. Returns list of lines from file."""

    def checkAlphabet(seq, alphabet=alphabet):
        """ Aborts program if 'seq' contains elements not in 'alphabet'"""
        for n in seq.rstrip():
            if n not in alphabet:
                logging.debug("Sequence contains non-nucleotide characters")
                raise TypeError("Sequence contains non-nucleotide characters")
                sys.exit(0)

    with open(infile) as data_file:
        logging.info('Succesfully opened %s' % infile)

        # Fasta parser
        while True:
            line = data_file.readline()
            if line == "":
                return  # Premature end of file, or just empty?
            if line[0] == ">":
                break

        while True:
            if line[0] != ">":
                raise ValueError("Records in Fasta files should start with '>' character")
            ID = line[1:].rstrip().split(" ")[0]
            seqs = []
            line = data_file.readline()
            while True:
                if not line:
                    break
                if line[0] == ">":
                    break
                checkAlphabet(line)
                seqs.append(line.rstrip())
                line = data_file.readline()

            yield ID, "".join(seqs).replace(" ", "").replace("\r", "")

            if not line:
                print("not line")
                return  # StopIteration

def viterbi(states, startProb, transitionProb, emissionProb, seq):
    """ Computes the Viterbi path of a sequence according to the given HMM model"""
    V = [{}]
    path = {}

    # Initialize base case (t == 0)
    for state in states:
        V[0][state] = startProb[state] * emissionProb[state][seq[0]]
        path[state] = [state]

    # Run Viterbi for t > 0
    for t in range(1, len(seq)):
        V.append({})
        newpath = {}
        for state in states:
            # calculate the probability of the new possible paths amd select the most probable
            prob, probState = max([(V[t-1][s] * transitionProb[s][state] * emissionProb[state][seq[t]], s) for s in states])
            # insert into Viterbi path
            V[t][state] = prob
            newpath[state] = path[probState] + [state]
        path = newpath

    finalState, pathProb = max((s, V[t][s]) for s in states)

    return path[finalState], pathProb

def writeData(data, outfile):
    """ Writes every element from the list 'data' into rows in file 'outfile' - tab-delimited."""
    with open(outfile, 'w') as tsvfile:
        wr = csv.writer(tsvfile, delimiter='\t')
        for i in range(0,len(data)):
            wr.writerow(data[i])
        tsvfile.close()
        logging.info('Succesfully wrote %s' % outfile)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Program canceled by user!\n")
        sys.exit(0)