#!/usr/bin/env python
"""
nucleoHMM v0.1
Andre Rendeiro <afrendeiro at gmail.com> - 2014

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.
"""

import os, sys, csv, logging, pickle, json
from optparse import OptionParser
from math import *


def main():
    # option parser
    usage = 'python NucleoHMM.py [-m -o -l] seq.tsv'
    parser = OptionParser(usage=usage)

    parser.add_option("-m", "--model",
    type="str", dest="model", default='single',
    help="Nucleotide model to use. Provided options: 'single'; 'dinucleotide'; Default: 'single'")

    parser.add_option("-t", "--train",
    type="str", dest="train",
    help="Tab-delimited data file to train model.")

    parser.add_option("--model-output",
    type="str", dest="modeloutput", default="user",
    help="File to ouput trained model. Default to: models/user")

    parser.add_option("--train-only",
    type="str", dest="trainonly",
    help="Train model only and exit.")

    parser.add_option("-o", "--outputDir",
    type="str", dest="outdir", default='output/',
    help="""Directory for output files. Default: to: output/
    Files will be named as the fasta identifier, if provided.
    If no will be named based on sequence number with the following prefix: NucleoHMM_seq<seqnumber>.tsv""")

    parser.add_option("-l", "--log",
    type="str", dest="log", default='logs/log.txt',
    help="directory/filename to store log file to: default = logs/log.txt")

    # read arguments and options
    (options, args) = parser.parse_args()
    if len(args) > 2 or len(args) == 0:
        # return help mesage if argument number is incorrect
        print(__doc__)
        parser.print_help()
        sys.exit(0)
    
    # log configuration
    logging.basicConfig(filename=options.log,level=logging.DEBUG)

    # prepare directories for logs and output
    if not os.path.exists(os.path.dirname(options.log)):
        os.makedirs(os.path.dirname(options.log))
    if not os.path.exists(os.path.dirname(options.outdir)):
        os.makedirs(os.path.dirname(options.outdir))

    # parse in/outfile argument
    infile = args[0]
    outdir = options.outdir

    ####### Train model on data
    if options.train:
        trainingData = openData(options.train)

        transitionProb, emissionProb = trainModel(trainingData)
        
        # Add initial parameters to model (assuming start always on 'Not-bound' state!)
        states = ('Bound', 'Not-bound')
        startProb = {'Bound': 0.5, 'Not-bound': 0.5}

        saveModel("models/" + options.modeloutput, states, startProb, transitionProb, emissionProb)
        options.model = "user"

        if options.trainonly:
            print("Program finished. Trained model with provided sequences")
            sys.exit(0)

    # choose model to work with
    if options.model == "single":
        states, startProb, transitionProb, emissionProb = loadModel("models/" + options.model)
        logging.info('Chose the single nucleotide model')
    elif options.model == "dinucleotide":
        states, startProb, transitionProb, emissionProb = loadModel("models/" + options.model)
        logging.info('Chose the single dinucleotide model')
    elif options.model == "user":
        states, startProb, transitionProb, emissionProb = loadModel("models/" + options.modeloutput)
        logging.info('Chose the nucleotide model trained by user')
    else:
        states, startProb, transitionProb, emissionProb = loadModel("models/" + options.model)
        logging.info('Chose a nucleotide model previously trained by user')

    ####### Start predicting
    seqs = openSeqs(infile)
    i = 1
    skip = 0
    
    for ID, SEQ in seqs:
        ## preprocess sequence
        SEQ = SEQ.strip().rstrip().upper()
                
        # make list with string sequence
        SEQ = list(SEQ)

        # process sequence according to nucleotide model
        if options.model == "dinucleotide":
            newseq = []
            for n in range(0,len(SEQ)-1):
                newseq.append(SEQ[n] + SEQ[n+1])
            SEQ = newseq

        ## run viterbi
        path, probability = viterbi(states, startProb, transitionProb, emissionProb, SEQ)

        logging.info('Found Viterbi path for sequence %s with probability %e.' % (ID, probability))

        ## write output to files
        # add header
        out = [["Position", "State"]]
        # add remaining positions and state
        for pos in range(len(path)):
            out.append([SEQ[pos], path[pos]])
        
        if ID == "":
            writeData(out, outdir + "NucleoHMM_seq" + str(i) + ".tsv")
        else:
            writeData(out, outdir + "/" + ID + ".tsv")
        
        i += 1

    logging.info('Program finished. Processed %d sequences' % (i - 1 - skip))
    print("Program finished. Processed %d sequences" % (i - 1 - skip))

def openData(infile):
    """ Opens data file 'infile' with nucleossome data to train. Returns list of tuples from all data lines in file."""
    with open(infile) as data_file:
        i = 0
        lines = []
        for line in data_file:
            # skip header line
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

def trainModel(data):
    """ Trains the model with provided data. Outputs the model transition and emission probabilities."""
    # Initialize model structure with 0 counts
    transitionProb = {'Bound' : {'Bound' : 0, 'Not-bound' : 0}, 'Not-bound' : {'Bound' : 0, 'Not-bound' : 0}}
    emissionProb = {'Bound' : {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}, 'Not-bound' : {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}}
    
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
    with open(infile) as data_file:
        logging.info('Succesfully opened %s' % infile)
        lines = data_file.readlines()
    data_file.close()
    # Fasta parser:
    alphabet = ["A", "C", "T", "G"]
    seqs = []
    found = False

    def checkAlphabet(seq, line, alphabet=alphabet):
        """ Aborts program if 'seq' contains elements not in 'alphabet'"""
        for n in seq.rstrip():
            if n not in alphabet:
                logging.debug("Sequence on line %d contains non-nucleotide characters" % line)
                raise TypeError("Sequence on line %d contains non-nucleotide characters" % line)
                sys.exit(0)

    for line in range(0, len(lines)):
        # skip header lines (not starting with ">")
        if (lines[line][0] != ">" and found == False):
            continue
        # get sequence identifier
        elif (lines[line][0] == ">"):
            # prevent two consecutive ID lines
            if (lines[line + 1][0] != ">"):
                # extract ID: split lines in spaces, get first (identifier) and get rid of the '>' sign
                ID = lines[line].rstrip().split(" ")[0].replace(">", "")
                # mark finding of ID line
                found = True
                # initialize empty seq
                SEQ = ""
            else:
                raise Exception("Input file malformated. Two consecutive identifier lines found at line %d." % (line + 1))
                logging.info("Input file malformated. Two consecutive identifier lines found at line %d." % (line + 1))
                sys.exit(0)
        # sequence lines
        elif found:
            # if not last line
            if (line != (len(lines) - 1)):
                # if next line empty or ID line, append to seqs
                if (lines[line + 1] != "" and lines[line + 1][0] != ">"):
                    # just add line to sequence
                    checkAlphabet(lines[line], line + 1)
                    SEQ += lines[line].rstrip()
                else:
                    # add line to sequence and append to seqs. Finish seq
                    checkAlphabet(lines[line], line + 1)
                    SEQ += lines[line].rstrip()
                    seqs.append((ID, SEQ))
                    found = False
            else:
                # If last line just add to sequence and append to seqs. Finish seq
                checkAlphabet(lines[line], line + 1)
                SEQ += lines[line].rstrip()
                seqs.append((ID, SEQ))
                found = False
    return seqs

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