#!/usr/bin/env/ python
import sys
import argparse
from collections import defaultdict

BIOMART_DIR = "/extra/wayne_scratch0/preserve/resnik/biomart/"
ENSEMBL_TO_ENTREZ = BIOMART_DIR + "{}_ensembl2entrez"
ORTHOLOG_FILE = BIOMART_DIR + "{}_{}"
SPECIES_IMPLEMENTED = ["SCerevisiae", "HSapiens", "DMelanogaster", "CElegans",
        "MMusculus", "RNorvegicus"]

def initParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g1", required=True, type=str)
    parser.add_argument("-g2", required=True, type=str)
    parser.add_argument("-a", "--alignment-file", required=True, type=str)
    return parser

def ensemble2entrez(filename):
    species = dict()
    try:
        f = open(ENSEMBL_TO_ENTREZ.format(filename), mode='r')
    except FileNotFoundError:
        print("{} does not exist. Please try again.", file=sys.stderr)
        sys.exit()
    else:
        for line in f:
            try:
                ensembleID,entrezID = line.strip().split()
            except ValueError:
                continue
            species[entrezID] = ensembleID
        f.close()
    return species

def loadOrtholog(speciesA, speciesB):
    ortholog = defaultdict(dict)
    try:
        f = open(ORTHOLOG_FILE.format(speciesA, speciesB), mode='r')
    except FileNotFoundError:
        print("{}-{} comparison does not exist.".format(speciesA, speciesB),
                file=sys.stderr)
        sys.exit()
    else:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 3:
                continue
            speciesA,speciesB,confidence = fields
            ortholog[speciesA][speciesB] = confidence
        f.close()
    return ortholog

def isOrtholog(geneA, geneB, speciesA, speciesB, ortholog):
    try:
        ensembleA = speciesA[geneA]
        ensembleB = speciesB[geneB]
    except KeyError:
        return None
    if ensembleB in ortholog[ensembleA]:
        return ortholog[ensembleA][ensembleB]
    else:
        return -1

if __name__ == '__main__':
    parser = initParser()
    args = parser.parse_args()
    speciesA = ensemble2entrez(args.g1.lower())
    speciesB = ensemble2entrez(args.g2.lower())
    ortholog = loadOrtholog(args.g1.lower(), args.g2.lower())

    geneListA = []
    nodeListB = []

    with open(args.alignment_file, mode='r') as f:
        for line in f:
            geneA,geneB = line.strip().split()
            geneListA.append(geneA)
            nodeListB.append(geneB)
    
    for geneA, geneB in zip(geneListA, nodeListB):
        score = isOrtholog(geneA, geneB, speciesA, speciesB, ortholog)
        print(geneA, geneB, score, sep="\t")

