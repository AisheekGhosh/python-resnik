#!/usr/bin/env/ python
from logging import root
import sys

# import pandas as pd
from pathlib import Path
from unittest.mock import DEFAULT


root_directory = Path.cwd()
if str(root_directory).endswith('python'):
    root_python_directory = root_directory
    root_project_directory = root_directory.parent
else:
    root_python_directory = Path.joinpath(root_directory, "python")
    root_project_directory = root_directory


print(root_directory)
fss_directory = Path.joinpath(root_python_directory, "fastsemsim-code", "fastsemsim")

fss_directory_str = str(fss_directory)
if fss_directory_str not in sys.path:
    sys.path.insert(0,fss_directory_str)


from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.SemSim.ObjSetSemSim import ObjSetSemSim
import fastsemsim.SemSim as SemSim


from collections import defaultdict
from decimal import Decimal

from subprocess import call
import os
import argparse



ONTOLOGY_CONVERT = {"BP" : "Process", "MF" : "Function", "CC" : "Component"}
ONTOLOGY_ROOT = {"BP" : 8150, "MF" : 3674, "CC" : 5575}
TAXON_CODES = {"AThaliana","CElegans","DMelanogaster","HSapiens", "MMusculus", "RNorvegicus",
               "SCerevisiae","SPombe"}
TAXON_PATH = str(Path.joinpath(root_python_directory, "taxons")) + '/'
TAXON_FILES = {key : TAXON_PATH + key + ".txt" 
                for key in TAXON_CODES }
GRAPH_FILE_PATH = str(Path.joinpath(root_python_directory, "networks")) + '/'
DEFAULT_SANA_OUT = str(Path.joinpath(root_python_directory, 'sana.out'))
GRAPH_A_NAME = 'CElegans'
GRAPH_B_NAME = 'MMusculus'

def getGOterms(ac, gene):
#    result = []
#    for geneName in gafCommon[gene]:
#        if geneName in ac.annotations:
#            result += list(ac.annotations[geneName].keys())
    try:
        return list(ac.annotations[gene].keys())
    except KeyError:
#    print(gene, "not in ac.annotations")
        return list

def geneList(fileName):
    with open(fileName,'r') as f:
        # Skips first 4 lines of a LEDA graph file
        for i in range(4):
            f.readline()
        count = int(f.readline().strip())
        return [f.readline().strip()[2:-2] for i in range(count)]

def MNE(geneListA, commonA, geneListB, commonB, ac):
    count = dict()
    total = 0
    for geneA, geneB in zip(geneListA, geneListB):
        total += 2
        a = getGOterms(ac, commonA, geneA)
        b = getGOterms(ac, commonB, geneB)
        for termA in a:
            try:
                count[termA] += 1
            except KeyError:
                count[termA] = 1
        for termB in b:
            try:
                count[termB] += 1
            except:
                count[termB] = 1
    frac = lambda x : Decimal(count[x]) / Decimal(total)
    result = sum(map(lambda x : frac(x) * frac(x).log10() , count)) * Decimal(len(count)).log10()
    return result 

def getAlignment(filename):
    with open(filename,'r') as f:
        # return map(int, f.readline().strip().split())
        return map (int, f.readline().strip().split())

def sortGAF(filename, aspects, output):
    sortString = "{{if ({}) print}}"
    arg = " || ".join('$8=="{}"'.format(aspect) for aspect in aspects)
    with open(output,'w') as f:
        call(["awk", "-F", "\t", sortString.format(arg), filename], stdout=f)

def initParser(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-g1", "--graph1", required=False, type=str, 
            help="graph1, which graph2 is aligned to", default=GRAPH_A_NAME + '.gw')
    parser.add_argument("-g2", "--graph2", required=False, type=str,
            help="graph2", default=GRAPH_B_NAME + '.gw')
    parser.add_argument("-t", "--taxon", type=str, 
            help="Taxonomies to include, separated by a '-'", default="AThaliana-SPombe")
    parser.add_argument("-so", "--sana-out", type=str,
            help="Output file", default=DEFAULT_SANA_OUT)
    parser.add_argument("-a", "--alignment-file", type=str,
            help="Two column alignment file"
            # , default="/github_repos/python-resnik/python/new_nets/AT_SC_path_3.align"
            )
    parser.add_argument("--mode", default="max", type=str,
            choices=["BMA","avg","max"], help="Evaluation mode")
    parser.add_argument("-ec", "--evidence-codes", type=str, help="Evidence Codes") 
    parser.add_argument("-i", "--include-ec", action="store_true", 
            help="Include evidence codes, default set to exclude")
    parser.add_argument("-oc", "--ontology-root", default="BP-MF-CC", type=str, 
		    help="Ontology root, defaults to all. Separate by '-'") 
    parser.add_argument("-o", "--output", type=str, 
            help="Output file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose")
    parser.add_argument("--semsim", default="Resnik", type=str,
            choices=[x for x in SemSim.term_SemSim_measures],
            help="Method of evaluating similarity")
    return parser

def getTaxons(name):
    with open(TAXON_FILES[name], mode='r') as f:
        return [line.strip() for line in f]



if __name__ == "__main__":

    GENE_ONTOLOGY_FILE = str(Path.joinpath(root_project_directory, "go.obo"))
    acSourceType = "ncbi"
    parser = initParser()
    args = parser.parse_args()
    ontologyRoots = {ONTOLOGY_ROOT[root] for root in args.ontology_root.split("-")}
    evidenceCodes = set(args.evidence_codes.split("-")) if args.evidence_codes else {}
    if hasattr(args, "ontology_codes"):
        if args.ontology_codes:
            ontologyCodes = set(ONTOLOGY_CONVERT[code] for code in args.ontology_codes.split("-"))
        else:
            ontologyCodes = ONTOLOGY_CONVERT.values()
    else:
        ontologyCodes = ONTOLOGY_CONVERT.values()
    taxonCodes = []

    if args.taxon:
        for taxon in args.taxon.split("-"):
            taxonCodes.extend(getTaxons(taxon)) 
    else:
        taxon1 = os.path.basename(args.graph1)[:-3]
        taxon2 = os.path.basename(args.graph2)[:-3]
        taxonCodes.extend(getTaxons(taxon1) if taxon1 in TAXON_CODES else [])
        taxonCodes.extend(getTaxons(taxon2) if taxon2 in TAXON_CODES else [])

    taxonCodes = set(taxonCodes)
    taxonInclude = True if len(taxonCodes) else False

    # Display all parameters
    if args.verbose:
        verboseFormat = "{:20s} :  {}"
        print(verboseFormat.format("Graph 1", args.graph1))
        print(verboseFormat.format("Graph 2", args.graph2))
        print(verboseFormat.format("Annotation Corpus", "gene2go"))
        print(verboseFormat.format("Alignment file", args.alignment_file))
        print(verboseFormat.format("Ontology file", GENE_ONTOLOGY_FILE))
        print(verboseFormat.format("Include taxons", taxonCodes))
        choice = "{} evid. code".format(("Exclude","Include")[args.include_ec])
        print(verboseFormat.format(choice, evidenceCodes))
        print(verboseFormat.format("Include ontologies", ontologyCodes))
        print(verboseFormat.format("Output File", args.output))
        print(verboseFormat.format("Semantic Similarity", args.semsim))
        print(verboseFormat.format("Mode", args.mode))

    acFile = str(Path.joinpath(root_project_directory, "gene2go")) 
    acParams = {}
    acParams["filter"] = {} # filter section to remove undesired annotations
    acParams["filter"]["EC"] = {"EC" : evidenceCodes, "inclusive" : args.include_ec}
    acParams["filter"]["taxonomy"] = {"taxonomy" : taxonCodes, "inclusive" : taxonInclude}

    if args.sana_out: 
        geneListA = geneList(GRAPH_FILE_PATH + os.path.basename(args.graph1)[:-3] + "/" + args.graph1) 
        geneListB = geneList(GRAPH_FILE_PATH + os.path.basename(args.graph2)[:-3] + "/" + args.graph2) 

        alignment = getAlignment(args.sana_out)
        # compareList = [geneListB[i] for i in alignment]
        compareList = [geneListB[i] for i,x in enumerate(alignment)]
    
    if args.alignment_file: 
        geneListA = []
        compareList = []
        with open(args.alignment_file, mode='r') as f:
            for line in f:
                geneA,geneB = line.strip().split()
                geneListA.append(geneA)
                compareList.append(geneB)

    ontology = ontologies.load(GENE_ONTOLOGY_FILE)
    ac = AnnotationCorpus.AnnotationCorpus(ontology)
    ac.parse(acFile, acSourceType, acParams)
    ac.isConsistent()
    scoreFunction = ObjSetSemSim(ontology, ac, args.semsim, args.mode)
    score_list = []

    if args.output:
        with open(args.output,'w') as f:
            for geneA,geneB in zip(geneListA, compareList):
                scores = [scoreFunction.SemSim(geneA, geneB, root) 
                        for root in ontologyRoots]
                maxScore =  max(score for score in scores if score != None)
                f.write("{}\t{}\t{}\t{}\r\n".format(args.mode,geneA,geneB,maxScore))
    else:
        for geneA,geneB in zip(geneListA, compareList):
            
            # StarStriker1478 = [scoreFunction.SemSim(geneA, geneB, root) for root in ontologyRoots]
            # score =  max(scoreFunction.SemSim(geneA, geneB, root) 
            #         for root in ontologyRoots)
            scores = [scoreFunction.SemSim(geneA, geneB, root) 
                    for root in ontologyRoots]
            try:
                score = max(x for x in scores if x != None)
            except:
                score = None
            # score =  max(scoreFunction.SemSim(geneA, geneB)) # 4 required arguments for SemSim function
            print(args.mode,geneA,geneB,score,sep="\t")
            score_list.append((args.mode,geneA,geneB,score))

    # # print(score_list)
    # score_df = pd.DataFrame(score_list, columns =['Mode', 'GeneA', 'GeneB', 'Score'])
    # # print(score_df)

    # filepath = Path('/github_repos/python-resnik/outputs/' + GRAPH_A_NAME + '_' + GRAPH_B_NAME +'.csv')  
    # filepath.parent.mkdir(parents=True, exist_ok=True)  
    # score_df.to_csv(filepath) 
