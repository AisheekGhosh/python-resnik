#!/usr/bin/env/ python

def gafCommonTerm(fileName):
    result = dict() 
    gafFile = open(fileName,"r")
    for line in gafFile:
        if line[0] == '!':
            continue
        params = line.strip('\n').split('\t')
        result[params[2]] = params[1]
        for synonym in params[4].split("|"):
            result[synonym] = params[1]
    gafFile.close()
    return result

def geneList(fileName):
    with open(fileName,'r') as f:
        # Skips first 4 lines of a LEDA graph file
        for i in range(4):
            f.readline()
        count = int(f.readline().strip())
        return [f.readline().strip()[2:-2] for i in range(count)]

if __name__ == '__main__':
    a = gafCommonTerm('/extra/wayne_scratch0/preserve/resnik/gene_info')
    b = gafCommonTerm('/extra/wayne_scratch0/preserve/resnik/Saccharomyces_cerevisiae.gene_info')
    human = geneList('/extra/wayne_scratch0/preserve/resnik/networks/human/human.gw')
    yeast = geneList('/extra/wayne_scratch0/preserve/resnik/networks/yeast/yeast.gw')
    with open('human_ncbi.txt', mode='w') as f:
        for node in human:
            if node in a:
                print(node,a[node],sep="\t", file=f)
            else:
                print(node,"None", sep="\t",file=f)
    with open('yeast_ncbi.txt', mode='w') as f:
        for node in yeast:
            if node in b:
                print(node,b[node],sep="\t",file=f)
            else:
                print(node,"None", sep="\t",file=f)

