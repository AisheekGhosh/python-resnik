# RESNIK

## What is Resnik?

Resnik is a scoring system used for aligning two different networks. The scoring returns a value from 0 to around 14, based on how similar the two alignments are. It can also return "none", if there is not enough information.

# REQUIREMENTS

Miniconda environment, with Python 3.6+

# Reading Output

A file is created in a folder called outputs/, which contains all the files of alignments and scorings. The output is a .csv file, with geneA as a column, geneB as a column, and the score.

# HOW TO SET UP

## Windows

set up this way:
* Download Miniconda 3
* Run this command on base environment:        

>     cd [location of your python-resnik]
>     conda env create -f py39_resnik.yml
* Activate conda environment `py39_resnik`
* Run 
        
>     python resnik4.py -g1 [graph 1] -g2 [graph 2] -so [sana.out file]

## Mac/Linux

set up this different way:
* altlist 1
* altlist 2
* altlist 3

--------
NOTE: the module used in this code, fastsemsim, can be updated to a 2019 version. However, for the purposes of this project, this module is not required as it is already implemented in the file through the use of directories

## Relevant Links
* https://docs.conda.io/en/latest/miniconda.html - Miniconda 3 Download
* https://arxiv.org/abs/1607.02642 - SANA Paper
* https://sana.ics.uci.edu - Official SANA Website
* https://github.com/waynebhayes/SANA/tree/SANA2/doc - SANA GitHub Page
* https://sourceforge.net/projects/fastsemsim/ - Fastsemsim Download
* https://sites.google.com/site/fastsemsim/ - Fastsemsim Official Site

<!-- Made by Aisheek Ghosh -->