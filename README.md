# TSAcore

The code in this folder contains the core TSAFinder algorithm. It has been tested in python 2 and 3. 

# Prerequisites
# Python packages:
os  \ 
numpy \
sys \
subprocess \
multiprocessing \
itertools \
linecache
# Software:
netMHCpan (tested with version 4.0)

# Inputs
Tumor fastq file \
Control fastq file

# Usage
# Arguments:
1. kmer length (currently supports 8-11) 
2. getmeta. Set to 'getmeta' to produce the metadata, which is a file with each tumor specific peptide and its associated rnaseq read. Set to 'FALSE' to exclude this step. 
3. nfiles. Number of parallel processes (cores) used. Minimum of 4. Tested with 20. 
4. patient2processn. Control sample name. 
5. patient2processt. Tumor sample name. 
6. hlatype. HLA genotype. e.g. ['HLA-A33:01', 'HLA-A24:01', 'HLA-B54:01', 'HLA-B58:01', 'HLA-C03:01'] 

Enter the neoantigen_toolbox_v4.3.py file and change the following: \
readlen_tumor. The tumor sample read length. \
readlen_normal. The control sample read length. \
tumor_depth_constant. Specifies the expression threshold for the tumor sample. Default = 100. \
normal_depth_constant. Specifies the expression threshold for the control sample.  Default = 2. \
tumor_fpath. Specifies the path to the tumor fastq file. \
normal_fpath. Specifies the path to the control fastq file. 
