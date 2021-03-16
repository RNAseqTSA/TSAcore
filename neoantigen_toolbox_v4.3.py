#!/bin/python
import time
print('Start')
start_time = time.time()
import os
import numpy as np
import sys
import subprocess as sp
import multiprocessing as mp
import itertools as itertools
import linecache

def uuid_map(fname, cancer):
    with open(fname) as f:
        lines = f.read().splitlines()
    lines = [line.rstrip('\n').split('\t') for line in lines]
    patients = list()
    for line1 in lines:
        if line1[3] == 'Aligned Reads' and line1[7] == 'Solid Tissue Normal':
            for line2 in lines:
                if line2[3] == 'Aligned Reads' and line2[7] != 'Solid Tissue Normal' and line1[5] == line2[5]:
                    patients.append([line2[6], line2[0], line2[1], line1[0], line1[1]])
    fout = open(cancer+'_uuids.txt','w')
    for patient in patients:
        fout.write(patient[0]+'\t'+patient[1]+'\t'+patient[2]+'\t'+patient[3]+'\t'+patient[4]+'\n')
    return('Done')

def bam2seq(args):
    print('Converting bams to sequences')
    fname = args[0]
    oname = args[1]
    print('module load samtools/1.2; samtools view '+fname+' > '+fname[0:fname.find('.')]+'.sam')
    os.system('module load samtools/1.2; samtools view '+fname+' > '+fname[0:fname.find('.')]+'.sam')
    print('awk \'{print $3 " " $4 " " $10}\' '+fname[0:fname.find('.')]+'.sam > '+oname)
    os.system('awk \'{print $3 " " $4 " " $10}\' '+fname[0:fname.find('.')]+'.sam > '+oname)
    print('rm '+fname[0:fname.find('.')]+'.sam')
    #os.system('rm '+fname[0:fname.find('.')]+'.sam')
    print('Done: '+str(round(time.time() - start_time,4))+'s')
    os.system('module unload samtools/1.2')
    return('Done')

def fastq2seq(args):
    if args[0] == 'paired-end':
        fname1 = args[1]
        fname2 = args[2]
        oname = args[3]
        os.system('awk "(NR%4==2)" '+fname1+' > tmp1.txt')
        os.system('awk "(NR%4==2)" '+fname2+' > tmp2.txt')
        os.system('cat tmp1.txt tmp2.txt > '+oname)
        os.system('rm tmp1.txt tmp2.txt')
    else:
        fname1 = args[1]
        oname = args[2]
        os.system('awk "(NR%4==2)" '+fname1+' > '+oname)
    return('Done')

def bam2fastq(args):
    # convert bam file to fastqs (MIND-NUMBINGLY SLOW)
    fpath = args[0]
    fname_prefix = args[1]
    foutpath = args[2]
    line1 = 'module load gcc;'
    line2 = 'module load bedtools;'
    line3 = 'module load samtools/1.2;'
    line4 = 'samtools sort -n '+fpath+'/'+fname_prefix+'.bam '+fpath+'/aln.qsort;'
    line5 = 'bedtools bamtofastq -i '+fpath+'/aln.qsort.bam -fq '+foutpath+'/'+fname_prefix+'_rd1.fq -fq2 '+foutpath+'/'+fname_prefix+'_rd2.fq;'
    line6 = 'module unload gcc;'
    line7 = 'module unload bedtools;'
    line8 = 'module unload samtools/1.2;'
    os.system(line1+line2+line3+line4+line5+line6+line7+line8)

def fastq2hla(args):
    # Generating HLA predictions
    fpath = args[0]
    fname_prefix = args[1]
    tool_path = args[2]
    pat_num = args[3]
    kmer_k = args[4]
    if args[5] == 'paired-end':
        line1 = 'module unload PrgEnv-cray;'
        line2 = 'module load PrgEnv-intel;'
        line3 = 'module load R;'
        line4 = 'module load bowtie1;'
        line5 = 'export PATH=$PATH:'+tool_path+'/seq2HLA;'
        line6 = 'chmod u+x '+tool_path+'/seq2HLA/seq2HLA.py;'
        line7 = 'chmod u+x '+tool_path+'/seq2HLA/command.R;'
        line8 = 'export PYTHONPATH=$PYTHONPATH:'+tool_path+'/seq2HLA/biopython-1.58;'
        line9 = 'cd '+tool_path+'/seq2HLA; outputPrefix="hlaPred'+str(pat_num)+'_'+str(kmer_k)+'";'
        line10 = 'python seq2HLA.py -1 '+fpath+'/'+fname_prefix+'_R1.fastq -2 '+fpath+'/'+fname_prefix+'_R2.fastq -r $outputPrefix -l 100;'
        line11 = 'module unload R;'
        line12 = 'module unload PrgEnv-intel;'
        line13 = 'module load PrgEnv-cray;'
        #os.system(line1+line2+line3+line4+line5+line6+line7+line8+line9+line10+line11+line12+line13)
        os.system(line3+line4+line5+line6+line7+line8+line9+line10)
        # Extracting HLA predictions from output file
        HLAalleles = []
        fin = open(tool_path+'/seq2HLA/hlaPred'+str(pat_num)+'_'+str(kmer_k)+'-ClassI.HLAgenotype','r')
        line = fin.readline()
        line = fin.readline()
        while line != '':
            line = line.rstrip('\n')
            if 'hoz' in line:
                line = line.split('\t')
                tmp = line[1]
                tmp = tmp.split('*')
                HLAalleles.append('HLA-'+tmp[0]+tmp[1]+':01')
            else:
                line = line.split('\t')
                tmp = line[1]
                tmp = tmp.split('*')
                HLAalleles.append('HLA-'+tmp[0]+tmp[1]+':01')
                tmp = line[3]
                tmp = tmp.split('*')
                HLAalleles.append('HLA-'+tmp[0]+tmp[1]+':01')
            line = fin.readline()
        os.system('rm '+tool_path+'/seq2HLA/hlaPred'+str(pat_num)+'_'+str(kmer_k)+'*')
    return HLAalleles

def split_file(args):
    fname = args[0]
    nfiles = int(args[1])
    print('Splitting '+fname+' into '+str(nfiles)+' files')
    #Begin by counting the lines of input sequence file
    print('Counting lines')
    if os.path.isfile(fname+'.txt'):
        os.system('wc -l '+fname+'.txt > '+fname+'_tmp.txt')
    elif os.path.isfile(fname):
        os.system('wc -l '+fname+' > '+fname+'_tmp.txt')
    else:
        print('Error: function: split_file(): file does not exist')
    with open(fname+'_tmp.txt') as f:
        N = f.readline()
    os.system('rm '+fname+'_tmp.txt')
    N = N.split(" ")
    N = int(N[-2])
    if os.path.isfile(fname+'.txt'):
        os.system('split -l '+str(N/nfiles+1)+' '+fname+'.txt '+fname)
    else:
        os.system('split -l '+str(N/nfiles+1)+' '+fname+' '+fname)
    print('Done splitting '+fname+' into '+str(nfiles)+' files:'+str(round(time.time() - start_time,4))+'s \tFiles: '+str(nfiles))
    return('Done')

def merge_files(args): #merges the split files into single file with cat
    fpath = args[0]
    fname_prefix = args[1]
    nfiles = args[2]
    exten = args[3]
    alphabet = args[4]
    rmstat = args[5]
    print('Merging ' + fname_prefix + ' files into one file')
    i = 0;
    j = 0;
    fnames = [''] * nfiles
    if nfiles > 100:
        for i in range(0,20):
            for j in range(0,20):
                fnames[i*20+j] = fpath+'/'+fname_prefix+alphabet[i]+alphabet[j]+exten
    elif nfiles < 100:
        for i in range(0,20):
            fnames[i] = fpath+'/'+fname_prefix+alphabet[0]+alphabet[i]+exten
    command = 'cat'
    for i in range(nfiles):
        command = command + ' ' + fnames[i]
    command = command + ' > ' + fpath + '/' + fname_prefix + exten
    os.system(command)
    if rmstat == 'yes':
        command2 = 'rm'
        for i in range(nfiles):
            command3 = command2 + ' ' + fnames[i]
            os.system(command3)
    print('Done merging ' + fname_prefix + ' files into one file:'+str(round(time.time() - start_time,4))+'s \tFiles: '+str(nfiles))
    return('Done')

def translate(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

def write_kmers(line,meta,outfile,k,AAlist,iteration,i):
    idx = 0;
    ll = len(line);
    while idx+k-1<ll:
        if len(meta) == 1 and (('_' in line[idx:idx+k]) == False):
            if iteration == 0 and AAlist.index(line[idx]) < 10:
                outfile_temp = outfile[AAlist.index(line[idx])*20+AAlist.index(line[idx+1])]
                outfile_temp.write(str(i)+' '+line[idx:idx+k]+'\n');
            elif iteration == 1 and AAlist.index(line[idx]) >= 10:
                outfile_temp = outfile[AAlist.index(line[idx])*20+AAlist.index(line[idx+1])-200]
                outfile_temp.write(str(i)+' '+line[idx:idx+k]+'\n');
        idx = idx+1;
    return('Done')

def translate_file(args):
    fname_prefix_sub = args[0]
    k = args[1]
    offset = args[2]
    print('Translating '+fname_prefix_sub+' to '+str(k)+'mers ')
    infile = open(fname_prefix_sub, 'r')
    AAlist = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    for iteration in range(0,2):
        outfile = range(0,200)
        infile = open(fname_prefix_sub, 'r')
        for outcount in range(0,10):
            for outcount2 in range(0,20):
                outfile[outcount*20+outcount2] = open(fname_prefix_sub[:-2]+AAlist[outcount+(iteration*10)]+AAlist[outcount2]+fname_prefix_sub[-2:]+'_aminodat', 'w')
        i=offset[fname_prefix_sub[-1]];
        line = infile.readline()
        line = line.rstrip('\n')
        while line != '':
            line = line.split(' ')
            if len(line) == 1:
                indexer = 0
            else:
                indexer = 2
            if 'N' not in line[indexer]:
                tmp = line[indexer][0:(len(line[indexer])-(len(line[indexer])%3))]
                aa1 = translate(tmp)
                tmp = line[indexer][1:(len(line[indexer])-((len(line[indexer])-1)%3))]
                aa2 = translate(tmp)
                tmp = line[indexer][2:(len(line[indexer])-((len(line[indexer])-2)%3))]
                aa3 = translate(tmp)
                write_kmers(aa1,line,outfile,k,AAlist,iteration,i)
                write_kmers(aa2,line,outfile,k,AAlist,iteration,i)
                write_kmers(aa3,line,outfile,k,AAlist,iteration,i)
            i = i+1;
            line = infile.readline()
            line = line.rstrip('\n')
        infile.close()
        for outcount in range(0,200):
            outfile[outcount].close()
    return('Done')

def get_col(args):
    fpath = args[0]
    fnameprefix = args[1]
    fnamepostfix = args[2]
    onamepostfix = args[3]
    col = args[4]
    os.system('awk \'{print $'+str(col)+'}\' '+fpath+'/'+fnameprefix+fnamepostfix+' > '+fpath+'/'+fnameprefix+onamepostfix)
    return('Done')

def sort_file(args):
    fpath = args[0]
    fname = args[1]
    parnum = args[2]
    col = args[3]
    #os.system('sort -T -k '+str(col)+','+str(col)+fpath+' --parallel='+str(parnum)+' -o '+fpath+'/'+fname+'_sorted '+fpath+'/'+fname)
    os.system('sort -k '+str(col)+','+str(col)+' --parallel='+str(parnum)+' -o '+fpath+'/'+fname+'_sorted '+fpath+'/'+fname)                             # make -uo to remove duplicates
    return('Done')

def tumor_specific(args):
    fpath = args[0]
    fname_prefix1 = args[1]
    fname_prefix2 = args[2]
    os.system('comm -23 '+fpath+'/'+fname_prefix1+' '+fpath+'/'+fname_prefix2+' > '+fpath+'/'+fname_prefix1+'_tsa')
    return('Done')

def expr_cutoff(args):
    fname = args[0]
    cut = args[1]
    cmd = 'awk \' $1 >= '+str(cut)+' \' '+fname+' > '+fname+'_hiexp'
    cmd1 = 'awk \'{if($1>'+str(cut)+'){print $2}}\' '+fname+' > '+fname+'_hiexpseqs'
    os.system(cmd)
    os.system(cmd1)
    return('Done')

def get_binding(args):
    fname = args[0]
    netMHCpan_path = args[1]
    HLAalleles = args[2]
    kmer_k = args[3]
    command_string = 'cat';
    for HLAallele in HLAalleles:
        tmp_dir = fname.split('/')
        tmp_dir = '/tmp_'+tmp_dir[len(tmp_dir)-2]+'_'+tmp_dir[len(tmp_dir)-1]
        print(tmp_dir)
        os.system(netMHCpan_path+'/netMHCpan_v2 '+tmp_dir+' -a '+HLAallele+' -l '+str(kmer_k)+' -p '+fname+' > '+fname+'_affinity_'+HLAallele)
        os.system('rm -rf '+netMHCpan_path+tmp_dir)
        command_string = command_string + ' ' + fname + '_affinity_' + HLAallele
    command_string = command_string + ' > ' + fname+'_affinity'
    os.system(command_string)
    for HLAallele in HLAalleles:
        os.system('rm ' + fname+'_affinity_'+HLAallele)
    return('Done')

def is_number(s):
    try:
        int(s)
        return(True)
    except ValueError:
        return(False)

def clean_binding(args):
    fname = args
    print(fname)
    fin = open(fname,'r')
    fout = open(fname+'_clean','w')
    line = fin.readline();
    while line != '':
        line = line.strip('\n')
        if line == '':
            pass
        elif line[0] == '#':
            pass
        elif line[0] =='-':
            line = fin.readline()
            line = fin.readline()
        elif line[0] == ' ' or line[0] == '\t':
            line = line.split()
            if is_number(line[0]):
                fout.write(",".join(line)+'\n')
            else:
                 pass;
        else:
            pass
        line = fin.readline()
    return('Done')

def subset_metadata(args):
    dbfname = args[0] #aminodat
    queryfname = args[1]; #TSA list
    seqname = args[2]; #seq file
    dbinfile = open(dbfname,'r') #open aminodat file
    kmer_k = int(args[3]);
    query_seqs = [line.rstrip('\n') for line in open(queryfname,'r')] #extract TSAs
    kskip = kmer_k+1
    for querycount in range(0,len(query_seqs)):
        os.system('grep '+str(query_seqs[querycount])+' '+dbfname+' > '+queryfname+'_'+str(querycount)+'_temp');
        if querycount == 1:
            os.system('cat '+queryfname+'_'+str(querycount-1)+'_temp '+queryfname+'_'+str(querycount)+'_temp > '+queryfname+'_temp')
        elif querycount > 1:
            os.system('cat '+queryfname+'_'+str(querycount)+'_temp >> '+queryfname+'_temp')
    if len(query_seqs) == 1:
        os.system('cp '+queryfname+'_0_temp '+queryfname+'_temp')
    if len(query_seqs) >= 1:
        os.system('awk \'{print $2,$1}\' '+queryfname+'_temp > '+queryfname+'_temp1')
        os.system('sort -k2,2 '+queryfname+'_temp1 | uniq -s '+str(kskip)+' | awk \'{print $2}\' > '+queryfname+'_temp_line')
        os.system('awk \'FNR==NR{line[$0]=$0; next} FNR in line\' '+queryfname+'_temp_line '+seqname+' > '+queryfname+'_metadata_temp_seq')
        os.system('sort -k2,2 '+queryfname+'_temp1 | uniq -s '+str(kskip)+' | awk \'{print $1}\' > '+queryfname+'_metadata_temp_aa')
        os.system('sort -k2,2 -r '+queryfname+'_temp1 | uniq -s '+str(kskip)+' -d > '+queryfname+'_temp')
        linecount = sp.check_output('wc -l < '+queryfname+'_temp',shell=True);
        linecount = int(linecount.rstrip('\n'))
        while linecount != 0:
            os.system('sort -k2,2 '+queryfname+'_temp | uniq -s '+str(kskip)+' | awk \'{print $2}\' > '+queryfname+'_temp_line')
            os.system('awk \'FNR==NR{line[$0]=$0; next} FNR in line\' '+queryfname+'_temp_line '+seqname+' > '+queryfname+'_metadata_temp_seq1')
            os.system('sort -k2,2 '+queryfname+'_temp | uniq -s '+str(kskip)+' | awk \'{print $1}\' > '+queryfname+'_metadata_temp_aa1')
            os.system('cat '+queryfname+'_metadata_temp_seq1 >> '+queryfname+'_metadata_temp_seq')
            os.system('cat '+queryfname+'_metadata_temp_aa1 >> '+queryfname+'_metadata_temp_aa')
            os.system('sort -k2,2 -r '+queryfname+'_temp | uniq -s '+str(kskip)+' -d > '+queryfname+'_temp1')
            os.system('mv '+queryfname+'_temp1 '+queryfname+'_temp')
            linecount = sp.check_output('wc -l < '+queryfname+'_temp',shell=True);
            linecount = int(linecount.rstrip('\n'))
        os.system('paste '+queryfname+'_metadata_temp_aa '+queryfname+'_metadata_temp_seq | column -s $\'\t\' -t > '+queryfname+'_metadata')        
        os.system('rm '+queryfname+'*temp*')
    print(dbfname+' done')
    return('Done')

def num2letters(num):
    rep = 1
    return([''.join(x) for x in itertools.product('ACDEFGHIKLMNPQRSTVWY', repeat=rep)])

def num2letters_alpha(num):
    rep = 0
    while num > 0:
        rep = rep+1
        num = num/20
    return([''.join(x) for x in itertools.product('abcdefghijklmnopqrstuvwxyz', repeat=rep)])

def gen_input_tuples(file_pos,file_prefix,file_suffix,other_args,num_proc,alphabettype):
    tuple_len = len(other_args)+1
    if alphabettype == 'alpha':
        letters = num2letters_alpha(num_proc)
    elif alphabettype == 'amino':
        letters = num2letters(num_proc)
    tuple_list = tuple()
    for i in range(num_proc):
        for j in range(tuple_len):
            if j==0:
                if j < file_pos:
                    tmp_tuple = other_args[j]
                else:
                    tmp_tuple = file_prefix+letters[i]+file_suffix
            elif j < file_pos:
                if j>1:
                    tmp_tuple = tmp_tuple+(other_args[j],)
                else:
                    tmp_tuple = (tmp_tuple,)+(other_args[j],)
            elif j == file_pos:
                if j>1:
                    tmp_tuple = tmp_tuple+(file_prefix+letters[i]+file_suffix,)
                else:
                    tmp_tuple = (tmp_tuple,)+(file_prefix+letters[i]+file_suffix,)
            else:
                if j>1:
                    tmp_tuple = tmp_tuple+(other_args[j-1],)
                else:
                    tmp_tuple = (tmp_tuple,)+(other_args[j-1],)
        if i==0:
            tuple_list = tmp_tuple
        elif i>1:
            tuple_list = tuple_list+(tmp_tuple,)
        else:
            tuple_list = (tuple_list,)+(tmp_tuple,)
    return(tuple_list)

def dynamic_map(func,input_tuples):
    pool = mp.Pool(processes=len(input_tuples))
    exec_str = 'pool.map('+func+',input_tuples)'
    eval(exec_str)

def multiproc_to_amino1(args):
    patseq_path = args[0]
    fname = args[1]
    min_processes = args[2]
    AAlist = args[3]
    for count in range(0,20):
        sort_file((patseq_path,fname+AAlist[count]+'_aminodat',1,2))
        os.system('rm '+patseq_path+'/'+fname+AAlist[count]+'_aminodat')
        get_col((patseq_path,fname+AAlist[count],'_aminodat_sorted','_amino_sorted',2))
        os.system('uniq -cd '+patseq_path+'/'+fname+AAlist[count]+'_amino_sorted > '+patseq_path+'/'+fname+AAlist[count]+'_amino_sorted_uniq') #tumorpep output file

def multiproc_to_amino2(args):
    patseq_path = args[0]
    fname = args[1]
    peplen = args[2]
    kmer_k = args[3]
    depth_constant = args[4]
    kmer_depth = args[5]
    AAlist = args[6]
    for count in range(0,20):
        #filter out low expressed aminos
        expr_cutoff((patseq_path+'/'+fname+AAlist[count]+'_amino_sorted_uniq',((peplen-kmer_k+1)/(peplen-8+1))*depth_constant*(kmer_depth/2e9)))

def multiproc_to_amino3(args):
    patseq_path = args[0]
    fname_tumor = args[1]
    fname_normal = args[2]
    AAlist = args[3]
    netMHCpan_path = args[4]
    hlaPreds = args[5]
    kmer_k = args[6]
    getmeta = args[7]
    for count in range(0,20):
        #filter out normally expressed aminos
        tumor_specific((patseq_path,fname_tumor+AAlist[count]+'_amino_sorted_uniq_hiexpseqs',fname_normal+fname_tumor[-1]+AAlist[count]+'_amino_sorted_uniq_hiexpseqs'))
        print('tumor_specific '+str(count)+' done')

def multiproc_to_amino4(args):
    patseq_path = args[0]
    fname_tumor = args[1]
    fname_normal = args[2]
    AAlist = args[3]
    netMHCpan_path = args[4]
    hlaPreds = args[5]
    kmer_k = args[6]
    getmeta = args[7]
    for count in range(0,20):
        get_binding((patseq_path+'/'+fname_tumor+AAlist[count]+'_amino_sorted_uniq_hiexpseqs_tsa',netMHCpan_path,hlaPreds,kmer_k))
        print('get_binding '+str(count)+' done')

def multiproc_to_amino5(args):
    patseq_path = args[0]
    fname_tumor = args[1]
    fname_normal = args[2]
    AAlist = args[3]
    netMHCpan_path = args[4]
    hlaPreds = args[5]
    kmer_k = args[6]
    getmeta = args[7]
    for count in range(0,20):
        subset_metadata((patseq_path+'/'+fname_tumor+AAlist[count]+'_aminodat_sorted',patseq_path+'/'+fname_tumor+AAlist[count]+'_amino_sorted_uniq_hiexpseqs_tsa',patseq_path+'/'+fname_tumor[:-1],kmer_k))

def multiproc_merge(args):
    patseq_path = args[0]
    patient = args[1]
    AAlist = args[2]
    post = args[3]
    alphabet = args[4]
    rmstat = args[5]
    for mergecount in range(0,20):
        merge_files((patseq_path,patient+AAlist[mergecount],20,post,alphabet,rmstat))

print('Initialization Done: '+str(round(time.time() - start_time,4))+'s')

def main():
    path = '/fs/project/PCON0005/nas2-data2/colab/laumont/'
    netMHCpan_path = '/users/PAS1203/osu1042/local/netMHCpan-4.0'
    convert_bam_to_fastq = False
    gen_uuid_mapping = False
    kmer_k = int(sys.argv[1])
    getmeta = str(sys.argv[2])
    nfiles = int(sys.argv[3])  #number of processes used, note that minimum is 4
    patient2processn = str(sys.argv[4])
    patient2processt = str(sys.argv[5])
    hlatype = str(sys.argv[6])
    min_processes = [8,4,2]
    if nfiles < 8:
        min_processes[0] = nfiles
    elif nfiles < 4:
        min_processes[0] = nfiles
        min_processes[1] = nfiles
    elif nfiles < 2:
        min_processes[0] = nfiles
        min_processes[1] = nfiles
        min_processes[2] = nfiles
    readlen_tumor = 100
    readlen_normal = 80
    peplen_tumor = round(readlen_tumor/3)
    peplen_normal = round(readlen_normal/3)
    tumor_depth_constant = 10000
    normal_depth_constant = 2
    AAlist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    tumor_fpath = '/fs/project/PCON0005/nas2-data2/colab/laumont/fastq/tumor/' #path+'tumor/'
    normal_fpath = '/fs/project/PCON0005/nas2-data2/colab/laumont/fastq/normal/' #path+'normal/'
    patseq_path = path+str(patient2processt)+'_'+str(kmer_k)+'mer'
    alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    print(tumor_fpath)
    print(normal_fpath)
    print(patseq_path)

    print('MAKING DIRECTORIES AND EXTRACTING SEQUENCES')
    #pool = mp.Pool(processes=min_processes[2])
    if not os.path.isdir(patseq_path):
        print('Creating new output directory')
        os.system('mkdir '+patseq_path)
        #stats = pool.map(bam2seq,((tumor_fpath,patseq_path+'/tumor.txt'),(normal_fpath,patseq_path+'/normal.txt')))
    else:
        print('Output directory already exists. Is this a duplicate process?')
        #print('Deleting old output directory and replacing with new one: '+patseq_path)
        #os.system('rm '+patseq_path+'/*')
        #stats_bam2seq = pool.map(bam2seq,((tumor_fpath,patseq_path+'/tumor.txt'),(normal_fpath,patseq_path+'/normal.txt')))

        print('GETTING FASTQ FOR HLA PREDICTION')
    #if convert_bam_to_fastq:
    #    bam2fastq((path+cancer_type+'/'+line[3],line[4].split('.',1)[0],patseq_path))

    print('GETTING HLA ALLELES')
    if hlatype == 'NULL':
        hlaPreds = fastq2hla((tumor_fpath,patient2processt,'/users/PAS1203/osu1042/local',patient2processt,kmer_k,'paired-end'))
    else:
        hlaPreds = ['HLA-A*11:01', 'HLA-A*23:01', 'HLA-B*35:01','HLA-B*44:03','HLA-C*04:01'] #murine HLA alleles for laumont data

    print('CONVERTING FASTQ FILES TO SEQ FILES')
    os.system('fastx_reverse_complement -i '+normal_fpath+patient2processn+'_R1.fastq -o '+normal_fpath+patient2processn+'_RC_R1.fastq')
    os.system('fastx_reverse_complement -i '+tumor_fpath+patient2processt+'_R1.fastq -o '+tumor_fpath+patient2processt+'_RC_R1.fastq')
    fastq2seq(('paired-end',normal_fpath+patient2processn+'_RC_R1.fastq',normal_fpath+patient2processn+'_R2.fastq',patseq_path+'/'+patient2processn))
    fastq2seq(('paired-end',tumor_fpath+patient2processt+'_RC_R1.fastq',tumor_fpath+patient2processt+'_R2.fastq',patseq_path+'/'+patient2processt))
    #os.system('cp '+patseq_path+'/../'+patient2processn+' '+patseq_path)
    #os.system('cp '+patseq_path+'/../'+patient2processt+' '+patseq_path)

    print('TRANSLATING FILES TO AMINO ACIDS')
    split_file((patseq_path+'/'+patient2processn,20))
    split_file((patseq_path+'/'+patient2processt,20))
    seq_offsett = [''] * nfiles
    seq_offsetn = [''] * nfiles
    seq_offsett[0] = 1
    seq_offsetn[0] = 1
    for offsetcount in range(1,nfiles):
        seq_offsett[offsetcount] = int(sp.check_output('wc -l < '+patseq_path+'/'+patient2processt+alphabet[0]+alphabet[offsetcount-1],shell=True));
        seq_offsetn[offsetcount] = int(sp.check_output('wc -l < '+patseq_path+'/'+patient2processn+alphabet[0]+alphabet[offsetcount-1],shell=True));
        seq_offsett[offsetcount] = seq_offsett[offsetcount]+seq_offsett[offsetcount-1]
        seq_offsetn[offsetcount] = seq_offsetn[offsetcount]+seq_offsetn[offsetcount-1]

    offset_dictt = {alphabet[i]: seq_offsett[i] for i in range(len(seq_offsett))}
    offset_dictn = {alphabet[i]: seq_offsetn[i] for i in range(len(seq_offsetn))}
    input_tuples = gen_input_tuples(0,patseq_path+'/'+patient2processn,'',(kmer_k,offset_dictn),20,'alpha')
    dynamic_map('translate_file',input_tuples)
    input_tuples = gen_input_tuples(1,patient2processn,'',(patseq_path,AAlist,'_aminodat',alphabet,'yes'),20,'amino')
    dynamic_map('multiproc_merge',input_tuples)
    os.system('rm '+patseq_path+'/'+patient2processn+'a*')
    
    input_tuples = gen_input_tuples(0,patseq_path+'/'+patient2processt,'',(kmer_k,offset_dictt),20,'alpha')
    dynamic_map('translate_file',input_tuples)
    input_tuples = gen_input_tuples(1,patient2processt,'',(patseq_path,AAlist,'_aminodat',alphabet,'yes'),20,'amino')
    dynamic_map('multiproc_merge',input_tuples)
    os.system('rm '+patseq_path+'/'+patient2processt+'a*')

    print('SORTING AMINODAT FILES, RETRIEVING AMINOS, & GETTING UNIQUE AMINOS')
    for count in range(0,20):
        for count2 in range(0,20):
            sort_file((patseq_path,patient2processt+AAlist[count]+AAlist[count2]+'_aminodat',20,2))
            os.system('rm '+patseq_path+'/'+patient2processt+AAlist[count]+AAlist[count2]+'_aminodat')
            get_col((patseq_path,patient2processt+AAlist[count]+AAlist[count2],'_aminodat_sorted','_amino_sorted',2))
            os.system('uniq -cd '+patseq_path+'/'+patient2processt+AAlist[count]+AAlist[count2]+'_amino_sorted > '+patseq_path+'/'+patient2processt+AAlist[count]+AAlist[count2]+'_amino_sorted_uniq')
            os.system('rm '+patseq_path+'/'+patient2processt+AAlist[count]+AAlist[count2]+'_amino_sorted')

    for count in range(0,20):
        for count2 in range(0,20):
            sort_file((patseq_path,patient2processn+AAlist[count]+AAlist[count2]+'_aminodat',20,2))
            os.system('rm '+patseq_path+'/'+patient2processn+AAlist[count]+AAlist[count2]+'_aminodat')
            get_col((patseq_path,patient2processn+AAlist[count]+AAlist[count2],'_aminodat_sorted','_amino_sorted',2))
            os.system('uniq -cd '+patseq_path+'/'+patient2processn+AAlist[count]+AAlist[count2]+'_amino_sorted > '+patseq_path+'/'+patient2processn+AAlist[count]+AAlist[count2]+'_amino_sorted_uniq') #tumorpep output file
            os.system('rm '+patseq_path+'/'+patient2processn+AAlist[count]+AAlist[count2]+'_amino_sorted')
            os.system('rm '+patseq_path+'/'+patient2processn+AAlist[count]+AAlist[count2]+'_aminodat_sorted')                                      
    
        print('MERGING AMINO FILES')
    pool = mp.Pool(processes=min_processes[2])
    stats_merge_file = pool.map(merge_files, ((patseq_path,patient2processn,400,'_amino_sorted_uniq',AAlist,'no'),(patseq_path,patient2processt,400,'_amino_sorted_uniq',AAlist,'no')))

    print('REMOVING SEQUENCES BELOW EXPRESSION CUTOFF, FILTER OUT NORMAL AMINOS, & PREDICT BINDING AFFINITY')
    tumor_kmer_depth = sum(1 for line in open(patseq_path+'/'+patient2processt+'_amino_sorted_uniq'))
    normal_kmer_depth = sum(1 for line in open(patseq_path+'/'+patient2processn+'_amino_sorted_uniq'))
    print('tumor cutoff'+str(tumor_depth_constant*(tumor_kmer_depth/2e9)))
    print('normal cutoff' +str(normal_depth_constant*(normal_kmer_depth/2e9)))
    #os.system('rm '+patseq_path+'/'+patient2processt+'_amino_sorted_uniq')
    os.system('rm '+patseq_path+'/'+patient2processn+'_amino_sorted_uniq')
    
    input_tuples = gen_input_tuples(1,patient2processt,'',(patseq_path,peplen_tumor,kmer_k,tumor_depth_constant,tumor_kmer_depth,AAlist),20,'amino')
    dynamic_map('multiproc_to_amino2',input_tuples)

    input_tuples = gen_input_tuples(1,patient2processn,'',(patseq_path,peplen_normal,kmer_k,normal_depth_constant,normal_kmer_depth,AAlist),20,'amino')
    dynamic_map('multiproc_to_amino2',input_tuples)

    print('MERGING AMINO FILES')
    pool = mp.Pool(processes=min_processes[2])
    stats_merge_file = pool.map(merge_files, ((patseq_path,patient2processn,400,'_amino_sorted_uniq_hiexp',AAlist,'no'),(patseq_path,patient2processt,400,'_amino_sorted_uniq_hiexp',AAlist,'no')))
    
    for rmcount in range(0,20):
        for rmcount2 in range(0,20):
            os.system('rm '+patseq_path+'/'+patient2processn+AAlist[rmcount]+AAlist[rmcount2]+'_amino_sorted_uniq')
            os.system('rm '+patseq_path+'/'+patient2processt+AAlist[rmcount]+AAlist[rmcount2]+'_amino_sorted_uniq')
            os.system('cp '+patseq_path+'/'+patient2processn+AAlist[rmcount]+AAlist[rmcount2]+'_amino_sorted_uniq_hiexpseqs ..')
            os.system('rm '+patseq_path+'/'+patient2processn+AAlist[rmcount]+AAlist[rmcount2]+'_amino_sorted_uniq_hiexp')
            os.system('rm '+patseq_path+'/'+patient2processt+AAlist[rmcount]+AAlist[rmcount2]+'_amino_sorted_uniq_hiexp')
    
    
    input_tuples = gen_input_tuples(1,patient2processt,'',(patseq_path,patient2processn,AAlist,netMHCpan_path,hlaPreds,kmer_k,getmeta,patient2processt),20,'amino')
    dynamic_map('multiproc_to_amino3',input_tuples)

    input_tuples = gen_input_tuples(1,patient2processt,'',(patseq_path,patient2processn,AAlist,netMHCpan_path,hlaPreds,kmer_k,getmeta,patient2processt),20,'amino')
    dynamic_map('multiproc_to_amino4',input_tuples)
    if getmeta == 'getmeta':
        input_tuples = gen_input_tuples(1,patient2processt,'',(patseq_path,patient2processn,AAlist,netMHCpan_path,hlaPreds,kmer_k,getmeta,patient2processt),20,'amino')
        dynamic_map('multiproc_to_amino5',input_tuples)
    
    #Merge and clean the netMHCpan output
    merge_files((patseq_path,patient2processt,400,'_amino_sorted_uniq_hiexpseqs_tsa_affinity',AAlist,'yes'))
    merge_files((patseq_path,patient2processt,400,'_amino_sorted_uniq_hiexpseqs_tsa_metadata',AAlist,'yes'))
    clean_binding(patseq_path+'/'+patient2processt+'_amino_sorted_uniq_hiexpseqs_tsa_affinity')

    print('CLEANING DIRECTORY AND WRITING FILES')
    if getmeta=='getmeta':
        os.system('mv '+patseq_path+'/'+patient2processt+'_amino_sorted_uniq_hiexpseqs_tsa_metadata '+path+'/'+str(patient2processt)+'_'+str(kmer_k)+'mer_neoantmetadata.txt')
    os.system('mv '+patseq_path+'/'+patient2processt+'_amino_sorted_uniq_hiexpseqs_tsa_affinity_clean '+path+'/'+str(patient2processt)+'_'+str(kmer_k)+'mer_neoantigens.txt')
    os.system('mv '+patseq_path+'/'+patient2processt+'_amino_sorted_uniq_hiexp '+path+'/'+str(patient2processt)+'_'+str(kmer_k)+'mer_tumorpep_hiexp.txt')
    os.system('mv '+patseq_path+'/'+patient2processt+'_amino_sorted_uniq '+path+'/'+str(patient2processt)+'_'+str(kmer_k)+'mer_tumorpep.txt')
    os.system('mv '+patseq_path+'/'+patient2processn+'_amino_sorted_uniq_hiexp '+path+'/'+str(patient2processt)+'_'+str(kmer_k)+'mer_normalpep.txt')
    meta_out = open(path+'/'+str(patient2processt)+'_'+str(kmer_k)+'mer_patmetadata.txt','w')
    meta_out.write(tumor_fpath+'\n')
    meta_out.write(normal_fpath+'\n')
    meta_out.write(patseq_path+'\n')
    meta_out.write(str(hlaPreds)+'\n')
    meta_out.close()
    os.system('rm -rf '+patseq_path)

if __name__ == "__main__":
    main()
    print('Finished: '+str(round(time.time()-start_time,6))+'s')
