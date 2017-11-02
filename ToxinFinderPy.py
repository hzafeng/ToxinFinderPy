# -*- coding: utf-8 -*-
import os
import argparse
import subprocess
what_i_do = "a script to find the toxin gene in BT's genome or metagenome"
parser = argparse.ArgumentParser(description=what_i_do)
parser.add_argument('-i', dest='input_files', type=str, nargs='+',
                required=True, help='input genome file in fasta(required)', default=None)
parser.add_argument('-db', dest='local_db', type=str, nargs=1,
                required=True, help='local database_files for blast in fasta(required)', default=None)                
parser.add_argument('-o', dest='output_files', type=str, nargs=1,
                required=True, help='the output file (required)', default=None)
parser.add_argument('-t', dest='threads_num', type=str, nargs=1,
                required=False, help='the num of threads used for blast', default=1)

def check_arguments(arguments):
   for file in arguments["input_files"]:
       try:
           open(file)
           close(file)
       except IOError:
           print "Could not open " + file
           sys.exit(1)
args = vars(parser.parse_args())

print "**********************************************************************"

order_makedb = 'makeblastdb'+' -in '+args['local_db'][0]+' -out '+'known_gene '+'-dbtype '+'prot'
subprocess.call(order_makedb,shell = True)
aa_files=[];blast_outs=[]
for f in args["input_files"]:
   
   genome_name = f.split('/')[-1].split('.fasta')[0]

   log = genome_name+'.log'
   protein_file = genome_name+'_aa.fasta'
   aa_files.append(protein_file)
   prodigal_order =  'prodigal'+' -i '+f+' -c'+' -m'+' -g'+' 11'+' -f'+' gbk'+' -q'+' -o'+' ' +log+' -a '+' '+protein_file

   subprocess.call(prodigal_order,shell = True)
   out_fn = genome_name +'.blast.out'    
   order1 = 'blastp'+' -query '+protein_file+' -db'+' known_gene'+' -evalue'+' 10'+' -num_alignments'+' 5'+' -outfmt'+' 6'+' -num_threads '+args['threads_num'][0]+' -out '+out_fn

   blast_outs.append(out_fn)
   subprocess.call(order1,shell = True)
dict_all =''
for fi in aa_files:
   of = open(fi,'r')
   dict_all += of.read()
   of.close()
dict_ = dict_all.split('>')
dict1_=dict_[1:]
dict1={}
for contig in dict1_:
    aaa=contig.split('\n',1)
    gene_id = aaa[0].split()[0]
    dict1[gene_id] = aaa[1].replace('\n','')

def excert():
    f4 = open(args['output_files'][0]+'_4.csv','w')
    f3 = open(args['output_files'][0]+'_3.csv','w')
    f2 = open(args['output_files'][0]+'_2.csv','w')
    f1 = open(args['output_files'][0]+'_1.csv','w')
    str1=''
    str2=''
    str3=''
    str4=''
    str_hit = ''
    str_excert = ''
    for f in blast_outs:
         files = open(f,'r')
         str_hit += files.read().replace('\n\n','')
         files.close()
    each_line_ = str_hit.split('\n')
    for each_line in each_line_:
        if each_line.split():
            each_ = each_line.split('\t')
            if float(each_[10]) < 0.00001 :
                leng_query = len(dict1[each_[0]])
                leng_target = each_[1].split('_len')[-1]
                bi = format(float(leng_query)/float(leng_target),'.2f')
                cov = format(float(each_line[3])/float(length_query),'.2f')
                if float(bi) > 0.6 and float(cov) > 0.5:
                    if float(each_[2]) > 95 :
                        str4 += '>'
                        str4 += each_[0]
                        str4 += '------------------------'
                        str4 += each_[1]
                        str4 += '\n'
                        str4 += dict1[each_[0]]+'\n'
                    if float(each_[2]) < 45 :
                        str1 += '>'
                        str1 += each_[0]
                        str1 += '------------------------'
                        str1 += each_[1]
                        str1 += '\n'
                        str1 += dict1[each_[0]]+'\n'
                    if float(each_[2]) >= 45 and float(each_[2]) <= 78:
                        str2 += '>'
                        str2 += each_[0]
                        str2 += '------------------------'
                        str2 += each_[1]
                        str2 += '\n'
                        str2 += dict1[each_[0]]+'\n'               
                    if float(each_[2]) >= 78 and float(each_[2]) <= 95:  
                        str3 += '>'
                        str3 += each_[0]
                        str3 += '------------------------'
                        str3 += each_[1]
                        str3 += '\n'
                        str3 += dict1[each_[0]]+'\n'  
    f1.write(str1)
    f2.write(str2)
    f3.write(str3)
    f4.write(str4)

excert()

print 'done^_^'

