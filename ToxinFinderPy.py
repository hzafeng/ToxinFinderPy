# -*- coding: utf-8 -*-
#!/usr/bin/env python
import os
import argparse
import subprocess
what_i_do = "a script to find the toxin gene in BT's genome or metagenome"
parser = argparse.ArgumentParser(description=what_i_do)
parser.add_argument('-i', dest='input_files', type=str, nargs=1,
                required=True, help='input genome file in fasta(required)', default=None)
parser.add_argument('-db', dest='local_db', type=str, nargs=1,
                required=True, help='local database_files for blast in fasta(required)', default=None)                
parser.add_argument('-o', dest='output_files', type=str, nargs=1,
                required=True, help='the output file (required)', default=None)

def check_arguments(arguments):
   for file in arguments["input_files"]:
       try:
           open(file)
           close(file)
       except IOError:
           print "Could not open " + file
           sys.exit(1)
args = vars(parser.parse_args())
print args

print "running prodigal ..."
for f in args["input_files"]:
    log = genome_name+'.log'
    protein_file = genome_name+'.fasta'
    prodigal_order =  'prodigal'+' -i '+f+' -c'+' -m'+' -g'+' 11'+' -f'+' gbk'+' -q'+' -o'+' ' +os.path.join(path_1,log)+' -a '+' '+os.path.join(path_2,protein_file)
    subprocess.call(prodigal_order,shell = True)
    print "prodigal done !"
              
dict_all =''

print "run makeblastdb ..."
order_makedb = 'makeblastdb'+' -in '+ os.path.join(database_path,database_filename)+' -out '+'known_gene '+'-dbtype '+'prot'
subprocess.call(order_makedb,shell = True)
print "makeblastdb done !"
all_file_1 = os.listdir(".")    #刷新当前工作目录下的文件
print "run blast ..."
for each_genome in all_file_1:
    if '.fasta' in each_genome:
        data_file = open(path+"/"+each_genome,'r')
        dict_all += data_file.read()
        data_file.close() 
        name = each_genome.split('.')[0]
        out_fn = name +'.out'    
        order1 = 'blastp'+' -query '+os.path_2.join(path_2,each_genome)+' -db'+' known_gene'+' -evalue'+' 5'+' -num_alignments'+' 1'+' -outfmt'+' 6'+' -num_threads'+' 8'+' -out '+os.path.join(out_path,out_fn)
        subprocess.call(order1,shell = True)
print 'blast done!'
#-----------excert blast by identity-----------
dict_ = dict_all.split('>')
dict1_=dict_[1:];dict1={}
for contig in dict1_:
    aaa=contig.split('\n',1)
    gene_id = aaa[0].split()[0]
    dict1[gene_id] = aaa[1]
#-------------建立每个基因长度和id对应的字典------------
def excert(paths):
    f4 = open(path+'/'+'4.csv','w')
    f3 = open(path+'/'+'3.csv','w')
    f2 = open(path+'/'+'2.csv','w')
    f1 = open(path+'/'+'1.csv','w')
    str1=''
    str2=''
    str3=''
    str4=''
    str_hit = ''
    str_excert = ''
    for root,dirs,files in os.walk(paths):
        for fn in files:
            files = open(paths+"/"+fn,'r')
            str_hit += files.read()

    each_line_ = str_hit.split('\n')
    for each_line in each_line_:
        if each_line.split():
            each_ = each_line.split('\t')
            if float(each_[10]) < 0.00001 :
                num_n = len(dict1[each_[0]])/60
                leng_query = len(dict1[each_[0]]) - num_n
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
                        str4 += dict1[each_[0]]
                    if float(each_[2]) < 45 :
                        str1 += '>'
                        str1 += each_[0]
                        str1 += '------------------------'
                        str1 += each_[1]
                        str1 += '\n'
                        str1 += dict1[each_[0]]
                    if float(each_[2]) >= 45 and float(each_[2]) <= 78:
                        str2 += '>'
                        str2 += each_[0]
                        str2 += '------------------------'
                        str2 += each_[1]
                        str2 += '\n'
                        str2 += dict1[each_[0]]                
                    if float(each_[2]) >= 78 and float(each_[2]) <= 95:  
                        str3 += '>'
                        str3 += each_[0]
                        str3 += '------------------------'
                        str3 += each_[1]
                        str3 += '\n'
                        str3 += dict1[each_[0]]       
    f1.write(str1)
    f2.write(str2)
    f3.write(str3)
    f4.write(str4)

excert(out_path)

print 'done^_^'
