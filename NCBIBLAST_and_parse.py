from __future__ import division
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
import os
import argparse
current_path = os.path.abspath(".")
files = os.listdir('.')
for each_file in files:
    if '.table' in each_file:
        work_name = each_file.split('.')[0]
        input_file = open(os.path.join(current_path,each_file),'r')
        input_data = input_file.read().split('\n')[1:-1]
        all_id_num = len(input_data)
        input_file.close()
        out_file_name = work_name+'.fasta'
        out_fasta = open(os.path.join(current_path,out_file_name),'w')
        for each_line in input_data:
            each__ = each_line.split('\t')
            out_fasta.write('>'+each__[0]+'\n'+each__[-2]+'\n')
        out_fasta.close()
    else:
        pass
#_____________________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________________
count = 1
for seq_record in SeqIO.parse(os.path.join(current_path,out_file_name), "fasta"):
         
    print 'Running',count,',totle',all_id_num
    count += 1
    blast_handle_nr = NCBIWWW.qblast("blastp", "nr",seq_record.seq,alignments=500,expect=0.00001)
    NCBI_blast_result_file_nr = seq_record.id+'_nr.xml'
    blast_result_nr = open(os.path.join(current_path,NCBI_blast_result_file_nr),'w')                    
    blast_result_nr.write(blast_handle_nr.read())
    blast_result_nr.close()
    result_handle_nr = open(os.path.join(current_path,NCBI_blast_result_file_nr),'r')
    blast_record_nr = NCBIXML.read(result_handle_nr)
    list_have_nr = []
    title_nr = []
    title_all_nr= []
    for i in range(0,len(blast_record_nr.alignments)):
        title_all_nr.append(blast_record_nr.descriptions[i].title)
        if 'ICP' in blast_record_nr.descriptions[i].title or 'Insecticidal' in blast_record_nr.descriptions[i].title or 'Crystal Protein' in blast_record_nr.descriptions[i].title or 'toxin'  in blast_record_nr.descriptions[i].title or 'Cry' in blast_record_nr.descriptions[i].title or 'Cyt' in blast_record_nr.descriptions[i].title or 'Vip' in blast_record_nr.descriptions[i].title:
            list_have_nr.append(i)
            title_nr.append(blast_record_nr.descriptions[i].title)
    a=1
    query_nr = []
    match_nr = []
    expect_nr = []
    identity_nr = []
    cov_nr = []
    
    all_txt_filename_nr = seq_record.id+'_nr_all.txt'
    all_txt_nr = open(os.path.join(current_path,all_txt_filename_nr),'w')
    title_count = 0
    for alignment in blast_record_nr.alignments:
        
        for hsp in alignment.hsps:
            all_txt_nr.write("*************Alignments****************"+'\n'+str(title_all_nr[title_count])+'\n'+'identity:'+str(int(hsp.identities)/int(hsp.align_length)*100)+'\n'+'coverage:'+str(int(hsp.align_length)/(int(hsp.query_end)-int(hsp.query_start)+1)*100)+'\n'+'expect:'+str(hsp.expect)+'query:'+str(hsp.query)+'\n'+'match:'+str(hsp.match)+'\n')

        title_count +=1 

    all_txt_nr.close()


    if len(list_have_nr) == 0:
        pass
    else:
        for alignment in blast_record_nr.alignments:
            
            if a in list_have_nr:
                for hsp in alignment.hsps:
                    identity_nr.append(int(hsp.identities)/int(hsp.align_length)*100)
                    cov_nr.append(int(hsp.align_length)/(int(hsp.query_end)-int(hsp.query_start)+1)*100)
                    query_nr.append(hsp.query)
                    match_nr.append(hsp.match)
                    expect_nr.append(hsp.expect) 
            a += 1
        txt_file_name_nr = seq_record.id+'_nr_toxin_hit.txt'
        csv_file_name_nr = seq_record.id+'_nr_toxin_hit.csv'
        txt_file_nr = open(os.path.join(current_path,txt_file_name_nr),'w')
        csv_file_nr = open(os.path.join(current_path,csv_file_name_nr),'w')
        csv_file_nr.write('id'+','+'identity'+','+'coverage'+','+'Except'+','+'title'+','+'query'+','+'hit'+'\n')
        for q in range(0,len(list_have_nr)):
            p = q+1
            txt_file_nr.write("*************Alignments****************"+'\n'+title_nr[q]+'\n'+'identity(%):'+str(identity_nr[q])+'\n'+'coverage(%):'+str(cov_nr[q])+'\n'+'expect:'+str(expect_nr[q])+'\n'+'query:'+str(query_nr[q])+'\n'+'match:'+str(match_nr[q])+'\n')
            csv_file_nr.write(str(p)+','+str(identity_nr[q])+','+str(cov_nr[q])+','+str(expect_nr[q])+','+str(title_nr[q])+','+str(query_nr[q])+','+str(match_nr[q])+'\n')
        txt_file_nr.close()
        csv_file_nr.close()
        
    blast_handle_swi = NCBIWWW.qblast("blastp", "swissprot",seq_record.seq,alignments=500,expect=0.00001)
    NCBI_blast_result_file_swi = seq_record.id+'_swi.xml'
    blast_result_swi = open(os.path.join(current_path,NCBI_blast_result_file_swi),'w')                    
    blast_result_swi.write(blast_handle_swi.read())
    blast_result_swi.close()
    result_handle_swi = open(os.path.join(current_path,NCBI_blast_result_file_swi),'r')
    blast_record_swi = NCBIXML.read(result_handle_swi)
    list_have_swi = []
    title_swi = []
    title_all_swi= []
    for i in range(0,len(blast_record_swi.alignments)):
        title_all_swi.append(blast_record_swi.descriptions[i].title)
        if 'ICP' in blast_record_swi.descriptions[i].title or 'Insecticidal' in blast_record_swi.descriptions[i].title or 'Crystal Protein' in blast_record_swi.descriptions[i].title or 'toxin'  in blast_record_swi.descriptions[i].title or 'Cry' in blast_record_swi.descriptions[i].title or 'Cyt' in blast_record_swi.descriptions[i].title or 'Vip' in blast_record_swi.descriptions[i].title:
            list_have_swi.append(i)
            title_swi.append(blast_record_swi.descriptions[i].title)
    b=1
    query_swi = []
    match_swi = []
    expect_swi = []
    identity_swi = []
    cov_swi = []
    
    all_txt_filename_swi = seq_record.id+'_swi_all.txt'
    all_txt_swi = open(os.path.join(current_path,all_txt_filename_swi),'w')
    title_count = 0
    for alignment in blast_record_swi.alignments:
        
        for hsp in alignment.hsps:
            all_txt_swi.write("*************Alignments****************"+'\n'+str(title_all_swi[title_count])+'\n'+'identity:'+str(int(hsp.identities)/int(hsp.align_length)*100)+'\n'+'coverage:'+str(int(hsp.align_length)/(int(hsp.query_end)-int(hsp.query_start)+1)*100)+'\n'+'expect:'+str(hsp.expect)+'query:'+str(hsp.query)+'\n'+'match:'+str(hsp.match)+'\n')

        title_count +=1 

    all_txt_swi.close()


    if len(list_have_swi) == 0:
        pass
    else:
        for alignment in blast_record_swi.alignments:
            
            if b in list_have_swi:
                for hsp in alignment.hsps:
                    identity_swi.append(int(hsp.identities)/int(hsp.align_length)*100)
                    cov_swi.append(int(hsp.align_length)/(int(hsp.query_end)-int(hsp.query_start)+1)*100)
                    query_swi.append(hsp.query)
                    match_swi.append(hsp.match)
                    expect_swi.append(hsp.expect) 
            b += 1
        txt_file_name_swi = seq_record.id+'_swi_toxin_hit.txt'
        csv_file_name_swi = seq_record.id+'_swi_toxin_hit.csv'
        txt_file_swi = open(os.path.join(current_path,txt_file_name_swi),'w')
        csv_file_swi = open(os.path.join(current_path,csv_file_name_swi),'w')
        csv_file_swi.write('id'+','+'identity'+','+'coverage'+','+'Except'+','+'title'+','+'query'+','+'hit'+'\n')
        for q in range(0,len(list_have_swi)):
            p = q+1
            txt_file_swi.write("*************Alignments****************"+'\n'+title_swi[q]+'\n'+'identity(%):'+str(identity_swi[q])+'\n'+'coverage(%):'+str(cov_swi[q])+'\n'+'expect:'+str(expect_swi[q])+'\n'+'query:'+str(query_swi[q])+'\n'+'match:'+str(match_swi[q])+'\n')
            csv_file_swi.write(str(p)+','+str(identity_swi[q])+','+str(cov_swi[q])+','+str(expect_swi[q])+','+str(title_swi[q])+','+str(query_swi[q])+','+str(match_swi[q])+'\n')
        txt_file_swi.close()
        csv_file_swi.close()

    print seq_record.id ,'done'


