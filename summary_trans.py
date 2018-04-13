#!/usr/bin/env python
#coding:utf-8
### 比对步骤
## cd /lustre/work/yongdeng/Project && python ~/align_clean.py -i 160761_test/QC/ -g 160761_test/genome/ -mf 160761_test/align_clean.mf -o 160761_test/align
import argparse,os,sys,re
import RNA
from collections import defaultdict

parser = argparse.ArgumentParser(description="This Script is used to align cleandata to reference genome by bowtie2 program.")
parser.add_argument("-i","--input_gtf",type=str,help="the gtf file for summary",required=True,nargs="*")
parser.add_argument("-o","--output_stat",type=str,help="the output summary file of all input gtf",required=True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
args = parser.parse_args()

def summary_gtf(gtfile):
    gene_id=set()
    transcriptDict=defaultdict(int)
    with open(gtfile) as gtf:
        for line in gtf:
            if line.startswith("#") or not line.split("\t")[2]=="exon":
                continue
            line_list = line.strip().split("\t")
            gene_id.add(re.search('gene_id "(.+?)";',line_list[-1]).group(1))
            transcript_key = re.search('transcript_id "(.+?)";',line_list[-1]).group(1)
            transcriptDict[transcript_key] += abs(int(line_list[3])-int(line_list[4])) + 1
    return gene_id,transcriptDict
    
    
def N50(seq_list):
    seq_list_sorted = sorted(seq_list,reverse=True)
    base_sum = sum(seq_list)
    N50_pos = base_sum / 2.0 
    ValueSum,n50 = 0,0 
    for value in seq_list_sorted:
        ValueSum += value
        if ValueSum >= N50_pos:
            n50 = value
            break
    return n50
    
def get_len_dis(start,end,seq):
    a = 0
    for i in seq:
        if i >= start and i <= end:
            a+=1
    return a    

def main():
    fo = open(args.output_stat,"w")
    fo.write("Sample\tGene Number\tTranscript Number\tExon Total length(bp)\tAverage Transcript length(bp)\tMax Transcript length(bp)\tMin Transcript length(bp)\tN50 length(bp,without intron)" + "\n")
    gtflist = args.input_gtf
    for gtfile in gtflist:
        gene_id,transcriptDict = summary_gtf(gtfile)
        gene_num = len(gene_id)
        transcript_num = len(transcriptDict)
        total_len = sum(transcriptDict.values())
        ave_len = "%.f"%(float(total_len)/transcript_num)
        max_T = max(transcriptDict.values())
        min_T = min(transcriptDict.values())
        N_50 = N50(transcriptDict.values())
        fo.write(os.path.basename(os.path.dirname(os.path.abspath(gtfile))) + "\t" + '\t'.join(map(str,[gene_num,transcript_num,total_len,ave_len,max_T,min_T,N_50])) + "\n")
        with open(gtfile+".length_stat","w") as lenstat:
            lenstat.write('Gene_Number\tTrans_Number\tTotal_Length(bp)\tAve_len\tMax\tMin\tN50\n')
            lenstat.write('\t'.join(map(str,[gene_num,transcript_num,total_len,ave_len,max_T,min_T,N_50])) + "\n")
        n = max_T/100
        with open(gtfile+".transcript.length_distribution","w") as lendis:
            lendis.write('Length\tNumber\n')
            for i in range(1,n+1):
                start = i*100
                if start >= 3000:
                    break
                end = (i+1)*100 - 1
                lendis.write(str(start)+ "-" + str(end) + "\t" + str(get_len_dis(start,end,transcriptDict.values())) + "\n")
            lendis.write(">=3000" + "\t" + str(get_len_dis(start,max_T,transcriptDict.values())) + "\n")
    fo.close()
            
if __name__ == "__main__":
    main()
    
