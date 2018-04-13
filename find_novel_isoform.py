#!/usr/bin/env python
#coding:utf-8

import os,argparse,re

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is get all novel trans and count the number of know and novel")
    parser.add_argument("-gtf","--gtf",type=str,help="The input gtf file",required = True, nargs="+")
    parser.add_argument("-ref","--ref_gtf",type=str,help="The ref_gtf file",required = True)
    parser.add_argument("-o","--outputdir",type=str,help="The output dir of novel isoform gtf file")
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()

def get_ref_trans_id(refgtf):
    reftransset = set()
    refgeneset = set()
    with open(refgtf) as fi:
        for line in fi:
            if not line or line.startswith('#'):
                continue
            if '#' in line:                        
                line = line.split('#')[0].strip()
            if re.search('transcript_id ".+?";',line):
                reftransset.add(re.search('transcript_id "(.+?)";',line).group(1))
            if re.search('gene_id ".+?";',line):
                refgeneset.add(re.search('gene_id "(.+?)";',line).group(1))
    return refgeneset,reftransset
            
def main():
    args = parseArg()
    refgene,reftrans = get_ref_trans_id(args.ref_gtf)
    if not os.path.isdir(args.outputdir):os.makedirs(args.outputdir)
    for gtf in args.gtf:
        total_trans = 0
        total_gene = set()
        novel_trans = 0
        know_gene = set()
        transcripts = []
        with open(gtf) as fi,open(os.path.join(args.outputdir,"novel_transcripts.gtf"),"w") as fo:
            for line in fi:
                if not line or line.startswith('#'): continue
                if re.search('gene_id ".+?";',line):
                    geneid = re.search('gene_id "(.+?)";',line).group(1)
                    total_gene.add(geneid)
                    if geneid in refgene:
                        know_gene.add(geneid)
                if "\ttranscript\t" in line:
                    total_trans += 1
                    if not len(transcripts):
                        transcripts = [line,]
                        continue 
                    trans = transcripts[0]
                    if re.search('transcript_id "(.+?)";',trans).group(1) not in reftrans:
                        fo.writelines(transcripts)
                        novel_trans += 1
                    transcripts = [line,]
                else:
                    transcripts.append(line)
            trans = transcripts[0]
            total_trans += 1
            if re.search('transcript_id "(.+?)";',trans).group(1) not in reftrans:
                fo.writelines(transcripts)
                novel_trans += 1
        with open(os.path.join(args.outputdir,"transcripts.gtf.known_and_novel_stat"),"w") as fo:
            fo.write('Total Assembled Transcripts\t' + str(total_trans) + "\n")
            fo.write('Total Assembled Loci\t' + str(len(total_gene)) + "\n")
            fo.write('In Known Transcripts\t' + str(total_trans-novel_trans) + "\n")
            fo.write('Corresponding Known Loci\t' + str(len(know_gene)) + "\n")
            fo.write('In Novel Transcripts\t' + str(novel_trans) + "\n")
            fo.write('Corresponding Novel Loci\t' +str(len(total_gene) - len(know_gene)) + "\n")
            
                
                
if __name__ == "__main__":
    main()
             
        
