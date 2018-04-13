#!/usr/bin/env python
#coding:utf-8

import argparse,re,os
from math import log

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is get all gene anno table from total gtf file and anno file")
    parser.add_argument("-gtf","--gtf",type=str,help="The input gtf file.",required = True)
    parser.add_argument("-id","--idlist",type=str,help="The transcripts id list file",required = True)
    parser.add_argument("-o","--output",type=str,help="The novel trans output dir")
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def get_trans_dict(gtf):
    transdict = {}
    with open(gtf) as fi:
        for line in fi:
            if "\texon\t" in line:
                transid = re.search('transcript_id "(.+?)";',line).group(1)
                exonpos = map(int,line.split("\t")[3:5])
                transdict.setdefault(transid,[]).append(exonpos)
    return transdict
    
def main():
    args = parseArg()
    transdict = get_trans_dict(args.gtf)
    with open(args.idlist) as fi,open(args.output,"w") as fo:
        h = fi.next()
        fo.write("Type\tTranscript_exon_num\tlog2Transcript_exon_num\tTranscript_length\tlog2Transcript_length\n")
        for line in fi:
            id = line.strip().split("\t")[0]
            if id in transdict:
                Transcript_exon_num = str(len(transdict[id]))
                log2Transcript_exon_num = str(log(int(Transcript_exon_num),2))
                trans_len = str(sum([abs(i[0]-i[1]+1) for i in transdict[id]]))
                log2trans_len = str(log(int(trans_len),2))
                fo.write("\t".join(["Transcript",Transcript_exon_num,log2Transcript_exon_num,trans_len,log2trans_len]) + "\n")
                
if __name__ == "__main__":
    main()
                
                
