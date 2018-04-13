#!/usr/bin/env python
#coding:utf-8

import argparse,os,sys,re,fileinput
from collections import defaultdict

parser = argparse.ArgumentParser(description="This Script is used to change th gene_id name in *.gtf generated by stringTie,the input gtf file will be covered by output gtf file and be renamed as *.gtf.orl. replaced id will be output in stderr")
parser.add_argument("-i","--input_gtf",type=str,help="the input transcripts.gtf file",required = True,nargs='+')
parser.add_argument("-g","--ref_gtf",type=str,help="The ref gtf file",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
args = parser.parse_args()

def get_replace_gene_dict(gtf,refgene_info):
    replace_gene_dict = {}   ## trackgene:refgene
    gene_dict = {}    ## trackgene:refgene
    gene_info = defaultdict(dict) ## trackgene:trackgene_pos
    with open(gtf) as fi:
        for line in fi:
            if line.strip() and not line.startswith("#"):
                line = line.strip().split("\t")
                description = line[-1]
                geneid = re.search('gene_id "(.+?)";',description).group(1)
                gene_info[geneid].setdefault("pos",[]).extend([int(line[3]),int(line[4])])
                gene_info[geneid].setdefault("chr",set()).add(line[0])
                gene_info[geneid].setdefault("strand",set()).add(line[6])
                #if re.search('gene_id ".+?";',description):
                if re.search('ref_gene_id ".+?";',description):
                    gene_dict.setdefault(geneid,set()).add(re.search('ref_gene_id "(.+?)";',description).group(1))
    for k,v in gene_dict.items():
        if len(v) == 1:
            v = list(v)[0]
            if k == v:continue
            trackinfo = gene_info[k]
            chrinfo = list(trackinfo["chr"])[0]
            posinfo = [min(trackinfo["pos"]),max(trackinfo["pos"])]
            strandinfo = list(trackinfo["strand"])[0]
            if [chrinfo,posinfo,strandinfo] == [refgene_info[v]["chr"],refgene_info[v]["pos"],refgene_info[v]["strand"]]:
                replace_gene_dict[k] = v
    return replace_gene_dict
   
def get_refgene_info(refgtf):
    refgene_info = defaultdict(dict)
    with open(refgtf) as fi:
        for line in fi:
            if "\tgene\t" in line:
                gname = re.search('gene_id "(.+?)";',line).group(1)
                line = line.split("\t")
                refgene_info[gname]["pos"] = [int(line[3]),int(line[4])]
                refgene_info[gname]["chr"] = line[0]
                refgene_info[gname]["strand"] = line[6]
    return refgene_info

def main():
    refgene_info = get_refgene_info(args.ref_gtf)
    for gtf in args.input_gtf:
        track_ref = get_replace_gene_dict(gtf,refgene_info)
        for k in track_ref:sys.stderr.write("%s\t%s\n"%(k,track_ref[k]))
        for line in fileinput.input(gtf,backup = ".orl",inplace = 1):
            if re.search('gene_id ".+?";',line):
                geneid = re.search('gene_id "(.+?)";',line).group(1)
                if geneid in track_ref:
                    print line.replace(geneid,track_ref[geneid]),
                else:
                    print line,
            else:
                print line,
                
if __name__ == "__main__":
    main()
