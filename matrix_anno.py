#!/usr/bin/env python
#coding:utf-8

import argparse,re,os
from math import log

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is get all gene anno table from total gtf file and anno file")
    parser.add_argument("-t","--table",type=str,help="The input count or fpkm table file.",required = True)
    parser.add_argument("-anno","--anno",type=str,help="The input ref anno file.",required = True)
    parser.add_argument("-gtf","--gtf",type=str,help="The input gtf file.",required = True)
    parser.add_argument("-type","--type",type=str,help="The type of input table file.",choices=["gene","trans"],required = True)
    parser.add_argument("-o","--output",type=str,help="The novel trans output dir",required = True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def get_gene_trans_dict(gtf):
    genedict = {}   ## trackgene:refgene
    transdict = {}  ## tracktrans:trackgene
    genestrand = {}
    with open(gtf) as fi:
      for line in fi:
        if re.search('transcript_id ".+?";',line):
            strand = line.split("\t")[6]
            geneid = re.search('gene_id "(.+?)";',line).group(1)
            transid = re.search('transcript_id "(.+?)";',line).group(1)
            genestrand[geneid] = strand
            if re.search('ref_gene_id ".+?";',line):
                genedict.setdefault(geneid,set()).add(re.search('ref_gene_id "(.+?)";',line).group(1))
                transdict[transid] = re.search('ref_gene_id "(.+?)";',line).group(1)
                genestrand[re.search('ref_gene_id "(.+?)";',line).group(1)] = strand
            else:
                genedict[geneid] = set([geneid])
    return genedict,transdict,genestrand

def get_anno_dict(anno):
    anno_dict = {}
    with open(anno) as fi:
        for line in fi:
            if line.strip() and not line.startswith("#"):
                k = line.split("\t")[0]
                anno_dict[k] = line.strip().split("\t")
    return anno_dict

def main():
    args = parseArg()
    genedict,transdict,genestrand = get_gene_trans_dict(args.gtf)
    anno_dict = get_anno_dict(args.anno)
    if args.type == "gene":
        with open(args.table) as fi, open(args.output,"w") as fo:
            header = fi.next()
            fo.write(header.strip() + '\tRef_gene\tRef_gene_info\tEntrez\tUniprot\tSymbol\tDescription\n')
            for line in fi:
                linelist = line.strip().split("\t")
                id = linelist[0]
                if id in genedict:
                    only_ref = [refg for refg in genedict[id] if refg in anno_dict]
                    Ref_gene = "|".join(only_ref)
                    Ref_gene_info = "|".join([anno_dict[g][3] + ":" + genestrand[id] + ":" + anno_dict[g][4] + "-" + anno_dict[g][5] for g in only_ref])
                    Entrez = "|".join([anno_dict[g][1] for g in only_ref])
                    Uniprot = "|".join([anno_dict[g][2] for g in only_ref])
                    Symbol = "|".join([anno_dict[g][-2] for g in only_ref])
                    Description = "|".join([anno_dict[g][-1] for g in only_ref])
                    if not only_ref:Ref_gene,Ref_gene_info,Entrez,Uniprot,Symbol,Description = list("-")*6
                else:
                    Ref_gene,Ref_gene_info,Entrez,Uniprot,Symbol,Description = list("-")*6
                fo.write(line.strip() + "\t" + "\t".join([Ref_gene,Ref_gene_info,Entrez,Uniprot,Symbol,Description]) + "\n")  
    else:
        with open(args.table) as fi, open(args.output,"w") as fo:
            header = fi.next()
            fo.write(header.strip() + '\tRef_gene\tRef_gene_info\tEntrez\tUniprot\tSymbol\tDescription\n')
            for line in fi:
                linelist = line.strip().split("\t")
                id = linelist[0]
                if id in transdict and transdict[id] in anno_dict:
                    Ref_gene = transdict[id]
                    Ref_gene_info = anno_dict[transdict[id]][3] + ":" + genestrand[transdict[id]] + ":" + anno_dict[transdict[id]][4] + "-" + anno_dict[transdict[id]][5]
                    Entrez = anno_dict[transdict[id]][1]
                    Uniprot = anno_dict[transdict[id]][2]
                    Symbol = anno_dict[transdict[id]][-2]
                    Description = anno_dict[transdict[id]][-1]
                else:
                    Ref_gene,Ref_gene_info,Entrez,Uniprot,Symbol,Description = list("-")*6
                fo.write(line.strip() + "\t" + "\t".join([Ref_gene,Ref_gene_info,Entrez,Uniprot,Symbol,Description]) + "\n")  
                    
                    
                
if __name__ == "__main__":
    main()

