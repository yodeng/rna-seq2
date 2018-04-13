#!/usr/bin/env python
#coding:utf-8
  
import os,argparse

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to annotate deg result file")
    parser.add_argument("-i","--input",type=str,help="The input deg result file from DEG analysis, like deg.down.txt or deg.up.txt files",required = True,nargs="+")
    parser.add_argument("-sg","--sample_group",type=str,help="The all samples and groups rep. like groupA=sample1,sample2,... groupB=sample3,sample4,...",required = True,nargs="+")
    parser.add_argument("-nc","--normalcount",type=str,help="The normalized count file",required = True)
    parser.add_argument("-anno","--anno",type=str,help="the track_ref anno file",required = True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def groupdict(samplegroup):
    d = {}
    for rep in samplegroup:
        g = rep.split("=")[0]
        s = rep.split("=")[-1].split(",")
        d[g] = s
    return d
    
def annofile(anno):
    d = {}
    with open(anno) as an:
        header  = an.next()
        for line in an:
            k = line.split("\t")[0]
            v = line.strip().split("\t")[1:]
            d[k] = v
    return d
    
def countfile(count):
    d = {}
    with open(count) as fc:
        samples = fc.next().strip().split("\t")[1:]
        for line in fc:
            k = line.split("\t")[0]
            td = dict(zip(samples,line.strip().split("\t")[1:]))
            d[k] = td
    return d
    
def main():
    args = parseArg()
    groupDict = groupdict(args.sample_group)
    countDict = countfile(args.normalcount)
    annDict = annofile(args.anno)
    for infile in args.input:
        with open(infile) as fi,open(os.path.splitext(infile)[0] + ".anno.txt","w") as fo:
            header = fi.next().strip().split("\t")
            group1 = header[2].partition("_")[-1]
            group2 = header[3].partition("_")[-1]
            fo.write("ID\t%s\tlog2FoldChange\tpvalue\tpadj"%("\t".join(groupDict[group1] + groupDict[group2])) + '\tGene_info\tRef_gene\tRef_gene_info\tEntrez\tUniprot\tSymbol\tDescription\n')
            for line in fi:
                linelist = line.strip("\n").split("\t")
                ID = linelist[0]
                group12_count = [countDict[ID][i] for i in groupDict[group1]+groupDict[group2]]
                logfc,pval,padj = linelist[4:7]
                annoline = "\t".join(annDict[ID])
                fo.write(ID + "\t" + "\t".join(group12_count) + "\t" + "\t".join(linelist[4:7]) + "\t" + annoline + "\n") 

if __name__ == "__main__":
    main()
            
            
