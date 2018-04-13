#!/usr/bin/env python
#coding:utf-8
import sys,re,argparse,os,glob
from collections import defaultdict
from math import ceil

parser = argparse.ArgumentParser(description="Generates tab-sep files containing the count matrices for genes or transcripts, using the coverage values found in the output of `stringtie -e`. IDs with all zero value in samples will not be output. if input a directory, only $input_dir/*/transcripts.gtf.filter.gtf will be match")
parser.add_argument('-i', '--input', required = True,help="the parent directory of the sample sub-directories or a textfile listing the paths to GTF files, if input a file, sample_name gtf_path info must be in one row")
parser.add_argument('-l', '--length', default=150, type=int, help="the average read length [default: 150]")
parser.add_argument('-g', '--gene', default='gene_count_matrix.txt', help="where to output the gene count matrix [default: gene_count_matrix.txt]")
parser.add_argument('-t', '--trans' ,default='transcript_count_matrix.txt', help="where to output the transcript count matrix [default: transcript_count_matrix.txt]")
parser.add_argument('-r', '--remove' ,default=False,action="store_true", help="weather to remove id with any zero value in samples, default: no_remove")
args = parser.parse_args()

RE_GENE_ID=re.compile('gene_id "([^"]+?)"')
RE_GENE_NAME=re.compile('gene_name "([^"]+?)"')
RE_TRANSCRIPT_ID=re.compile('transcript_id "([^"]+?)"')
RE_COVERAGE=re.compile('cov "([\-\+\d\.]+?)"')

def getGeneID(s, ctg, tid):
    r=RE_GENE_ID.search(s)
    if r: return r.group(1)
    r=RE_GENE_NAME.search(s)
    if r: return ctg+'|'+r.group(1)
    return tid
  
def getCov(s):
    r=RE_COVERAGE.search(s)
    if r:
        v=float(r.group(1))
        if v<0.0: v=0.0
        return v
    return 0.0

def get_sample_info(input):
    samples = [] # List of tuples. If sample list, (first column, path). Else, (subdirectory name, path to gtf file in subdirectory)
    if (os.path.isfile(input)):
        # gtfList = True
        try:
            fin = open(input, 'r')
            for line in fin:
                if line.strip() and not line.startswith("#"):
                    lineLst = tuple(line.strip().split()[0:2])
                    if (len(lineLst) != 2):
                        print "Error: Text file with sample ID and path invalid (%s)" % (line.strip())
                        exit(1)
                    if lineLst[0] in samples:
                        print "Error: Sample ID duplicated (%s)" % (lineLst[0])
                        exit(1)
                    if not os.path.isfile(lineLst[1]):
                        print "Error: GTF file not found (%s)" % (lineLst[1])
                        exit(1)
                    samples.append(lineLst)
        except IOError:
            print "Error: List of .gtf files, %s, doesn't exist" % (input)
            sys.exit(1)
    else:
        # gtfList = False
        ## Check that opts.input directory exists
        if not os.path.isdir(input):
            parser.print_help()
            print " "
            print "Error: sub-directory '%s' not found!" % (input)
            sys.exit(1)
        samples = [(i,glob.glob(os.path.join(input,i,"transcripts.gtf.filter.gtf"))[0]) for i in next(os.walk(input))[1] if len(glob.glob(os.path.join(input,i,"transcripts.gtf.filter.gtf")))]
    return sorted(samples)
    
def get_trans_dict(samples):    
    geneIDs={}     ### key=transcriptID, values=geneID
    for i,s in enumerate(samples):
        with open(s[1]) as f:
            split = [l.split("\t") for l in f.readlines()]
        for v in split:
            if len(v)>2 and v[2]=="transcript":
                t_id=RE_TRANSCRIPT_ID.search(v[-1]).group(1)
                g_id=getGeneID(v[-1], v[0], t_id)
                geneIDs.setdefault(t_id,g_id)
    return geneIDs
    
def get_count(samples,read_len,geneIDs):
    transDict={}   ### key=transcriptID, values=sum-count
    geneDict={}
    for q,s in enumerate(samples):
        print q, s[0]
        f = open(s[1])
        transcript_len=0
        for l in f:
            if l.startswith("#"): continue
            v=l.split('\t')
            if v[2]=="transcript":
                if transcript_len>0:
                    transDict.setdefault(t_id, defaultdict(int))
                    transDict[t_id][s[0]] = int(ceil(coverage*transcript_len/read_len))
                t_id=RE_TRANSCRIPT_ID.search(v[-1]).group(1)
                g_id=getGeneID(v[-1], v[0], t_id)
                coverage=getCov(v[-1])
                transcript_len=0
            if v[2]=="exon":transcript_len+=int(v[4])-int(v[3])+1
        transDict.setdefault(t_id, defaultdict(int))
        transDict[t_id].setdefault(s[0], int(ceil(coverage*transcript_len/read_len)))
        for i,v in transDict.iteritems():
            geneDict.setdefault(geneIDs[i],defaultdict(int))
            geneDict[geneIDs[i]][s[0]] += v[s[0]]     
    return transDict,geneDict

def main():
    samples = get_sample_info(args.input)
    geneIDs = get_trans_dict(samples)
    transDict,geneDict = get_count(samples,args.length,geneIDs)
    if args.trans:
        with open(args.trans, 'w') as t_outfile:
            t_outfile.write("transcript_id\t" + "\t".join([x for x,y in samples]) +"\n")
            if args.remove:
              for i in transDict:
                if all([transDict[i][s] for s,v in samples]):           ## remove any zero value ids
                    t_outfile.write(i + "\t" + "\t".join([str(transDict[i][s]) for s,v in samples]) + "\n")
            else:
              for i in transDict:
                  if any(([transDict[i][s] for s,v in samples])):   ##  remove all zero value ids
                    t_outfile.write(i + "\t" + "\t".join([str(transDict[i][s]) for s,v in samples]) + "\n")
    if args.gene:
        with open(args.gene, 'w') as g_outfile:
            g_outfile.write("gene_id\t" + "\t".join([x for x,y in samples]) +"\n")
            if args.remove:
              for i in geneDict:
                if all([geneDict[i][s] for s,v in samples]):            ## remove any zero value ids
                    g_outfile.write(i + "\t" + "\t".join([str(geneDict[i][s]) for s,v in samples]) + "\n")
            else:
              for i in geneDict:
                  if any([geneDict[i][s] for s,v in samples]):      ##  remove all zero value ids
                    g_outfile.write(i + "\t" + "\t".join([str(geneDict[i][s]) for s,v in samples]) + "\n")
                
if __name__ == "__main__":
    main()
        
