#!/usr/bin/env python
#coding:utf-8

import sys,re,os

gtf = sys.argv[1]
ref = set()
for line in open(gtf):
    if re.search('gene_id ".+?";',line):
        ref.add(re.search('gene_id "(.+?)";',line).group(1))
    if re.search('transcript_id ".+?";',line):
        ref.add(re.search('transcript_id "(.+?)";',line).group(1))
        
# for line in fileinput.input(sys.argv[2:],backup = ".all",inplace = 1)
    # if fileinput.isfirstline():
        # print line,
    # else:
        # if line.split("\t")[0] in ref:
            # print line,

if sys.argv[-1] == "-na":
    for filename in sys.argv[2:-1]: 
        outfile = os.path.splitext(filename)[0]
        with open(filename) as fi,open(outfile,"w") as fo:
            header = fi.next()
            fo.write(header)
            for line in fi:
                if line.split("\t")[0] in ref:
                    fo.write(line)
else:
    for filename in sys.argv[2:]:
        outfile = os.path.splitext(filename)[0]
        with open(filename) as fi,open(outfile,"w") as fo:
            fo.write(fi.read())