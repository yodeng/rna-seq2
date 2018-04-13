#!/usr/bin/env python
#coding:utf-8

import sys,os
fo = open(sys.argv[-1],"w")
fo.write("Class\tTotal Assembled Transcripts\tTotal Assembled Loci\tIn Known Transcripts\tCorresponding Known Loci\tIn Novel Transcripts\tCorresponding Novel Loci\n")
for f in sys.argv[1:-1]:
    s = os.path.basename(os.path.dirname(f))
    with open(f) as fs:
        l = fs.readlines()
        fo.write(s+"\t"+"\t".join([i.strip().split("\t")[-1] for i in l]) + "\n")
fo.close()
    

