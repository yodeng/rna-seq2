#!/usr/bin/env python
#coding:utf-8
import sys,re
#  fpkm>0.01且长度大于100的转录本保留
for filename in sys.argv[1:]:
    with open(filename) as fi,open(filename+".filter.gtf","w") as fo:
        transcripts = []
        for line in fi:
            if line.startswith("#"):
                fo.write(line)
                continue
            if line.split("\t")[2] == "transcript":
                if not len(transcripts):
                    transcripts = [line,]
                    continue
                trans = transcripts[0]
                fpkm = re.search('FPKM "(.+?)";',trans).group(1)
                length = abs(int(trans.split("\t")[3]) - int(trans.split("\t")[4]))
                if float(fpkm) > 0.01 and length >= 100:
                    fo.writelines(transcripts)
                transcripts = [line,]
            else:
                transcripts.append(line)
        trans = transcripts[0]
        fpkm = re.search('FPKM "(.+?)";',trans).group(1)
        length = abs(int(trans.split("\t")[3]) - int(trans.split("\t")[4]))
        if float(fpkm) >= 0.01 and length >= 100: fo.writelines(transcripts)
                
                
