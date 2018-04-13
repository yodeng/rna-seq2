#!/usr/bin/env python
#coding:utf-8

import argparse,re,os
from math import log

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is draw violon plot of genes or trans fpkm plot")
    parser.add_argument("-i","--input",type=str,help="The input fpkm table",required = True)
    #parser.add_argument("-t","--type",type=str,help="The type of input fpkm table, choices=['trans','genes']",,choices=["trans","gene"],required = True)
    parser.add_argument("-o","--outfile",type=str,help="The output file",required = True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def draw_plot(fpkm,out):
    violindict = {}
    with open(fpkm) as fi:
        header = fi.next().strip().split("\t")[1:]
        for line in fi:
            line = line.strip().split("\t")[1:]
            for n,i in enumerate(line):
                violindict.setdefault(header[n],[]).append(float(i))
    with open(fpkm+".violin.txt","w") as fo:
        fo.write("Sample\tlog2Fpkm\n")
        for s in violindict:
            if s == "length":continue
            for f in violindict[s]:
                if float(f) < 1:continue
                fo.write(s + "\t" + str(log(float(f),2)) + "\n")
    os.system("Rscript /lustre/work/yongdeng/software/eukaryon/flow/draw_plot_violin_log_sample.r %s %s"%(fpkm+".violin.txt",out))
            
if __name__ == "__main__" :
    args=parseArg()
    draw_plot(args.input,args.outfile)
    
            
        
    
