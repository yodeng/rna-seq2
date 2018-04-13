#!/usr/bin/env python
#coding:utf-8

import argparse,re,os

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def parseArg():
    parser = argparse.ArgumentParser(description="This Script used for novel transcrtpts analysis")
    parser.add_argument("-i","--input",type=str,help="The gtf info. must be like 'sampleanme=/path/to/gtf ...'",required = True,nargs="+")
    parser.add_argument("-r","--refgtf",type=str,help="The reference gtf file",required = True)
    parser.add_argument("-o","--outputdir",type=str,help="The novel trans output dir")
    parser.add_argument("-mf","--makeflow",type=str,help="The output makeflow file")
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def get_gtf_info(info):
    d = {}
    for g in info:
        k = g.strip().split("=")[0].strip()
        v = g.strip().split("=")[-1].strip()
        d[k] = os.path.abspath(v)
    return d
    
    
def write_novel_trans_makeflow(makeflow,infodict,refgtf,outdir):
    fo = open(makeflow,"a+")
    copyinfo = infodict.copy()
    mergename = [i for i in infodict if "merge" in i][0]
    mergepath = copyinfo.pop(mergename)
    
    for s in copyinfo:
        fo.write("CATEGORY=change_id_name_%s\n"%s)
        fo.write(os.path.abspath(infodict[s]) + ".orl : " + os.path.abspath(infodict[s]) + "\n")
        fo.write("\tpython /lustre/work/yongdeng/software/eukaryon/flow/change_gene_id_in_transcript.py -i %s -g %s &> %s\n\n"%(infodict[s],refgtf,os.path.join(outdir,s,"change_gene_id.log")))
        fo.write("CATEGORY=find_novel_%s\n"%s)
        fo.write(os.path.join(outdir,s,"novel_transcripts.gtf ") + os.path.join(outdir,s,"transcripts.gtf.known_and_novel_stat : ") + infodict[s] + ".orl\n")
        fo.write("\tpython /lustre/work/yongdeng/software/eukaryon/flow/find_novel_isoform.py -gtf %s -ref %s -o %s &> %s\n\n"%(infodict[s],refgtf,os.path.join(outdir,s),os.path.join(outdir,s,"novel_transcripts.gtf.log")))
    
    fo.write("CATEGORY=find_novel_merge\n")
    fo.write(os.path.join(outdir,"merge","novel_transcripts.gtf ") + os.path.join(outdir,"merge","transcripts.gtf.known_and_novel_stat : ") + infodict["merge"] + ".orl\n")
    fo.write("\tpython /lustre/work/yongdeng/software/eukaryon/flow/find_novel_isoform.py -gtf %s -ref %s -o %s &> %s\n\n"%(infodict["merge"],refgtf,os.path.join(outdir,"merge"),os.path.join(outdir,"merge","novel_transcripts.gtf.log")))

    fo.write("CATEGORY=stat_novel_trans\n")
    fo.write(os.path.join(outdir,"novel_transcripts.state.txt : ") + " ".join([os.path.join(outdir,s,"transcripts.gtf.known_and_novel_stat") for s in infodict]) + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/eukaryon/flow/stat_novel.py " + " ".join([os.path.join(outdir,s,"transcripts.gtf.known_and_novel_stat") for s in sorted(copyinfo.keys())]) + " " + os.path.join(outdir,mergename,"transcripts.gtf.known_and_novel_stat") + " " + os.path.join(outdir,"novel_transcripts.state.txt") + "\n\n")   
    fo.close()
    
def main():
    args = parseArg()
    samplegtf = get_gtf_info(args.input)
    mkdir(args.outputdir,samplegtf.keys())
    write_novel_trans_makeflow(args.makeflow,samplegtf,os.path.abspath(args.refgtf),os.path.abspath(args.outputdir))
    
if __name__ == "__main__":
    main()
    

        
