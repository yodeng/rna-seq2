#!/usr/bin/env python
#coding:utf-8
### 转录本组装及count值统计
### 
## 
import argparse,os,sys,glob
import RNA

parser = argparse.ArgumentParser(description="This Script is used to assemble transcripts and find novel trans")
parser.add_argument("-i","--input_dir",type=str,help="The input align directory",required = True)
parser.add_argument("-g","--gtf",type=str,help="the genome reference gtf file used for assemble transcripts",required = True)
parser.add_argument("-o","--output_dir",type=str,help="The output assemble directory",required = True)
parser.add_argument("-p","--parallel",type=int,help="the cpu number used, default: 8",required = False, default = 8)
parser.add_argument("-na","--noassemble",action = "store_false",default=True,help="whether assemble or not, default: assemble")
parser.add_argument("-s","--strand",action = "store_true",default= False,help="whether the data is from a strand-specific assay, default: False(fr-unstrand). if True, fr-firststrand will be done")
#parser.add_argument("-n","--novel_dir",type=str,help="the output of novel transcripts, if set.")
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
parser.add_argument("-t","--type",type=str,help="The expression quality value type. 'TPM','FPKM','FPKM_UQ' can be chosen. default: 'TPM FPKM FPKM_UQ'",default = ["TPM","FPKM","FPKM_UQ"],nargs="*")
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue  
            
def write_assemble_makeflow(makeflow,samples,aligndir,assembledir,gtf,thread,strand):
    samtools = RNA.get_bin_abspath("samtools")
    stringtie = RNA.get_bin_abspath("stringtie")
    fo = open(makeflow,"a+")
    for sn in samples:
        if not os.path.exists(os.path.join(aligndir,sn,"accepted_hits.bam")):
            fo.write("CATEGORY=Sam-sort_%s\n"%sn)
            fo.write("%s : %s\n" %(os.path.join(aligndir,sn,"accepted_hits.bam"),os.path.join(aligndir,sn,"accepted_hits.sam")))
            fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
            fo.write("\t%s sort -@ %d -o %s -O bam -T %s %s " %(samtools,thread,os.path.join(aligndir,sn,"accepted_hits.bam"),os.path.join(aligndir,sn,"temp"),os.path.join(aligndir,sn,"accepted_hits.sam")))
            fo.write("> %s 2> %s && rm -fr %s\n\n" %(os.path.join(aligndir,sn,"sort.out"),os.path.join(aligndir,sn,"sort.err"),os.path.join(aligndir,sn,"accepted_hits.sam")))
        fo.write("CATEGORY=assemble_%s\n"%sn)
        fo.write(os.path.join(assembledir,sn) + "/assemble.gtf : " + gtf + " " + os.path.join(aligndir,sn,"accepted_hits.bam") + "\n@BATCH_OPTIONS = -l h_vmem=25G\n")
        if strand:
            fo.write("\t%s %s -rf -p %d -G %s -o %s "%(stringtie,os.path.join(aligndir,sn,"accepted_hits.bam"),thread,gtf,os.path.join(assembledir,sn,"assemble.gtf")))
        else:
            fo.write("\t%s %s -p %d -G %s -o %s "%(stringtie,os.path.join(aligndir,sn,"accepted_hits.bam"),thread,gtf,os.path.join(assembledir,sn,"assemble.gtf")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(assembledir,sn,"stringtie.out"),os.path.join(assembledir,sn,"stringtie.err")))
    fo.write("CATEGORY=merge_all_gtf\n")
    fo.write(os.path.join(assembledir,"merge/merged.gtf : ") + " ".join([os.path.join(assembledir,s,"assemble.gtf") for s in samples]) + "\n")
    fo.write("\t%s --merge -m 100 -F 0.01 -G %s -p %d -o %s %s "%(stringtie,gtf,thread,os.path.join(assembledir,"merge/merged.gtf")," ".join([os.path.join(assembledir,s,"assemble.gtf") for s in samples])))
    fo.write("> %s 2> %s\n\n"%(os.path.join(assembledir,"merge/merged.out"),os.path.join(assembledir,"merge/merged.err")))
    for sn in samples:
        fo.write("CATEGORY=assemble_merge_%s\n"%sn)
        fo.write(os.path.join(assembledir,sn,"transcripts.gtf ") + os.path.join(assembledir,sn,"gene_abundance.txt : ") + os.path.join(assembledir,"merge/merged.gtf.orl") + " " + os.path.join(aligndir,sn,"accepted_hits.bam") + "\n@BATCH_OPTIONS = -l h_vmem=25G\n")
        if strand:
            fo.write("\t%s %s -rf -e -p %d -A %s -C %s -b %s -G %s -o %s "%(stringtie,os.path.join(aligndir,sn,"accepted_hits.bam"),thread,os.path.join(assembledir,sn,"gene_abundance.txt"),os.path.join(assembledir,sn,"transcripts_covered_by_reads.txt"),os.path.join(assembledir,sn,"Ballgown_table"),os.path.join(assembledir,"merge/merged.gtf"),os.path.join(assembledir,sn,"transcripts.gtf")))
        else:
            fo.write("\t%s %s -e -p %d -A %s -C %s -b %s -G %s -o %s "%(stringtie,os.path.join(aligndir,sn,"accepted_hits.bam"),thread,os.path.join(assembledir,sn,"gene_abundance.txt"),os.path.join(assembledir,sn,"transcripts_covered_by_reads.txt"),os.path.join(assembledir,sn,"Ballgown_table"),os.path.join(assembledir,"merge/merged.gtf"),os.path.join(assembledir,sn,"transcripts.gtf")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(assembledir,sn,"stringtie.out"),os.path.join(assembledir,sn,"stringtie.err")))
    fo.close()
    
def write_filter_gtf_makeflow(makeflow,assembledir,samples):
    filter_gtf = RNA.get_bin_abspath('filter_gtf.py')
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=filter_gtf\n")
    fo.write(" ".join([os.path.join(assembledir,sn,"transcripts.gtf.filter.gtf") for sn in samples]) + " : " + " ".join([os.path.join(assembledir,sn,"transcripts.gtf") for sn in samples]) + "\n")
    fo.write("\tpython %s %s\n\n"%(filter_gtf," ".join([os.path.join(assembledir,sn,"transcripts.gtf") for sn in samples])))
    fo.close()
    
def write_get_count_matrix_makeflow(makeflow,samples,assembledir):
    prepDE = RNA.get_bin_abspath("get_stringTie_count.py")
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_count_matrix\n")
    fo.write(os.path.join(assembledir,"gene_count_matrix.txt ") + os.path.join(assembledir,"transcript_count_matrix.txt : ") + " ".join([os.path.join(assembledir,sn,"transcripts.gtf.filter.gtf") for sn in samples]) + "\n")
    fo.write("\tpython %s -i %s -g %s -t %s &> %s\n\n"%(prepDE,assembledir,os.path.join(assembledir,"gene_count_matrix.txt"),os.path.join(assembledir,"transcript_count_matrix.txt"),os.path.join(assembledir,"get_count_matrix.log")))
    fo.close()

def write_change_gene_id_makeflow(makeflow,samples,gtf,ref_gtf):
    change_gene_id = RNA.get_bin_abspath("change_gene_id_in_transcript.py")
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=change_gene_id\n")
    fo.write(gtf + ".orl : " + gtf + " " + ref_gtf + "\n")
    fo.write("\tpython %s -i %s -g %s &> %s\n\n"%(change_gene_id,gtf,ref_gtf,os.path.join(os.path.dirname(gtf),"change_gene_id_in_transcript.log")))
    fo.close()

def write_get_expression_matrix_makeflow(makeflow,samples,assembledir):
    get_expression_matrix = RNA.get_bin_abspath("get_expression_matrix.py")
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_gene_expression_matrix\n")
    fo.write(os.path.join(assembledir,"gene.expression.txt : ") + " ".join([os.path.join(assembledir,sn,"gene_abundance.txt") for sn in samples]) + "\n")
    fo.write("\tpython %s -i %s -o %s -t TPM &> %s\n\n"%(get_expression_matrix," ".join([os.path.join(assembledir,sn,"gene_abundance.txt") for sn in samples]),os.path.join(assembledir,"gene.expression.txt"),os.path.join(assembledir,"get_gene_expression.log")))
    fo.write("CATEGORY=get_isoform_expression_matrix\n")
    fo.write(os.path.join(assembledir,"isoform.expression.txt : ") + " ".join([os.path.join(assembledir,sn,"transcripts.gtf.filter.gtf") for sn in samples]) + "\n")
    fo.write("\tpython %s -i %s -o %s -t TPM &> %s\n\n"%(get_expression_matrix," ".join([os.path.join(assembledir,sn,"transcripts.gtf.filter.gtf") for sn in samples]),os.path.join(assembledir,"isoform.expression.txt"),os.path.join(assembledir,"get_isoform_expression.log")))
    fo.close()
        
def write_summary_trans_stat_makeflow(makeflow,samples,assembledir):
    summary_trans = RNA.get_bin_abspath("summary_trans.py")
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=summary_transcripts\n")
    fo.write(os.path.join(assembledir,"assembly.state.txt ") + " ".join([os.path.join(assembledir,sn,name) for sn in samples for name in ["transcripts.gtf.filter.gtf.length_stat","transcripts.gtf.filter.gtf.transcript.length_distribution"]]) + " " + os.path.join(assembledir,"merge","merged.gtf.length_stat ") + os.path.join(assembledir,"merge","merged.gtf.transcript.length_distribution") +   " : " + " ".join([os.path.join(assembledir,sn,'transcripts.gtf.filter.gtf') for sn in samples]) + " " +  os.path.join(assembledir,"merge","merged.gtf") + "\n")
    fo.write("\tpython %s -i %s -o %s &> %s\n\n"%(summary_trans," ".join([os.path.join(assembledir,sn,'transcripts.gtf.filter.gtf') for sn in samples]) + " " +  os.path.join(assembledir,"merge","merged.gtf"),os.path.join(assembledir,"assembly.state.txt"),os.path.join(assembledir,"summary_trans.log")))
    fo.close()
    
# def write_filter_feature_makeflow(makeflow,assembledir,gtf,assemble):
    # filter_new_feature = RNA.get_bin_abspath("filter_new_feature.py")
    # fo = open(makeflow,"a+")
    # fo.write("CATEGORY=filter_new_feature\n")
    # fo.write(" ".join([os.path.join(assembledir,i) for i in ["gene_count_matrix.txt","gene.expression.txt","isoform.expression.txt","transcript_count_matrix.txt"]]) + " : " + " ".join([os.path.join(assembledir,i) for i in ["gene_count_matrix.txt","gene.expression.txt","isoform.expression.txt","transcript_count_matrix.txt"]]) + "\n")
    # if assemble:
        # fo.write("\tpython %s %s %s\n\n"%(filter_new_feature,gtf," ".join([os.path.join(assembledir,i) for i in ["gene_count_matrix.txt","gene.expression.txt","isoform.expression.txt","transcript_count_matrix.txt"]])))
    # else:
        # fo.write("\tpython %s %s %s -na\n\n"%(filter_new_feature,gtf," ".join([os.path.join(assembledir,i) for i in ["gene_count_matrix.txt","gene.expression.txt","isoform.expression.txt","transcript_count_matrix.txt"]])))
    # fo.close()

def write_get_expression_matrix_by_count_makeflow(makeflow,samples,assembledir,type,gtf):
    get_exp_by_count = RNA.get_bin_abspath("get_exp_by_count.py")
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_gene_expression_matrix_by_count\n")
    fo.write(" ".join([os.path.join(assembledir,"gene_"+i+"_matrix.txt") for i in type]) + " : " + os.path.join(assembledir,"gene_count_matrix.txt ") + gtf + "\n")
    fo.write("\tpython %s -c %s -r %s -t gene -v %s &> %s\n\n"%(get_exp_by_count,os.path.join(assembledir,"gene_count_matrix.txt"),gtf," ".join(type),os.path.join(assembledir,"get_gene_by_count.log")))
    fo.write("CATEGORY=get_transcript_expression_matrix_by_count\n")
    fo.write(" ".join([os.path.join(assembledir,"transcript_"+i+"_matrix.txt") for i in type]) + " : " + os.path.join(assembledir,"transcript_count_matrix.txt ") + gtf + "\n")
    fo.write("\tpython %s -c %s -r %s -t transcript -v %s &> %s\n\n"%(get_exp_by_count,os.path.join(assembledir,"transcript_count_matrix.txt"),gtf," ".join(type),os.path.join(assembledir,"get_trans_by_count.log")))
    fo.close()
    
def write_noassemble_makeflow(makeflow,samples,aligndir,assembledir,ref_gtf,thread,strand):
    samtools = RNA.get_bin_abspath("samtools")
    stringtie = RNA.get_bin_abspath("stringtie")
    fo = open(makeflow,"a+")
    for sn in samples:
        if not os.path.exists(os.path.join(aligndir,sn,"accepted_hits.bam")):
            fo.write("CATEGORY=Sam-sort_%s\n"%sn)
            fo.write("%s : %s\n" %(os.path.join(aligndir,sn,"accepted_hits.bam"),os.path.join(aligndir,sn,"accepted_hits.sam")))
            fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
            fo.write("\t%s sort -@ %d -o %s -O bam -T %s %s " %(samtools,thread,os.path.join(aligndir,sn,"accepted_hits.bam"),os.path.join(aligndir,sn,"temp"),os.path.join(aligndir,sn,"accepted_hits.sam")))
            fo.write("> %s 2> %s && rm -fr %s\n\n" %(os.path.join(aligndir,sn,"sort.out"),os.path.join(aligndir,sn,"sort.err"),os.path.join(aligndir,sn,"accepted_hits.sam")))
    for sn in samples:
        fo.write("CATEGORY=assemble_ref_%s\n"%sn)
        fo.write(os.path.join(assembledir,sn,"transcripts.gtf ") + os.path.join(assembledir,sn,"gene_abundance.txt : ") + ref_gtf + " " + os.path.join(aligndir,sn,"accepted_hits.bam") + "\n@BATCH_OPTIONS = -l h_vmem=25G\n")
        if strand:
            fo.write("\t%s %s -rf -e -p %d -A %s -C %s -b %s -G %s -o %s "%(stringtie,os.path.join(aligndir,sn,"accepted_hits.bam"),thread,os.path.join(assembledir,sn,"gene_abundance.txt"),os.path.join(assembledir,sn,"transcripts_covered_by_reads.txt"),os.path.join(assembledir,sn,"Ballgown_table"),ref_gtf,os.path.join(assembledir,sn,"transcripts.gtf")))
        else:
            fo.write("\t%s %s -e -p %d -A %s -C %s -b %s -G %s -o %s "%(stringtie,os.path.join(aligndir,sn,"accepted_hits.bam"),thread,os.path.join(assembledir,sn,"gene_abundance.txt"),os.path.join(assembledir,sn,"transcripts_covered_by_reads.txt"),os.path.join(assembledir,sn,"Ballgown_table"),ref_gtf,os.path.join(assembledir,sn,"transcripts.gtf")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(assembledir,sn,"stringtie.out"),os.path.join(assembledir,sn,"stringtie.err")))
    fo.close()
    
def main():
    aligndir = os.path.abspath(args.input_dir)
    gtf = os.path.abspath(args.gtf)
    samples = sorted(next(os.walk(aligndir))[1])
    assembledir = os.path.abspath(args.output_dir)
    assemble = args.noassemble
    if assemble:
        mkdir(assembledir,samples + ["merge"])
        write_assemble_makeflow(args.makeflow,samples,aligndir,assembledir,gtf,args.parallel,args.strand)
        write_filter_gtf_makeflow(args.makeflow,assembledir,samples)
        write_summary_trans_stat_makeflow(args.makeflow,samples,assembledir)
        write_change_gene_id_makeflow(args.makeflow,samples,os.path.join(assembledir,"merge/merged.gtf"),gtf)
        write_get_count_matrix_makeflow(args.makeflow,samples,assembledir)
        write_get_expression_matrix_makeflow(args.makeflow,samples,assembledir)
        write_get_expression_matrix_by_count_makeflow(args.makeflow,samples,assembledir,args.type,os.path.join(assembledir,"merge/merged.gtf"))
    else:
        mkdir(assembledir,samples)
        write_noassemble_makeflow(args.makeflow,samples,aligndir,assembledir,gtf,args.parallel,args.strand)
        write_filter_gtf_makeflow(args.makeflow,assembledir,samples)
        write_get_count_matrix_makeflow(args.makeflow,samples,assembledir)
        write_get_expression_matrix_makeflow(args.makeflow,samples,assembledir)
        write_get_expression_matrix_by_count_makeflow(args.makeflow,samples,assembledir,args.type,gtf)
    
if __name__ == "__main__":
    main()


    
