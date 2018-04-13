#!/usr/bin/env python
#coding:utf-8
### 比对步骤
## cd /lustre/work/yongdeng/Project && python ~/align_clean.py -i 160761_test/QC/ -g 160761_test/genome/ -mf 160761_test/align_clean.mf -o 160761_test/align
import argparse,os,sys,glob,re
import RNA

parser = argparse.ArgumentParser(description="This Script is used to align cleandata to reference genome by hisat2 program.")
parser.add_argument("-i","--input_dir",type=str,help="The input cleandata directory",required = True)
#parser.add_argument("-g","--genome_dir",type=str,help="the genome reference dir, 'fa' and 'gtf' file must under this dir, and align index will be generate in this dir",required=True)
parser.add_argument("-fa","--fa",type=str,help="the genome reference fasta file",required=True)
parser.add_argument("-gtf","--gtf",type=str,help="the genome reference gtf file",required=True)
parser.add_argument("-p","--parallel",type=int,help="the cpu number used",required = False, default = 8)
parser.add_argument("-o","--output_dir",type=str,help="The output ailgned directory",required = True)
parser.add_argument("-a","--assemble",action = "store_true",default= False,help="whether assemble or not, '--assemble' means assemble will be done in the next, default: no assemble if not set")
parser.add_argument("-s","--strand",action = "store_true",default= False,help="whether the data is from a strand-specific assay, default: False if no set")
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue           

def get_hisat_version():
    hisat2 = RNA.get_bin_abspath("hisat2")
    version = os.popen(hisat2 + " --version").readline().split()[-1]
    return version

def write_align_makeflow(samples,makeflow,input_d,output_d,thread,assemble,strand,fa,gtf):
    hisat2_extract_splice_sites = RNA.get_bin_abspath("hisat2_extract_splice_sites.py")
    hisat2_extract_exons = RNA.get_bin_abspath("hisat2_extract_exons.py")
    hisat2_build = RNA.get_bin_abspath("hisat2-build")
    hisat2 = RNA.get_bin_abspath("hisat2")
    reflat = RNA.get_bin_abspath("create_refFlat.pl")
    #if not os.path.exists(genome_dir):
    #    print "%s dir not exists, please check!" %genome_dir
    #    sys.exit(0)
    #gtf = glob.glob(os.path.join(genome_dir,"*.gtf"))[0]
    #fa = glob.glob(os.path.join(genome_dir,"*.fa"))[0]
    #if not fa or not gtf:
    #    print "please check the fa and gtf file in %s" %genome_dir
    #    sys.exists(0)
    if not os.path.exists(fa) or not os.path.exists(gtf):
        print "fa or gtf file not exists"
        sys.exit(1)
    fa = os.path.abspath(fa)
    gtf = os.path.abspath(gtf)
    gtfdir = os.path.dirname(gtf)
    version = get_hisat_version()
    indexdir = os.path.join(os.path.dirname(fa),"hisat2_" + version + "_index")
    if not os.path.isdir(indexdir ): 
        try:
            os.makedirs(indexdir)
        except:
            print "Error: permission denied, could not make fasta index under fa dir, could not makedir %s"%indexdir
            sys.exit(1)
    genome_prefix = os.path.join(indexdir,os.path.splitext(os.path.basename(fa))[0] + '.hisat')
    fo =open(makeflow,"a+")
    if not os.path.exists(indexdir+"/hisat_build-genome.log"):
        fo.write("CATEGORY=extract_splice_sites\n")
        fo.write(gtf+".ss : " + gtf + "\n")
        fo.write("\tpython %s %s > %s.ss 2> %s\n\n"%(hisat2_extract_splice_sites,gtf,gtf,gtfdir+"/hisat2_extract_splice_sites.err"))
        fo.write("CATEGORY=hisat2_extract_exons\n")
        fo.write(gtf+".exon : " + gtf + "\n")
        fo.write("\tpython %s %s > %s.exon 2> %s\n\n"%(hisat2_extract_exons,gtf,gtf,gtfdir+"/hisat2_extract_exons.err"))
        fo.write("CATEGORY=hisat2-build\n")
        fo.write("%s : %s %s %s\n"%(indexdir+"/hisat_build-genome.log",fa,gtf+".ss", gtf+".exon"))
        #fo.write("\t%s -p %d --ss %s --exon %s %s %s &> %s\n\n" %(hisat2_build,thread,gtf+".ss",gtf+".exon",fa,genome_prefix,genome_dir+"/hisat_build-genome.log  || rm -fr " + genome_dir + "/hisat_build-genome.log"))
        fo.write("\t%s -p %d %s %s &> %s\n\n" %(hisat2_build,thread,fa,genome_prefix,indexdir+"/hisat_build-genome.log || mv " + indexdir + "/hisat_build-genome.log " + indexdir + "/hisat_build-genome.log.err"))
    if not os.path.exists(fa + ".hdrs"):
        spec = "_".join((re.split("[^a-zA-Z0-9]",os.path.basename(fa))[:2]))
        fo.write("CATEGORY=create_hdrs\n")
        fo.write(fa + ".hdrs : " + fa + "\n")
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/flow/create_hdrs.pl %s %s %s &> %s\n\n"%(fa,spec,fa+".hdrs",os.path.join(os.path.dirname(fa),"create_hdrs.log")))
    if not os.path.exists(gtf + ".refFlat"):
        fo.write("CATEGORY=create_refFlat\n")
        fo.write(gtf + ".refFlat : " + gtf + "\n")
        fo.write("\tperl %s %s %s &> %s\n\n"%(reflat,gtf,gtf + ".refFlat",gtfdir+"/refflat.err"))
    for sn in samples:
        r1 = os.path.join(input_d,sn,"R1.clean.fq.gz")
        r2 = os.path.join(input_d,sn,"R2.clean.fq.gz")
        od = os.path.join(output_d,sn)
        fo.write("CATEGORY=hisat2-align_%s\n"%sn)
        #fo.write(os.path.join(od,"hisat2.err ") + os.path.join(od,"accepted_hits.sam : ") + r1 + " " + r2 + " " + gtf+".ss " + indexdir + "/hisat_build-genome.log\n")
        fo.write(os.path.join(od,"accepted_hits.sam : ") + r1 + " " + r2 + " " + gtf+".ss " + indexdir + "/hisat_build-genome.log\n")
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        if strand:
            fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/checkpaired.py %s %s &> %s && %s --fr --rna-strandness RF %s-p %d --known-splicesite-infile %s -x %s -1 %s -2 %s -S %s "%(r1,r2,os.path.join(od,"hisat2.err"),hisat2,assemble,thread,gtf+".ss",genome_prefix,r1,r2,os.path.join(od,"accepted_hits.sam")))
        else:
            fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/checkpaired.py %s %s &> %s && %s --fr %s-p %d --known-splicesite-infile %s -x %s -1 %s -2 %s -S %s "%(r1,r2,os.path.join(od,"hisat2.err"),hisat2,assemble,thread,gtf+".ss",genome_prefix,r1,r2,os.path.join(od,"accepted_hits.sam")))
        fo.write("> %s 2> %s\n\n" %(os.path.join(od,"hisat2.out"),os.path.join(od,"hisat2.err")))
    fo.close()

def main():
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    fa = os.path.abspath(args.fa)
    gtf = os.path.abspath(args.gtf)
    samples = sorted(next(os.walk(input_dir))[1])
    mkdir(output_dir,samples)
    assemble = "--dta " if args.assemble else ""
    if not os.path.exists(os.path.abspath(os.path.dirname(args.makeflow))):
        os.makedirs(os.path.dirname(args.makeflow))
    write_align_makeflow(samples,args.makeflow,input_dir,output_dir,args.parallel,assemble,args.strand,fa,gtf)

if __name__=="__main__":
    main()
    
