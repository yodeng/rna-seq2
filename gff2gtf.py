#!/usr/bin/env python
#coding:utf-8
## gffread -T GCF_000313675.1_ASM31367v1_genomic.gff -o GCF_000313675.1_ASM31367v1_genomic.gtf -F
### python scripts.py gff gtf

### 考虑所有的情况之后，对真核生物gff文件转化为可用的gtf文件适用。

import sys,re,os
from collections import OrderedDict

gene_dict = OrderedDict()
no_gene_dict = {}
gene_locus = {}
tran_locus = {}
genetype = {}

def change_product(product):
    return product.replace('%3B',";").replace('%20'," ").replace('%28',"(").replace('%29',")").replace('%2F',"/").replace('%2C',',')

for line in open(sys.argv[1]):
    line_list = line.strip().split("\t")
    if re.search('ID=gene.+?;',line_list[-1]):
        gene_locus[re.search('ID=(gene.+?);',line_list[-1]).group(1)] = line_list[3:8]
        genetype[re.search('ID=(gene.+?);',line_list[-1]).group(1)] = re.search('gene_biotype=(.+?);',line_list[-1]+";").group(1)
    if re.search('ID=rna.+?;',line_list[-1]):
        tran_locus[re.search('ID=(rna.+?);',line_list[-1]).group(1)] = line_list[3:8]

for line in os.popen('gffread -T -F %s -o-'%sys.argv[1]):
    if re.search('gene_id "(gene.+?)";',line):
        gene_dict.setdefault(re.search('gene_id "(gene.+?)";',line).group(1),[]).append(line)       
    else:
        no_gene_dict.setdefault(re.search('transcript_id "(rna.+?)";',line).group(1),[]).append(line)   ### 该基因在gff文件中没有gene_id

gene_id_list = []

fo = open(sys.argv[-1],"w")
n = 0    
for i in no_gene_dict:
    chr = no_gene_dict[i][0].split("\t")[0]    ### 只考虑长度为1的值列表
    source = no_gene_dict[i][0].split("\t")[1]
    type = re.search('gbkey "(.+?)";',no_gene_dict[i][0]).group(1)
    product = re.search('product ".+?";',no_gene_dict[i][0]).group(0) + " " if re.search('product ".+?";',no_gene_dict[i][0]) else ""
    product = change_product(product)
    note = re.search('Note ".+?";',no_gene_dict[i][0]).group(0) + " " if re.search('Note ".+?";',no_gene_dict[i][0]) else ""
    # if product:
        # gene_id = 'gene_id "%s";'%product.split('"')[1] + " " 
    # else:
        # n += 1
        # gene_id = 'gene_id "%s_%d";'%(type,n) + " "
    n += 1
    gene_id = 'gene_id "%s_%d";'%(type,n) + " "

    transcript_id = 'transcript_id "' + gene_id.split('"')[1] + '_T";' + " "
    gene_biotype = 'gene_biotype "%s";'%type + " "
    product = re.search('product ".+?";',no_gene_dict[i][0]).group() + " " if re.search('product ".+?";',no_gene_dict[i][0]) else ""
    product = change_product(product)
    gene_pos = tran_locus[re.search('transcript_id "(.+?)";',no_gene_dict[i][0]).group(1)]
    trans_pos = gene_pos
    entrez_id = 'Entrez_id "' + re.search('GeneID:(.+?)\W',no_gene_dict[i][0]).group(1) + '";' + " " if re.search('GeneID:(.+?)\W',no_gene_dict[i][0]) else ""
    locus = no_gene_dict[i][0].split("\t")[3:8]
    fo.write("\t".join([chr,source,"gene"]) + "\t" + "\t".join(gene_pos) + "\t" + (gene_id + gene_biotype + entrez_id).strip()  + "\n"  ) 
    fo.write("\t".join([chr,source,type]) + "\t" + "\t".join(gene_pos) + "\t" + (gene_id + transcript_id + gene_biotype + entrez_id).strip() + "\n"      ) 
    fo.write("\t".join([chr,source,"exon"]) + "\t" + "\t".join(locus) + "\t" + (gene_id + transcript_id + gene_biotype + product + note + entrez_id + 'exon_number "1";').strip()  + "\n")
    
for i in gene_dict:
    gene_dict[i] = sorted(gene_dict[i],key=lambda x:x.split("\t")[3]) if gene_dict[i][0].split("\t")[6] == "+" else sorted(gene_dict[i],key=lambda x:x.split("\t")[3],reverse = True)
    trans = set([t for t in re.findall('transcript_id "(.+?)";',str(gene_dict[i])) if "rna" in t])
    if len(trans) == 1:   # 一个基因只有一个转录本情况
        n = 0
        chr = gene_dict[i][0].split("\t")[0]
        source = gene_dict[i][0].split("\t")[1]
        type = re.findall('gbkey "(.+?)";',str(gene_dict[i]))
        # if "mRNA" in type or "CDS" in type:
            # type = "mRNA"
        # else:
            # type = "transcript" 
        # #gene_biotype = re.search('gene_biotype ".+?";',str(gene_dict[i])).group() + " "
        type = "mRNA" if "mRNA" in type or "CDS" in type else genetype[re.search('gene_id "(.+?)";',str(gene_dict[i])).group(1)]
        entrez_id = 'Entrez_id "' + re.search('GeneID:(.+?)\W',str(gene_dict[i])).group(1) + '";' + " " if re.search('GeneID:(.+?)\W',str(gene_dict[i])) else ""

        id = re.search('gene_name "(.+?)";',str(gene_dict[i])).group(1)  ### 以gene_name作为gene_id,即原gff文件中的"gene=(.+?);"或"Name=(.+?);"捕获到的内容，若需要修改，可以修改此处
        if id in gene_id_list:
            new_id = id + "_" + str(gene_id_list.count(id))
        else:
            new_id = id
        gene_id_list.append(id)
        gene_id = 'gene_id "' + new_id + '";' + " "  
        
        gene_biotype = 'gene_biotype "' + genetype[re.search('gene_id "(.+?)";',str(gene_dict[i])).group(1)] + '";' + " "
        transcript_id = re.findall('transcript_id ".+?";',str(gene_dict[i]))[-1] + " " 
        if 'rna' in transcript_id:
            transcript_id = 'transcript_id "' + new_id + '_T";' + " "
        gene_pos = gene_locus[re.search('gene_id "(.+?)";',str(gene_dict[i])).group(1)]
        trans_pos = tran_locus[re.search('transcript_id "(.+?)";',str(gene_dict[i])).group(1)] if tran_locus.has_key(re.search('transcript_id "(.+?)";',str(gene_dict[i])).group(1)) else gene_pos
        fo.write("\t".join([chr,source,"gene"]) + "\t" + "\t".join(gene_pos) + "\t" + (gene_id + gene_biotype + entrez_id).strip() + "\n")  
        fo.write("\t".join([chr,source,type]) + "\t" + "\t".join(trans_pos) + "\t" + (gene_id + transcript_id + gene_biotype + entrez_id).strip() + "\n")  
        for line in gene_dict[i]:
            pre = line.split("\t")[:8] 
            if pre[2] == "exon":
                n +=1
                product = re.search('product ".+?";',line).group() + " " if re.search('product ".+?";',line) else ""
                product = change_product(product)
                exon_num = 'exon_number "%d";' %n + " "
                fo.write("\t".join(pre) + "\t" + (gene_id + transcript_id + gene_biotype + product + entrez_id + exon_num).strip() + "\n")
            else:
                product = re.search('product ".+?";',line).group() + " " if re.search('product ".+?";',line) else ""
                product = change_product(product)
                protein_id = re.search('protein_id ".+?";',line).group() + " " if re.search('protein_id ".+?";',line) else 'protein_id "' + gene_id + '_P";' + " "
                fo.write("\t".join(pre) + "\t" + (gene_id + transcript_id + gene_biotype + protein_id + product + entrez_id).strip() + "\n")
    if len(trans) > 1:               # 一个基因多个转录本情况
        chr = gene_dict[i][0].split("\t")[0]
        source = gene_dict[i][0].split("\t")[1]
        gene_pos = gene_locus[re.search('gene_id "(.+?)";',str(gene_dict[i])).group(1)]
        
        id = re.search('gene_name "(.+?)";',str(gene_dict[i])).group(1)  ### 以gene_name作为gene_id,即原gff文件中的"gene=(.+?);"或"Name=(.+?);"捕获到的内容，若需要修改，可以修改此处
        if id in gene_id_list:
            new_id = id + "_" + str(gene_id_list.count(id))
        else:
            new_id = id
        gene_id_list.append(id)
        gene_id = 'gene_id "' + new_id + '";' + " " 
        
        gene_biotype = 'gene_biotype "' + genetype[re.search('gene_id "(.+?)";',str(gene_dict[i])).group(1)] + '";' + " "
        entrez_id = 'Entrez_id "' + re.search('GeneID:(.+?)\W',str(gene_dict[i])).group(1) + '";' + " " if re.search('GeneID:(.+?)\W',str(gene_dict[i])) else ""
        fo.write("\t".join([chr,source,"gene"]) + "\t" + "\t".join(gene_pos) + "\t" + (gene_id + gene_biotype + entrez_id).strip() + "\n")  
        m = 0
        for rna in trans:     ### m记录转录本数目
            rna_list = [t for t in gene_dict[i] if rna in t]
            rna_list = sorted(rna_list,key=lambda x:x.split("\t")[3]) if rna_list[0].split("\t")[6] == "+" else sorted(rna_list,key=lambda x:x.split("\t")[3],reverse = True)
            n = 0      ### n记录exon数目
            type = re.findall('gbkey "(.+?)";',str(gene_dict[i]))
            type = "mRNA" if "mRNA" in type or "CDS" in type else genetype[re.search('gene_id "(.+?)";',str(rna_list)).group(1)]
            transcript_id = re.findall('transcript_id ".+?";',str(rna_list))[-1] + " " 
            if 'rna' in transcript_id:
                m += 1
                transcript_id = 'transcript_id "' + new_id + '_T%d";'%m + " "
            #fo.write("\t".join([chr,source,"gene"]) + "\t" + "\t".join(gene_pos) + "\t" + (gene_id + gene_biotype).strip() + "\n")  
            trans_pos = tran_locus[re.search('transcript_id "(.+?)";',str(rna_list)).group(1)] if tran_locus.has_key(re.search('transcript_id "(.+?)";',str(rna_list)).group(1)) else gene_pos
            fo.write("\t".join([chr,source,type]) + "\t" + "\t".join(trans_pos) + "\t" + (gene_id + transcript_id + gene_biotype + entrez_id).strip() + "\n")
            for line in rna_list:
                pre = line.split("\t")[:8] 
                if pre[2] == "exon":
                    n +=1
                    product = re.search('product ".+?";',line).group() + " " if re.search('product ".+?";',line) else ""
                    product = change_product(product)
                    exon_num = 'exon_number "%d";' %n + " "
                    fo.write("\t".join(pre) + "\t" + (gene_id + transcript_id + gene_biotype + product + entrez_id + exon_num).strip() + "\n")
                else:
                    product = re.search('product ".+?";',line).group() + " " if re.search('product ".+?";',line) else ""
                    product = change_product(product)
                    protein_id = re.search('protein_id ".+?";',line).group() + " " if re.search('protein_id ".+?";',line) else 'protein_id "' + gene_id + '_P";' + " "
                    fo.write("\t".join(pre) + "\t" + (gene_id + transcript_id + gene_biotype + protein_id + product + entrez_id).strip() + "\n")
    
    trans_gene = [t for t in re.findall('transcript_id ".+?";',str(gene_dict[i])) if "gene" in t]  ### 某个基因没有exon记录的情况(例如只有CDS记录)
    if len(trans_gene) > 0:
        n = 0
        chr = gene_dict[i][0].split("\t")[0]
        source = gene_dict[i][0].split("\t")[1]
        gene_pos = gene_locus[re.search('gene_id "(.+?)";',str(gene_dict[i])).group(1)]

        id = re.search('gene_name "(.+?)";',str(gene_dict[i])).group(1)  ### 以gene_name作为gene_id,即原gff文件中的"gene=(.+?);"或"Name=(.+?);"捕获到的内容，若需要修改，可以修改此处
        if id in gene_id_list:
            new_id = id + "_" + str(gene_id_list.count(id))
        else:
            new_id = id
        gene_id_list.append(id)
        gene_id = 'gene_id "' + new_id + '";' + " " 
        
        gene_biotype = 'gene_biotype "' + genetype[re.search('gene_id "(.+?)";',str(gene_dict[i])).group(1)] + '";' + " "
        entrez_id = 'Entrez_id "' + re.search('GeneID:(.+?)\W',str(gene_dict[i])).group(1) + '";' + " " if re.search('GeneID:(.+?)\W',str(gene_dict[i])) else ""
        transcript_id =  'transcript_id "' + new_id + '_T";' + " "
        fo.write("\t".join([chr,source,"gene"]) + "\t" + "\t".join(gene_pos) + "\t" + (gene_id + gene_biotype + entrez_id).strip() + "\n")
        fo.write("\t".join([chr,source,"mRNA"]) + "\t" + "\t".join(gene_pos) + "\t" + (gene_id + transcript_id + gene_biotype + entrez_id).strip() + "\n")  
        for cds in trans_gene:
            cds_list = [c for c in gene_dict[i] if cds in c]
            for line in cds_list:
                precds = line.split("\t")[:8]
                prexon = precds[:2] + ["exon",]  + precds[3:]
                protein_id = re.search('protein_id ".+?";',line).group() + " " if re.search('protein_id ".+?";',line) else 'protein_id "' + gene_id + '_P";' + " "
                n+=1
                product = re.search('product ".+?";',line).group() + " " if re.search('product ".+?";',line) else ""
                product = change_product(product)
                exon_num = 'exon_number "%d";' %n + " "
                fo.write("\t".join(prexon) + "\t" + (gene_id + transcript_id + gene_biotype + product + entrez_id + exon_num).strip() + "\n")
                fo.write("\t".join(precds) + "\t" + (gene_id + transcript_id + gene_biotype + product + protein_id + entrez_id).strip() + "\n")
fo.close()                   
 
