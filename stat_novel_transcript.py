#!/usr/bin/python
import os
import re
import sys
usge = """
usge :
        python TO.py -dir filepath -o out.file
"""
def argvlist(argv):
        argvlist = []
        for arg in argv:
                argvlist.append(arg)
        return argvlist
def samplename(dir):
        samplelist = os.listdir(dir)
        return samplelist
def Assembled_Transcripts(samplepathname):
	filename = samplepathname
	fileopen = open(filename).readlines()
	strin = fileopen[0].split()[3]
	return strin
def Assembled_Loci(samplepath):
        filename = samplepath
        fileopen = open(filename).readlines()
        strin = fileopen[1].split()[3]
	return strin
def Known_Transcripts(samplepath):
        filename = samplepath
        fileopen = open(filename).readlines()
        strin = fileopen[2].split()[3]
	return strin 
def Corresponding_Known_Loci(samplepath):
        filename = samplepath
        fileopen = open(filename).readlines()
        strin = fileopen[3].split()[3]
	return strin
def In_Novel_Transcripts(samplepath):
        filename = samplepath
        fileopen = open(filename).readlines()
        strin = fileopen[4].split()[3]
	return strin
def Correspondin_Novel_Loci(samplepath):
        filename = samplepath
        fileopen = open(filename).readlines()
        strin = fileopen[5].split()[3]
	return strin
def machname(list):
	listall = []
	for itme in list:
		if re.match(r'.+stat$',itme):
			listall.append(itme)
	return listall
def isdir(listdir,path):
        newlistdir = []
        for dir in listdir:
                dirpath = path+"/"+dir
                if os.path.isdir(dirpath):
                        newlistdir.append(dir)
        return newlistdir
def main():
	mainlist = argvlist(sys.argv)[1:]
        if len(mainlist)<1:
                print usge
                sys.exit()
        elif mainlist[0] == "-h" or mainlist[0] == "--help":
                print usge
                sys.exit()
        elif len(mainlist)!=4:
                print usge
                sys.exit()
        else:
		samplelist = samplename(mainlist[1])
		samplelist = isdir(samplelist,mainlist[1])
		if "cuffmerge" in samplelist:
			samplelist.remove("cuffmerge")
		else:
			pass
		samplelist = sorted(samplelist)
		filesave = open(mainlist[3],"w+")
		filesave.write("Class\tTotal Assembled Transcripts\tTotal Assembled Loci\tIn Known Transcripts\tCorresponding Known Loci\tIn Novel Transcripts\tCorresponding Novel Loci\n")
		for sampledirname in samplelist:
			
			samplenamestat = mainlist[1]+"/"+sampledirname+"/transcripts.gtf.known_and_novel_stat"
			filesave.write(sampledirname)
			filesave.write("\t"+Assembled_Transcripts(samplenamestat))
			filesave.write("\t"+Assembled_Loci(samplenamestat))
			filesave.write("\t"+Known_Transcripts(samplenamestat))
			filesave.write("\t"+Corresponding_Known_Loci(samplenamestat))
			filesave.write("\t"+In_Novel_Transcripts(samplenamestat))
			filesave.write("\t"+Correspondin_Novel_Loci(samplenamestat)+"\n")
		cuffmergestatpath = mainlist[1]+"/cuffmerge/"
		cuffmergestatlist =  machname(samplename(cuffmergestatpath))
		cuffmerge1 = cuffmergestatpath+cuffmergestatlist[0]
		filesave.write("cuffmerge"+"\t")
		filesave.write(Assembled_Transcripts(cuffmerge1)+"\t")
		filesave.write(Assembled_Loci(cuffmerge1)+"\t")
		filesave.write(Known_Transcripts(cuffmerge1)+"\t")
		filesave.write(Corresponding_Known_Loci(cuffmerge1)+"\t")
		filesave.write(In_Novel_Transcripts(cuffmerge1)+"\t")
		filesave.write(Correspondin_Novel_Loci(cuffmerge1)+"\n")
		filesave.close()
if __name__ == '__main__':
        main()
