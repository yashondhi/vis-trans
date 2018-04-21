#!/usr/bin/env python

#Script written by Jesse Breinholt (jessebreinholt@gmail.com) and Chandra Earl (sunray1@ufl.edu)
#Python script to assemble AHE data loci by loci. Use ./IBA_trans.py -h to get a list of options. It has been modified to assemble transciptoßmes with single end illumina data
#The script requires Concurrent Futures Python module for Python v. 2.5 - 2.7.9 (https://pypi.python.org/pypi/futures),
#Usearch (http://www.drive5.com/usearch/), and Bridger (https://sourceforge.net/projects/rnaseqassembly/).
#The program can use multiple threads but all threads must be on the same node. You must use the combination
#of threads and processes to equal the number of available threads. Each process will use the number of threads
#given so to get the total multiply -p by -t to get the total thread use. If only 8 threads are available you
#could use (-t 1 -p 8) or (-t 2 -p 4). We generally used nodes with 32 threads and use either (-p 16 -t 2) or (-t 1 -p 32).
#The more processes (-p) will mean more loci are being concurrently assembled and it is generally faster to set -p as
#high as possible and use one thread per process. In some instances when there is a lot of on target data Bridger
#can run out of memory with only one thread and the threads need to be increased (-t 2).

import os, sys
import concurrent.futures

#Define shell commands
blast_command_ntaa = "usearch -ublast %s -db %s -evalue 1e-5 -userout %s -userfields query -threads %s -strand both -maxaccepts 1"
blast_command_ntnt = "usearch -ublast %s -db %s -evalue 1e-15 -userout %s -userfields query -threads %s -strand both -maxaccepts 1"
blast_command_final = "usearch -ublast %s -db %s -evalue 1e-5 -target_cov 0.80 -qsegout %s -blast6out %s -threads %s -strand both"
get_single = "usearch -fastx_getseq %s -label %s -label_substr_match -fastaout %s"
cat_command = "cat %s %s > %s"
getseqs_command = "usearch -fastx_getseqs %s -labels %s -fastqout %s"
bridger_command_cov = "Bridger.pl --seqType fq --single %s --CPU %s --min_kmer_coverage %s --output %s --kmer_length %s"
bridger_command = "Bridger.pl --seqType fq --single %s --CPU %s --output %s --kmer_length %s"
blast_filter = "usearch -ublast %s -db %s -evalue 1e-5 -target_cov 0.80 -threads %s -matched %s -strand both"
remove_command = "rm -r %s"

arguments = sys.argv

namelist = []
name = []

if len(arguments) == 1:
	arguments.append("-h")
try:
	hflag = arguments.index("-h")
except:
	hflag = None

if hflag:
	sys.exit("\n\n#################################################################\n\n\t\tIteritive Baited Assembly IBA\n\n\n\t-raw1\tfirst raw data file\n\t-d\tdirectory where fasta files are located\n\t-n\tthe number of iterations\n\t-t\tNumber of threads for usearch and Bridger\n\t-p\tNumber of loci to concurrently assemble\n\t-g\tPair gap length bridger input\n\t-c\tkmer coverage depth for final assembly\n\t-taxa\ttaxa name of the raw data used\n\t-label\tlabel to search for when pulling out single sequences\n\t-k\tkmer length for Bridger assembly 19-32**\n\t  **note all reads must be bigger then the set kmer\n#################################################################\n\n")

try:
	raw1flag = arguments.index("-raw1")
	raw1_file = arguments[raw1flag+1]
except:
	sys.exit("Error: need raw data")


try:
	dflag = arguments.index("-d")
	fastadir = arguments[dflag+1]
except:
	sys.exit("Error: need directory where fasta files are located")

try:
	nflag = arguments.index("-n")
	iterations = int(arguments[nflag+1])
except:
	iterations = 3
try:
	taxaflag = arguments.index("-taxa")
	taxaname = arguments[taxaflag+1]
except:
	sys.exit("Error: need taxa name")
try:
	labelflag = arguments.index("-label")
	labelname = arguments[labelflag+1]
except:
	labelname = "BMORI"

try:
	cflag = arguments.index("-c")
	coverage = int(arguments[cflag+1])
except:
	coverage = 2
try:
	tflag = arguments.index("-t")
	Threads = arguments[tflag+1]
except:
	Threads = 4
try:
	pflag = arguments.index("-p")
	processes = arguments[pflag+1]
except:
	processes = 1
try:
	gflag = arguments.index("-g")
	gaplength = arguments[gflag+1]
except:
	gaplength = 200
try:
	kflag = arguments.index("-k")
	kmer = arguments[kflag+1]
except:
	kmer = 25


#define functions - nt/nt
def blast_iteration(rawdata1, bridgerin, blastout, listname, getseqname, T):
	#Blast
	os.system(blast_command_ntnt % (rawdata1, bridgerin, blastout, T))
	#get seqs
	os.system(getseqs_command % (rawdata1, blastout, getseqname))
			
def IAB(i):
#Blast data against alignments- nt/aa
	#blast
	print("\n\n\n\n\n########Assembling " + name[i] + "########\n\n\n")
	os.system(blast_command_ntaa % (raw1_file, fastadir + "/" + fastalist[i], name[i] + "_B1_r1.m8", Threads))
	#cat together
	os.system(cat_command % (name[i] + "_B1_r1.m8", name[i] + "_B2_r1.m8", name[i] + "_list_r1.txt"))
	#get seqs
	os.system(getseqs_command % (raw1_file, name[i] + "_B1_r1.m8", name[i] + "_B1_r1.tfastq"))
	#bridger
	os.system(bridger_command % (name[i] + "_B1_r1.tfastq", Threads, name[i] + "_1_out", kmer))
	killfq1=name[i] + "_B1_r1.tfastq"
	os.system(remove_command % (killfq1))
	
	os.system(get_single % (fastadir + "/" + fastalist[i], labelname, name[i] + "_singleref.fa"))
	os.system(blast_filter % (name[i] + "_1_out/Bridger.fasta", name[i] + "_singleref.fa", Threads, name[i] + "_filter.fa"))
#use function to loop through - blast assembly against raw data
	if os.path.getsize(name[i] + "_filter.fa") > 0:
		namelist.append(name[i])
		for n in range(iterations):
			if n == 0:
				blast_iteration(raw1_file, name[i] + "_filter.fa", name[i] + "_B1_r" + str(n+2) + ".m8", name[i] + "_list_r" + str(n+2) + ".txt", name[i] + "_B1_r" + str(n+2) + ".tfastq", Threads)
				os.system(bridger_command % (name[i] + "_B1_r" + str(n+2) + ".tfastq", Threads, name[i] + "_" + str(n+2) + "_out", kmer))
				killfq1=name[i] + "_B1_r" + str(n+2) + ".tfastq"
				os.system(remove_command % (killfq1))
				print("\n\n\nIteration " + str(n+2) + " done\n\n\n")
			elif n == iterations-1:
				blast_iteration(raw1_file, name[i] + "_filter.fa", name[i] + "_B1_r" + str(n+2) + ".m8", name[i] + "_list_r" + str(n+2) + ".txt", name[i] + "_B1_r" + str(n+2) + ".tfastq", Threads)
				os.system(bridger_command_cov % (name[i] + "_B1_r" + str(n+2) + ".tfastq", Threads, coverage, name[i] + "_" + str(n+2) + "_out", kmer))
				killfq1=name[i] + "_B1_r" + str(n+2) + ".tfastq"
				os.system(remove_command % (killfq1))
				print("\n\n\nFinal Iteration " + str(n+2) + " done\n\n\n")
			else:
				blast_iteration(raw1_file, name[i] + "_" + str(n+1) + "_out" + "/Bridger.fasta", name[i] + "_B1_r" + str(n+2) + ".m8", name[i] + "_list_r" + str(n+2) + ".txt", name[i] + "_B1_r" + str(n+2) + ".tfastq", Threads)
				os.system(bridger_command % (name[i] + "_B1_r" + str(n+2) + ".tfastq", Threads, name[i] + "_" + str(n+2) + "_out", kmer))
				killfq1=name[i] + "_B1_r" + str(n+2) + ".tfastq"
				os.system(remove_command % (killfq1))
				print("\n\n\nIteration " + str(n+2) + " done\n\n\n")
#final blast against alignments
		os.system(blast_command_final % (name[i] + "_" + str(n+2) + "_out" + "/Bridger.fasta", name[i] + "_singleref.fa", name[i] + "_final.fa", name[i] + "_blasttable.txt", Threads))
	
		print(name[i] + " finished\n\n")
	if os.path.getsize(name[i] + "_filter.fa") == 0:
		print("############################## Filter for " + name[i] + " was empty ######################")
		pass

#to allow multiple processors
def IBA_process():
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(processes)) as executor:
		fs = [executor.submit(IAB, i) for i in range(len(fastalist))]


fastalist = os.listdir(fastadir)
for i in range(len(fastalist)):
	try:
		prename = fastalist[i].split(".")
		name0 = prename[0]
	except:
		name0 = fastalist[i]
	name.append(name0)
print(name)
print(fastalist)
print(processes)
IBA_process()



print("##############################\n\nCompiling loci..............\n\n\n")
rename_final_out = open(taxaname +  "_finalseqs.fasta", "w")
final_table_out = open(taxaname + "_finaltable.table", "w")
#write to files
for name in namelist:	
	rename_final = open(name + "_final.fa", "r")
	line = rename_final.readline()
	while line:
		if line[0] == ">":
			line2 = line.lstrip(">").strip("\n")
			rename_final_out.write(">" + name + "_" + taxaname + "__" + line2 + "\n")
			line = rename_final.readline()
		else:
			rename_final_out.write(line)
			line = rename_final.readline()
	rename_final.close()
	table_in = open(name + "_blasttable.txt", "r")
	line = table_in.readline()
	while line:
		tablelist = line.split()
		tablestring = "\t".join(tablelist[1:])
		final_table_out.write(name + "_" + tablelist[0]+ "__" + taxaname + "\t" + tablestring + "\n")
		line = table_in.readline()
	table_in.close()
rename_final_out.close()
final_table_out.close()

#remove intermediate files
os.system(remove_command % ("*.m8"))
os.system(remove_command % ("*.txt"))
os.system(remove_command % ("*.fa"))
os.system(remove_command % ("*out"))

print("#################\nIBA Complete beep.. beep.. beep...\n\n\n\t\t ^ ^\n\t\t @ @\n\t\t  $\n\t       MMMMMMM\n\t       M <O> M\n\t       M     M\n\n\t   Have a Nice Day!\n\n")
sys.exit()
