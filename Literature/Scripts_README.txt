IBA.py: python script to assemble AHE data loci by loci
	Use ./IBA.py -h to get a list of options. The script requires Concurrent Futures Python module for Python v. 2.5 - 2.7.9 (https://pypi.python.org/pypi/futures), Usearch (http://www.drive5.com/usearch/), and Bridger (https://sourceforge.net/projects/rnaseqassembly/). The program can use multiple threads but all threads must be on the same node. You must use the combination of threads and processes to equal the number of available threads. Each process will use the number of threads given so to get the total multiply -p by -t to get the total thread use. If only 8 threads are available you could use (-t 1 -p 8) or (-t 2 -p 4). We generally used nodes with 32 threads and use either (-p 16 -t 2) or (-t 1 -p 32). The more processes (-p) will mean more loci are being concurrently assembled and it is generally faster to set –p as high as possible and use one thread per process. In some instances when there is a lot of on target data Bridger can run out of memory with only one thread and the threads need to be increased (-t 2). Note: All reads must be larger than the Kmer used in the assembly.

	-raw1	first raw data file
	-raw2	second raw data file
	-d	directory where reference fasta files are located (a directory where these is a single file per loci you want to assemble)
	-n	the number of iterations
	-t	Number of threads for usearch and Bridger
	-p	Number of loci to concurrently assemble
	-g	Pair gap length bridger input (see bridger readme file --pair_gap_length is gap length of paired reads)
	-c	kmer coverage depth for final assembly
	-taxa	taxa name of the raw data used
	-label	label of reference taxa, must be in every reference fasta file
	-k	kmer length for Bridger assembly 19-32**

example: ./IBA.py -raw1 taxa1_read1.fastq -raw2 taxa1_read1.fastq -d ./ref -n 3 -t 2 -p 4 -g 200 -c 25 -taxa IDnumber_Family_Genus_species -label BMORI -K 25


IBA_trans.py: python script to assemble AHE data loci by loci
	Use ./IBA.py -h to get a list of options. The script requires Concurrent Futures Python module for Python v. 2.5 - 2.7.9 (https://pypi.python.org/pypi/futures), Usearch (http://www.drive5.com/usearch/), and Bridger (https://sourceforge.net/projects/rnaseqassembly/). The program can use multiple threads but all threads must be on the same node. You must use the combination of threads and processes to equal the number of available threads. Each process will use the number of threads given so to get the total multiply -p by -t to get the total thread use. If only 8 threads are available you could use (-t 1 -p 8) or (-t 2 -p 4). We generally used nodes with 32 threads and use either (-p 16 -t 2) or (-t 1 -p 32). The more processes (-p) will mean more loci are being concurrently assembled and it is generally faster to set –p as high as possible and use one thread per process. In some instances when there is a lot of on target data Bridger can run out of memory with only one thread and the threads need to be increased (-t 2).Note: All reads must be larger than the Kmer used in the assembly.


	-raw1	first raw data file
	-raw2	second raw data file
	-d	directory where reference fasta files are located (a directory where these is a single file per loci you want to assemble)
	-n	the number of iterations
	-t	Number of threads for usearch and Bridger
	-p	Number of loci to concurrently assemble
	-g	Pair gap length bridger input (see bridger readme file --pair_gap_length is gap length of paired reads)
	-c	kmer coverage depth for final assembly
	-taxa	taxa name of the raw data used
	-label	label of reference taxa, must be in every reference fasta file
	-k	kmer length for Bridger assembly 19-32**

example: ./IBA.py -raw1 taxa1_read1.fastq -raw2 taxa1_read1.fastq -d ./ref -n 3 -t 2 -p 4 -g 200 -c 25 -taxa IDnumber_Family_Genus_species -label BMORI -K 25

extract_probe_region.py: python script to split alignments into head, probe, and tail regions base on the beginning and end of a reference sequence in the alignment for a list of single line alignment fasta files.
	To use the script the aligned fasta file must contain sequences data on a single line. The script will not work on a multiple line alignment fasta file. You must give the script a list of files you want to process, name of the reference sequence you want to use to trim the alignment, and a out directory (./BMORI.py list.txt refname outdir). The list must contain the name of a single fasta file per line and the files must be located in the directory where you submit the script unless the full path to the alignments are provided in the list. The name of the reference sequences will define the probe region with the first base (ATGC) and last BASE (ATGC) in the alignment. 
	example: ./extract_probe_region.py inlist BMORI outdir
	
s_hit_checker.py: python script to process the output of BLAST to find sequences that fit the single hit criteria.
	This script process a blast formatted table in blast -outfmt 6 from a blast analysis to find the best hits on 3 different contains and allows 3 of the best hits per contig (for example blast output from this: blastn -task blastn -query infile -db genomedatabasename -out outfile -outfmt 6 -max_target_seqs 3 -max_hsps 3). To run you need the blast table and a decimal values that is the percent of the second best hit bitscore/ best hit bitscore that you consider to close to determine if the sequence is single copy. For example if the best hit bit score is 100 and the second best hit is 90 then the second best hit bit score is >= to 90 and is too close to confidently consider it as single hit.  We have found 90 to be a good indication that the sequences is single copy but setting this lower will increases the probability the sequences is single copy.
	example: ./s_hit_checker.py inblasttable .90

ortholog_filter.py: python script to process the output of BLAST to find if the location of the best hit on the genome is the same location as the probe target from that genome.
	This script process a blast formatted table in blast -outfmt 6 from a blast analysis to find the single best hit location of a sequence (for example blast output from this: blastn -task blastn -query infile -db genomedatabasename -out outfile -outfmt 6 -max_target_seqs 1 -max_hsps 1). To run you need the blast table and the name of reference taxa that must be in the name of sequences of each loci from the reference genome. The fasta file you use to blast must contain sequence data from the reference genome for each loci you what to test. It will then compare the location where the data from the reference hit its genome to where the others sequences hit the genome to check for orthology. 
	example: ./ortholog_filter.py blasttable BMORI
	
split.py: python script to split a single line fasta file with many loci in to locus specific fasta files
	the script splits a fasta file where data for each sequences is on a single line from multiple loci in a file for each locus. The script will not work on a multiple line fasta file nor will it work on a file with blank lines. The script needs a single line formatted fasta file containing sequences from many loci and taxa and a character as a delimiter to split the loci from the taxa name. I suggest you uses an underscore. The name of the loci must be in the first part of the sequences name then the delimiter then what ever else you want (do not use dashes"-" or spaces in sequences names) for example example L001_Bombyx_mori. It works with an underscore as a delimiter. Other delimiters have not been tested but should work. example: ./split.py infile _
	example: ./split.py infile _

alignment_DE_trim.py: python script to trim alignments by density (1-100) and entropy (0-2)
	To use the script the processed fasta file must contain sequences data on a single line. The scripts takes infile, outfile, density, entropy, 1 (print pdf graph) or 0 (do not print pdf), optional taxa name that spans probe region. Density is per column in the alignment and equals the number of valid bases (ATGC) divided by the number of taxa so 75 would mean 75% percent of the taxa have data for that column in the alignment. Entropy is based on nucleotide entropy that ranges from 0-2, estimated with equation 1 of Xia et al. (2003). If the name of sequences for the probe region is given the first base (ATGC) and last base (ATGC) for the sequences is used to define the probe region. Within the probe region no trimming will be done based on density or entropy and the location of the probe region will be plotted on the graph. 
	example: python3 alignment_DE_trim.py infile outfile 75 1.75 1 BMORI

flank_dropper.py: python script to remove poorly aligned sequences in the flanking head and tail regions
	To use the script the processed fasta file must contain sequences data on a single line and all sequence data must be upper case and requires biopython. To turn sequence data to upper case you can use this sed command in unix (sed  -i '/^>/! y/acgtn/ACGTN/' filename, or for multiple files that end in fasta like this sed  -i '/^>/! y/acgtn/ACGTN/' *.fasta). The script takes an infile, outfile, taxa name for the probe region, desired head flanking standard deviation, and desired tail flanking standard deviation (./flank_dropper.py INfile OUTfile taxaNameForProbe StdevHead StdevTail).
	example: ./flank_dropper.py INfile OUTfile BMORI 2 2

counting_monster.py:python script to count the loci per taxa and put into a tab separated matrix
	To use the script the processed fasta file must contain sequences data on a single line. The script will not work on a multiple line fasta file. . The script needs a single line formatted fasta file containing sequences from many loci and taxa and a character as a delimiter to split the loci from the taxa name. I suggest you uses an underscore. The name of the loci must be in the first part of the sequences name then the delimiter then what ever else you want (do not use dashes"-" or spaces in sequences names) for example example L001_Bombyx_mori. It works with an underscore as a delimiter. Other delimiters have not been tested but should work.The script also expects each sequences to end in either _R being a reference sequences or contain the word comp which is added in the IBA assembly for each sequences.
	example: ./counting_monster.py fastafile _ 

removelist.py: python script to remove list of sequences from a fasta file
	To use you must have a fasta file and a list of sequences you want to remove from the fasta file and requires biopython. Sequence names must match 100%. Put the name of each desired sequence in a text file with one name per line. 
	example:  ./removelist.py fastafile wantedseqnamelist outputfile

getlist.py: python script to get list of sequences from a fasta file
	To use you must have a fasta file and a list of sequences you want from the fasta file and requires biopython. Sequence names must match 100%. Put the name of each desired sequence in a text file with one name per line. 
	example:  ./getlist.py fastafile wantedseqnamelist outputfile
	
contamination_filter.py: python script to process blast results of blasting sequences from each loci against themselves using usearch to identify contamination.
	This script process the results of selfblasting loci using usearch. The scripts needs BLAST TABLE (the blast output for all loci, blastfile must have an extention like .txt) _(delimitation character) 3(postion taxonomy you want to check for contamination example loci_DJM10A349_Nymphalidae_Limenitidinae_Limenitidiini_Cymothoe_lurida [loci_Code_Family_subfamily_tribe_genus_species, 0_1_2_3_4_5_6] subfamily(name of the taxonomic postion). 
	You can check any taxonomic postion it will tell you what sequences are 99% identical from different taxonomic ranks such as families or subfamilies.
	I suggest to look at the output pairfeq_list.txt to see if you have any taxa from different taxonomic groups with high amounts of identical sequences.
	Also look at the the output Locifeq_list.txt if there are a lot of sequnces for a certain loci it is probably just not variable between the chosen taxonomic group
	I suggest you alter the _del_list.txt by removing sequences that may have been identified due to little variation in the loci or if you identify a taxa pair with high contamination. Any seq name in the _del_list.txt will be deleted so here is your chance to decide if you really want to delete it.
	Example: ./contamination_filter.py BLASTTABLE.txt _ 3 subfamily
	
remove_duplicates.py: python script to identify and remove sequences for each taxon that had more than one sequence per locus. 
	This script processes a fasta file and deletes sequences when there remains two copies for the same loci and taxon. It requires Biopython.
	example: ./remove_duplicates.py fastafile
