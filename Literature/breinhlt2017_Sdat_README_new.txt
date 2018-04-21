README
This DRYAD package contains files that related to:

Jesse W. Breinholt, Chandra Earl, Alan R. Lemmon, Emily Moriarty Lemmon, Lei Xiao, Akito Y. Kawahara (2017) Resolving relationships among the megadiverse butterflies and moths with a novel pipeline for Anchored Phylogenomics.

Please cite this publication, the DRYAD accession, and Genbank SRA accession number when using any of the data or the LEP1 enrichment kit.

This DRYAD data package will includes 42 items, including addition to this README file:

Breinholt_et_al_Supplementary_Figure_S1.pdf: Supplementary Figure S1 from Breinholt et al. (2017)
Breinholt_et_al_Supplementary_Figure_S2.pdf: Supplementary Figure S2 from Breinholt et al. (2017)
Breinholt_et_al_Supplementary_Figure_S3.pdf: Supplementary Figure S3 from Breinholt et al. (2017)
Breinholt_et_al_Supplementary_Figure_S4.pdf: Supplementary Figure S4 from Breinholt et al. (2017)
Breinholt_et_al_Supplementary_Figure_S5.pdf: Supplementary Figure S5 from Breinholt et al. (2017)
Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx: Microsoft excel document including Supplementary Table S1-S11 from Breinholt et al. (2017)
Breinholt_et_al_Supplementary_File_2_Lep1.probe_set: Specification file for the Lep1 prob set used to order probes from Agilent Technologies (http://www.agilent.com/)
Breinholt_et_al_Supplementary_File_3.docx: Word document that expands discussion of  Breinholt et al. (2017) and discusses Lepidopteran relationships in more details
Lep1_ref.tar.gz: data from the references taxa for each loci to make the probes as well as used in the IBA assembly

JAVA_SourceCode.tar.gz: compressed directory holding A.R.L (alemmon@evotutor.org) java source code
	This directory contains readme and instructions for use and to compile the java code for IdentifySpacedKmers7, QuickScan5, and ShallowMapper4. It also contain the Lep1_ProbeDesign directory used with the java programs to design the Lep1 probe set
	(IdentifySpacedKmers7, IdentifySpacedKmers7_readme.txt,Lep1_ProbeDesign,LepRefFiles.txt,QuickScan5_readme.txt,QuickScan5.java,ShallowMapper4_readme.txt,ShallowMapper4.java)
	ShallowMapper4: java script by A.R.L used to identify intron boundaries in genes for five reference taxa by mapping raw genomic reads to the corresponding transcriptomic sequences 
	QuickScan5: java script by A.R.L used to scan the additional 23 transcriptomes and ESTs by generating reference kmers using the 5-species alignments and using those kmers to map contig sequences from the transcriptomes to the candidate locus set 

Breinholt_et_al_LOG_COMMANDS.log: Set of commands used to run the bioinformatic pipeline to generate data for Breinholt et al. 2017
Scripts_README.txt: Description of the python scripts and direction how to run them.

IBA.py: python script to assemble AHE data loci by loci
IBA_trans.py: python script to assemble AHE data loci by loci
extract_probe_region.py: python script to split alignments into head, probe, and tail regions base on the beginning and end of a reference sequence in the alignment for a list of single line alignment fasta files.
s_hit_checker.py: python script to process the output of BLAST to find sequences that fit the single hit criteria
ortholog_filter.py: python script to process the output of BLAST to find if the location of the best hit on the genome is the same location as the probe target from that genome.
split.py: python script to split a single line fasta file with many loci in to locus specific fasta files
alignment_DE_trim.py: python script to trim alignments by density (1-100) and entropy (0-2)
flank_dropper.py: python script to remove poorly aligned sequences in the flanking head and tail regions
counting_monster.py:python script to count the loci per taxa and put into a tab separated matrix
removelist.py: python script to remove list of sequences from a fasta file
getlist.py: python script to get list of sequences from a fasta file
contamination_filter.py: python script to process blast results of blasting sequences from each loci against themselves using usearch to identify contamination.
remove_duplicates.py: python script to identify and remove sequences for each taxon that had more than one sequence per locus. 

taxa_list.txt: List of Sample ID's used in nexus files, Raw AHE data file name and corresponding species names in tab-delimited text
Breinholtetal_RAWDATA.tar.gz: compressed file containing the raw Illumina (2X100) AHE data
final_soap_FG120036B.fa: Assembly of Apatelodes pithala from Genbank SRA accession #SRR1794032, using multiple kmers (13,23,33,43,63) with SOAPdenovo-Trans v1.01. Different Kmer assemblies were combined with cd-hit-est and processed with the fastx toolkit. See Breinholt et al.  (2017) for more details.
final_soap_calo2.fa: Assembly of Caloptilia triadicae from Genbank SRA accession #SRR1794032, using multiple kmers (13,23,33,43,63) with SOAPdenovo-Trans v1.01. Different Kmer assemblies were combined with cd-hit-est and processed with the fastx toolkit. See Breinholt et al.  (2017) for more details.
final_soap_GV120010B.fa: Assembly of Urbanus proteus from Genbank SRA accession #SRR1794082 , using multiple kmers (13,23,33,43,63) with SOAPdenovo-Trans v1.01. Different Kmer assemblies were combined with cd-hit-est and processed with the fastx toolkit. See Breinholt et al.  (2017) for more details.
Breinholt_et_al_acrossLep_full_assemblies_all_loci.fa: Fasta formatted sequence file containing sequences that pass pipeline step 1-6 for all loci and taxa in dataset 1-3. This file can be split using the split.py to separate into fasta files of individual loci.
Breinholt_et_al_shallow_full_assemblies_all_loci.fa: Fasta formatted sequence file containing sequences that pass pipeline step 1-6 for all loci and taxa in dataset 4-6. This file can be split using the split.py to separate into fasta files of individual loci.
Breinholt_et_al_allcodonpostion123_acrossLep.nex: Nexus file containing codon position 1 & 2  & 3 for 557 loci and 75 taxa used to make dataset 1-3. See taxa_list.txt for species names of each taxon, this is a nucleotide nexus file with a CHARSET that defines each gene that starts with codon position 1. For further information see Breinholt et al. (2017) and Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx in this Dryad package for more details.
Breinholt_et_al_degen12_DS1.nex: Dataset 1 (acrossLEP_AHE). Nexus file containing codon position 1 & 2 for 557 loci and 23 taxa. See taxa_list.txt for species names of each taxon, this is a nucleotide nexus file with a CHARSET that defines each gene that starts with codon position 1. Synonymous signal was removed using degen v1.4 Perl script (http://www.phylotools.com), and the third codon has been removed. Loci names correspond to Loci numbers in the Lep1 enrichment kit included in this DRAYD package. For further information see Breinholt et al. (2017) and Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx in this Dryad package for more details.
Breinholt_et_al_aminoacid_DS1.nex: Dataset 1 (acrossLEP_AHE). Nexus file containing amino acid data for 557 loci and 23 taxa. See taxa_list.txt for species names of each taxon, this is an amino acid nexus file with a CHARSET that defines each loci. Loci names correspond to Loci numbers in the Lep1 enrichment kit included in this DRAYD package.
Breinholt_et_al_degen12_DS2.nex: Dataset 2 (acrossLEP_AHE+PARTtrans). Nexus file containing codon position 1 & 2 for 557 loci and 75 taxa. See taxa_list.txt for species names of each taxon, this is a nucleotide nexus file with a CHARSET that defines each gene that starts with codon position 1. Synonymous signal was removed using degen v1.4 Perl script (http://www.phylotools.com), and the third codon has been removed. Loci names correspond to Loci numbers in the Lep1 enrichment kit included in this DRAYD package. For further information see Breinholt et al. (2017) and Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx in this Dryad package for more details.
Breinholt_et_al_aminoacid_DS2.nex: Dataset 2 (acrossLEP_AHE). Nexus file containing amino acid data for 557 loci and 75 taxa. See taxa_list.txt for species names of each taxon, this is an amino acid nexus file with a CHARSET that defines each loci. Loci names correspond to Loci numbers in the Lep1 enrichment kit included in this DRAYD package. For further information see Breinholt et al. (2017) and Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx in this Dryad package for more details.
Breinholt_et_al_degen12_DS3.nex: Dataset 3 (acrossLEP_AHE+ALLtrans ). Nexus file consists of both AHE and the transcriptomic data of Kawahara and Breinholt 2015. The file contains codon position 1 & 2 for 2948 loci and 76 taxa. See taxa_list.txt for species names of each taxon, this is a nucleotide nexus file with a CHARSET that defines each gene that starts with codon position 1. Synonymous signal was removed using degen v1.4 Perl script (http://www.phylotools.com), and the third codon has been removed. Loci names correspond to Loci numbers in the Lep1 enrichment kit included in this DRAYD package. For further information see Breinholt et al. (2017) and Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx in this Dryad package for more details.
Breinholt_et_al_DS4.nex: Dataset 4 (shallow_probe+flanks). Nexus file containing 749 loci and 48 taxa. Alignments were trimmed with a density of 60% and entropy of 1.5 using alignment_DE_trim.py and flacking regions were processed with the flank_dropper.py to remove head or tail sequences using 2 standard deviations for both the head and tail. See taxa_list.txt for species names of each taxon, this is a nucleotide nexus file with a CHARSET that defines each gene. For further information see Breinholt et al. (2017) and Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx in this Dryad package for more details.
Breinholt_et_al_DS5.nex: Dataset 5 (shallow_probe). Nexus file containing 749 loci and 48 taxa. The Extract_probe_region.py script was used on Dataset 4 to isolate data coming from the probe region. See taxa_list.txt for species names of each taxon, this is a nucleotide nexus file with a CHARSET that defines each gene. For further information see Breinholt et al. (2017) and Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx in this Dryad package for more details.
Breinholt_et_al_DS6.nex: Dataset 6 (shallow_flanks). Nexus file containing 749 loci and 35 taxa. The Extract_probe_region.py script was used on Dataset 4 to isolate data coming from the flanking regions region. See taxa_list.txt for species names of each taxon, this is a nucleotide nexus file with a CHARSET that defines each gene. For further information see Breinholt et al. (2017) and Breinholt_et_al_Supplementary_File_1_S1-S11.xlsx in this Dryad package for more details.

For questions regarding these data, contact Jesse Breinholt at jessebreinholt@gmail.com or Akito Kawahara at kawahara@flmnh.ufl.edu. If your question is about the Java scripts in JAVA_SourceCode.tar.gz contact Alan Lemmon (alemmon@evotutor.org) 
