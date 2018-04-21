#First call perl script which reads trees and writes 
#makeRtrees.pl must be in path!
command -v /home/galaxy/galaxy-dist/tools/osiris_phylogenetics/phylographics/makeRtrees.pl >/dev/null 2>&1 || { echo > $6 "ERROR: You must change the absolute path in the tab2trees.sh file to match your local system."; exit 1; }


#$1 infile
#$2 outfile
#$3 tree type (ie phylogram)
#$4 yes|no exclude tips
#$5 yes|no label taxa
#$6 name of Rfile "Rfile" by default in xml
#$7 yes|no to label OTUs with QUERY in title
#$8 yes|no to conduct midpoint rooting

/home/galaxy/galaxy-dist/tools/osiris_phylogenetics/phylographics/makeRtrees.pl $1 $2 $3 $4 $5 $7 $8 > $6 2>log.txt

R --vanilla < $6 2>log.txt
