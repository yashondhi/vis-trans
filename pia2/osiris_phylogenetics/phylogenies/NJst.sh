#First call perl script which reads trees and writes 
/home/galaxy/galaxy-dist/tools/Rtools/makeNJst.pl $1 $2 > Rnjst.R 2>log.txt

R --vanilla < Rnjst.R 2>log.txt
