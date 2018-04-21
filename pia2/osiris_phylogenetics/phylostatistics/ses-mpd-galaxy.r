#!/sw/math/R-2.15.3-shlib/bin/Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOTE: since picante is licensed under the GPL, and this program relies on 
# picante, the program is licensed under the GPL regardless
#
# See: http://cran.r-project.org/web/packages/picante/picante.pdf,
# http://www.gnu.org/licenses/old-licenses/gpl-2.0-faq.html#IfLibraryIsGPL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('picante') 

args <- commandArgs(trailingOnly = TRUE)

sample <- read.table(file = args[1])
tree <- read.tree(file = args[2])

# get community data matrix of sample
comm <- sample2matrix(sample)
# get phylogenetic distance matrix of tree
phydist <- cophenetic(tree)

# finally, run the processed info through ses.mpd to get the result we want
result <- ses.mpd(comm, phydist)

# capture result and output to file
out <- capture.output(result)
cat(out, file = args[3], sep = "\n")
