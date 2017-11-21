# I would use something like:
# /share/apps/bin/R --slave --no-save --no-restore --no-environ --silent --args inTreeFile outTreeFile < <code dir>/MidPoint_Rootingp.R"

library(ape)
library(phangorn)

args <- commandArgs(trailingOnly = TRUE)
inFile = args[1]
outFile = args[2]


in_tree <- read.tree(inFile)
root_tree = midpoint(in_tree)
write.tree(root_tree,file=outFile)
