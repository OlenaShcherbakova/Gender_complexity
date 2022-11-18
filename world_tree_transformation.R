#Transforming newick tree into nexus tree

wtree <- read.newick("data/phylogenies/world/world.tre")

for(i in 1:length(wtree$tip.label)){
  wtree$tip.label[i] <- sub(pattern = '^[a-zA-Z_]*\\.[a-zA-Z_]*\\.', replacement = "", x = wtree$tip.label[i])
}

write.nexus(wtree, file = "data/phylogenies/world/posterior.trees")