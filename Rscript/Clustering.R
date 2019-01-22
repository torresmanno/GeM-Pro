args <- commandArgs(trailingOnly = T)

pdf("Clustering.pdf")
pan <- readRDS(args[1])
nb <- as.numeric(args[2])
if (!require("pvclust")) install.packages("pvclust")
if (!require("ape")) install.packages("ape")
Clustering <- function(pangenome, nboot) {
  pangenome[pangenome != 0] <- 1
  Cluster <- pvclust::pvclust(pangenome, method.dist = "manhattan", method.hclust = "average", nboot = nboot)
  return(Cluster)
}
clust <- Clustering(pan, nb)
plot(clust)
saveRDS(clust, "Clustering.RDS")

ape::write.tree(ape::as.phylo(clust$hclust), file = "Clustering.tre")
