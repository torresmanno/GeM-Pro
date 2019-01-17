args = commandArgs(trailingOnly = T)

ref.str.df <- read.csv(args[1], stringsAsFactors = F)
strains <- readLines(args[2])
ANIs.dir <- args[3]

ANI_FUN <- function(ANIs.dir, ref.strain.df, strList){

    # Matrix to put the ANI values
  ani.table <-matrix(nrow = length(strList), ncol =  length(ref.strain.df[,1]), dimnames = list(strList, ref.strain.df[,1]))

  # Fill the matrix with the ANI values
  for(ref in colnames(ani.table)){
    for(str in rownames(ani.table)){
      ani.calc <- readLines(paste0(ANIs.dir,ref, "_VS_", str, ".ANI"))[14]
      ani.table[str, ref] <- as.numeric(strsplit(ani.calc, split = " ")[[1]][3])
    }
  }

  # making the group by the max ANI
  ani.group <- sapply(colnames(ani.table),function(x) NULL)
  for(str in rownames(ani.table)){
    if(max(ani.table[str,])>=94){
      group <- names(which.max(ani.table[str,]))
      ani.group[[group]] <- c(ani.group[[group]], str)
    }
  }
  ani.group <- ani.group[ref.strain.df[,1]]
  # Adding the ref strain to the group
  names(ani.group) <- ref.strain.df[,2]


  return(list(ANI.groups = ani.group, ANI.table = ani.table))
}
ANI <- ANI_FUN(ANIs.dir, ref.str.df, strains)

# dir.create(out.path, showWarnings = F, recursive = T)

write.csv(as.data.frame(ANI[["ANI.table"]]), paste0( "ANI_table.csv"),row.names = T,quote = F)
ani.group <- ANI$ANI.groups
ani.group.df <- data.frame(Strains = unlist(ani.group), Groups = rep(names(ani.group), sapply(ani.group, length)))

write.table(ani.group.df, paste0( "Phylogenetic_groups.csv"),row.names= F, col.names = T, quote = F, sep = ",")
# saveRDS(ANI, paste0(out.path, "ANI.RDS"))
