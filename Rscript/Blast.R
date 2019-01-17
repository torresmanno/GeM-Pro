args = commandArgs(trailingOnly = T)

main.strain <- args[1]
strains <- readLines(args[2])
subset.genes <- args[3]
fasta.path <- args[4]
outpath <- args[5]

pattern <- "_protein.faa"
blast <- function(main.strain, strains, pattern, fasta.path, output){
  cmd.blast <- function(){
    paste0("blastp -query ", fasta.path, main.strain, pattern," -evalue 1e-5 -subject ",fasta.path, strains, pattern," -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp' -max_target_seqs 1 -out ", output, main.strain, "_", strains, ".tsv")
  }
  cmds <- cmd.blast()
  dir.create(output, showWarnings = F)

null <- sapply(cmds, system)
}

if (subset.genes == "ALL"){

  blast(main.strain, strains, pattern, fasta.path, outpath)

} else if(file.exists(subset.genes)){

  genes.table <- read.csv(subset.genes, header = T, stringsAsFactors = F)
  genes <- genes.table$Acc_Ver
  if (!require("seqinr")) install.packages("seqinr")
  fasta <- seqinr::read.fasta(paste0(fasta.path, main.strain, pattern), seqtype = "AA", as.string = T)
  sub.fasta <- fasta[genes]
  write.fasta(sub.fasta, names=gsub(">", "", getAnnot(sub.fasta)), file.out = "sub_fasta.faa")

  cmd.blast <- function(){
    paste0("blastp -query sub_fasta.faa -evalue 1e-5 -subject ",fasta.path, strains, pattern," -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp' -max_target_seqs 1 -out ", outpath, main.strain, "_", strains, ".tsv")
  }
  cmds <- cmd.blast()
  dir.create(outpath, showWarnings = F)

  null <- sapply(cmds, system)

}

