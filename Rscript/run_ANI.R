args = commandArgs(trailingOnly = T)

ref.str <- args[1]
strains <- readLines(args[2])
OAT.path <- args[3]
blast.path <- args[4]
threads <- args[5]
fasta.paths <- args[6]
out.path <- args [7]

run_ANI <- function(OAT.path, blast_path, threads, ref.strain, strain, fasta.path, out.path){
  cmd.ani <- function(ref, str){
    paste0(
      "java -jar ",
      OAT.path,
      "OAT_cmd.jar -blastplus_dir ",
      blast_path,
      " -method ani -num_threads ",
      threads,
      " -fasta1 ",
      fasta.path,
      ref,
      "_genomic.fna -fasta2 ",
      fasta.path,
      str,
      "_genomic.fna > ",
      out.path,
      ref,
      "_VS_",
      str,
      ".ANI"
    )
  }
  dir.create(out.path, showWarnings = F, recursive = T)
  cmds <- unlist(lapply(ref.strain, cmd.ani, str = strain))
  null <- sapply(cmds, system)
}
run_ANI(OAT.path, blast.path, threads, ref.str, strains, fasta.paths, out.path)
