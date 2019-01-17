# Functions ---------------------------------------------------------------
GetSubjectFeature <-
  function(strain.list,
           feature.path,
           feature.pattern,
           proteome.path,
           proteome.pattern) {
    Read_Query <-
      function(strain,
               featurepath,
               featurepattern,
               proteomepath,
               proteomepattern) {
       # if no feature is present we use the fasta names that prodigal make to get the gene position
        Create_Feature_Table <- function(ProtFasta) {
          fasta <- seqinr::read.fasta(file = ProtFasta, seqtype = "AA")



          headers <- unlist(seqinr::getAnnot(fasta))
          headers <- strsplit(headers, split = "#")

          product_accession <- sapply(headers, function(h){
            trimws(gsub(pattern=">", replacement = "", h))[[1]]
          })
          start <- as.integer(sapply(headers, function(h){
            trimws(h[2])
          }))
          end <- as.integer(sapply(headers, function(h){
            trimws(h[3])
          }))
          genomic_accession <- sapply(product_accession, function(h){
            gsub("_.*", "",h )
          })
          Feature.Subject.NF <-
            data.frame(product_accession,
                       start,
                       end,
                       genomic_accession,
                       stringsAsFactors = F)
          return(Feature.Subject.NF)
        }

        # SUBJECT
        if (file.exists(paste0(featurepath , strain, featurepattern))) {
          Feature.Subject <-
            read.table(
              file = paste0(featurepath , strain, featurepattern),
              header = T,
              sep = "\t",
              comment.char = "",
              quote = "",
              stringsAsFactors = F
            )
          Feature.Subject <-
            Feature.Subject[Feature.Subject[, 1] == "CDS", c(11, 8, 9, 7)]
        } else{

          Feature.Subject <-
            data.frame(Create_Feature_Table(ProtFasta = paste0(
              proteome.path, strain, proteome.pattern
            )),
            stringsAsFactors = F)
        }
        return(Feature.Subject)
      }

    Feature.Subject <- sapply(strain.list, function(s){
      data.frame(Read_Query(strain = s,
                            featurepath = feature.path,
                            featurepattern = feature.pattern,
                            proteomepath = proteome.path,
                            proteomepattern = proteome.pattern),
                 stringsAsFactors = F)
    },simplify = F, USE.NAMES = T)

    return(Feature.Subject)
  }

# Script ------------------------------------------------------------------
args = commandArgs(trailingOnly = T)

if(!require("seqinr")) install.packages("seqinr")

strains <- readLines(args[1])
f.path <- args[2]
if (is.na(args[3])){
  f.pattern <- "_feature_table.txt"
}else{
  f.pattern <- args[3]
}
p.path <- args[4]
if (is.na(args[5])){
  p.pattern <- "_protein.faa"
}else{
  p.pattern <- args[5]
}

featu.list <- GetSubjectFeature(strains, f.path, f.pattern, p.path, p.pattern)

saveRDS(featu.list, "Subject_Feature_list.RDS")
