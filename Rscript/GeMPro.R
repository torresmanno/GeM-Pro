#### ---------------------------------------Functions ----------------------------------------
ReadBlastTables <- function(query.strain,
                            strain.list,
                            inpath){
  strains <- strain.list[!strain.list==query.strain]
  paths <- paste0(inpath, query.strain, "_", strains, ".tsv")
  blast.table.list <- sapply(strains, function(s){
    read.table( paste0(inpath, query.strain, "_", s, ".tsv"),
                header = F,
                fill = T,
                stringsAsFactors = F)
  }, simplify = F, USE.NAMES = T)

  return(blast.table.list)
}

OrthoTables <- function(blast.table.list,
                        cov.cutoff = 50,
                        cutoff.calc = "mean",
                        excluded.genes = NULL,
                        selected.genes = NULL) {



  column.names <- c("Gen",
                    "Subject.ID",
                    "P.Ident",
                    "Length",
                    "Mismatch",
                    "GapOpen",
                    "Qstart",
                    "Qend",
                    "Sstart",
                    "Send",
                    "Evalue",
                    "BitScore",
                    "Coverage" )

  # it is put in a list so it can be formated all together
  # Naming the tables
  blast.table.list <- lapply(blast.table.list, setNames, column.names)

  # Si hay que sacar una lista de genes query
  if (!is.null(excluded.genes)) {
    blast.table.list <- lapply(blast.table.list, function(df){
      df[-grep(paste(excluded.genes, collapse = "|"),
               x = df$Gen),]
    })
  }
  # Si hay que hacer un subset lo hago
  if (!is.null(selected.genes)) {
    blast.table.list <- lapply(blast.table.list, function(df){
      df[grep(paste(selected.genes, collapse = "|"),
              x = df$Gen), ]
    })
  }

  # Orther the tables by Gen and BitScore
  blast.table.list <- lapply(blast.table.list, function(df){
    df[order(df$Gen, -df$BitScore),]
  })

  # Select the uniq values
  tablaII.list <- lapply(blast.table.list, function(df){
    df[!duplicated(df$Gen), ]
  })

  # Orther y Subject and bitscore
  tablaII.list <- lapply(tablaII.list, function(df){
    df[order(df$Subject.ID, -df$BitScore),]
  })

  # Filter result by coverage
  tablaII.list.cov <- lapply(tablaII.list, function(df){
    subset(df, Coverage > cov.cutoff)
  })

  # Select the duplicated result in order to calculate the cutoff value
  tablaII.list.dup <- lapply(tablaII.list.cov, function(df){
    df[duplicated(df$Subject.ID),]
  })

  # Calc the cutoff with the desired method used as parameter
  cutoff <- sapply(tablaII.list.dup, function(df){
    if(nrow(df)!=0){
      if (cutoff.calc == "mean"){
        return(mean(df[, "P.Ident"]))
      }else if(cutoff.calc =="median"){
        return(median(df[, "P.Ident"]))
      }else if(cutoff.calc == "max"){
        return(max(df[, "P.Ident"]))
      }else {
        print("insert mean, median or max in cutoff_calc parameter")
      }
    }else {
      return(NA)
    }
  })

  # If the cutoff is NaN use the mean of the others cutoffs

  cutoff[is.na(cutoff)] <- mean(cutoff[!is.na(cutoff)])

  # Select the uniq result by Subject

  tablaII.list.uniq <- lapply(tablaII.list.cov, function(df){
    df[!duplicated(df$Subject.ID), ]
  })

  # Create the Ortho tables and those who were filter by the cutoff
  TablaIV= mapply(function(t, c){
    subset(t, P.Ident >= c)
  }, tablaII.list.uniq, cutoff, SIMPLIFY = F)

  Out.Cutoff= mapply(function(t, c){
    subset(t, P.Ident < c)
  }, tablaII.list.uniq, cutoff, SIMPLIFY = F)


  # make the Paralogues tables
  Paralogos <- mapply(rbind, tablaII.list.dup, Out.Cutoff, SIMPLIFY = F)

  orthoTablesList <- list(
    Orthologs = TablaIV,
    Paralogs = Paralogos,
    All = tablaII.list.cov,
    Cutoffs = cutoff
  )

  return(orthoTablesList)
}

Get.GeMPro.DB <-
  function(phylo.group.name,
           ortho.tables.list,
           selected.genes = NULL,
           excluded.genes = NULL) {
    # Creation of tables and mean value for BF calculation
    strainList <- names(ortho.tables.list$Orthologs)
    paral.tables <- list()
    ortho.tables <- list()
    paral.tables <- ortho.tables.list$Paralogs
    ortho.tables <- ortho.tables.list$Orthologs

    for (strain in strainList) {
      paral.tables[[strain]]$Strain <- strain
      paral.tables[[strain]]$Group <- phylo.group.name
      ortho.tables[[strain]]$Strain <- strain
      ortho.tables[[strain]]$Group <- phylo.group.name
    }
    paral.tables.unified <- data.table::rbindlist(paral.tables)
    ortho.tables.unified <- data.table::rbindlist(ortho.tables)

    # if we want to remove some genes we can use this option
    if (!is.null(excluded.genes)) {
      paral.tables.unified <-
        paral.tables.unified[-grep(paste(excluded.genes, collapse = "|"),
                                   x = paral.tables.unified$Gen), ]
      ortho.tables.unified <-
        ortho.tables.unified[-grep(paste(excluded.genes, collapse = "|"),
                                   x = ortho.tables.unified$Gen), ]
    }

    # make fun of density
    densidad.ortho = density(
      ortho.tables.unified$P.Ident / 100,
      from = 0.1,
      to = 1,
      bw = .02
    )
    densidad.paral = density(
      paral.tables.unified$P.Ident / 100,
      from = 0.1,
      to = 1,
      bw = .02
    )
    fun.den.ortho = approxfun(densidad.ortho)
    fun.den.paral = approxfun(densidad.paral)

    # Calculate the fraction of orth and paral prior
    total.hits <- sum(nrow(ortho.tables.unified), nrow(paral.tables.unified))

    ortho.fraction <-
      nrow(ortho.tables.unified) / total.hits
    paral.fraction <-
      nrow(paral.tables.unified) / total.hits

    #Create bayes curves
    Bayes.curve.orth <-
      as.data.frame(curve(fun.den.ortho(x) * ortho.fraction, 0.1, 1))
    Bayes.curve.paral <-
      as.data.frame(curve(fun.den.paral(x) * paral.fraction, 0.1, 1))
    # remove NAs
    Bayes.curve.orth$y[is.na(Bayes.curve.orth$y)] <- 0.0
    Bayes.curve.paral$y[is.na(Bayes.curve.paral$y)] <- 0.0

    # Create probability curves

    Prob.orth <-
      as.data.frame(curve((fun.den.ortho(x) * ortho.fraction) / (
        fun.den.ortho(x) * ortho.fraction + fun.den.paral(x) * paral.fraction
      ),
      from = 0.1,
      to = 1
      ))
    Prob.paral <-
      as.data.frame(curve((fun.den.paral(x) * paral.fraction) / (
        fun.den.ortho(x) * ortho.fraction + fun.den.paral(x) * paral.fraction
      ),
      from = 0.1,
      to = 1
      ))

    Prob.orth$y[is.na(Prob.orth$y)] <- 0.0
    Prob.paral$y[is.na(Prob.paral$y)] <- 0.0

    colnames(selected.genes)[c(1,2)] <- c("Gen" , "Description")
    # Select genes to analyse
    all.select <- rbind(ortho.tables.unified, paral.tables.unified)
    Bayes.Table.all <- all.select[, c(1, 2, 3, 14, 15)]
    Bayes.Table.all$Prob_Orth <- (
      fun.den.ortho(Bayes.Table.all$P.Ident / 100) * ortho.fraction
    ) / (
      fun.den.ortho(Bayes.Table.all$P.Ident / 100) * ortho.fraction +
        fun.den.paral(Bayes.Table.all$P.Ident / 100) * paral.fraction)
    Bayes.Table.all$Prob_Paral <- (
      fun.den.paral(Bayes.Table.all$P.Ident / 100) * paral.fraction
    ) / (
      fun.den.ortho(Bayes.Table.all$P.Ident / 100) * ortho.fraction +
        fun.den.paral(Bayes.Table.all$P.Ident / 100) * paral.fraction
    )
    Bayes.Table.all$Ortolog.by.Prob <-
      Bayes.Table.all$Prob_Orth > Bayes.Table.all$Prob_Paral

    Pre.Bayes.Table <-
      merge.data.frame(all.select[, c(1, 2, 3, 13, 14, 15)],
                       selected.genes, by = 1)

    Bayes.Table <-
      data.frame(
        ProteinID_Query = Pre.Bayes.Table[, 1],
        symbol = Pre.Bayes.Table[, 7],
        Pathway = Pre.Bayes.Table[, 8],
        Target_Org = Pre.Bayes.Table[, 5],
        Target_group = Pre.Bayes.Table[, 6],
        ProteinID_Target = Pre.Bayes.Table[, 2],
        ID_perc = Pre.Bayes.Table[, 3],
        Query_Coverage = Pre.Bayes.Table[, 4],
        Prob_Orth = (
          fun.den.ortho(Pre.Bayes.Table$P.Ident / 100) * ortho.fraction
        ) / (
          fun.den.ortho(Pre.Bayes.Table$P.Ident / 100) * ortho.fraction +
            fun.den.paral(Pre.Bayes.Table$P.Ident / 100) * paral.fraction
        ),
        Prob_Paral = (
          fun.den.paral(Pre.Bayes.Table$P.Ident / 100) * paral.fraction
        ) / (
          fun.den.ortho(Pre.Bayes.Table$P.Ident / 100) * ortho.fraction +
            fun.den.paral(Pre.Bayes.Table$P.Ident / 100) * paral.fraction
        ),
        BF_Orth = fun.den.ortho(Pre.Bayes.Table$P.Ident / 100) * ortho.fraction,
        BF_Paral = fun.den.paral(Pre.Bayes.Table$P.Ident / 100) * paral.fraction,
        stringsAsFactors = F
      )
    Bayes.Table$Ortolog.by.Prob <-
      Bayes.Table$Prob_Orth > Bayes.Table$Prob_Paral
    Bayes.Table$Ortolog.by.BF <-
      Bayes.Table$BF_Orth > Bayes.Table$BF_Paral
    Bayes.Table$Score <-
      (Bayes.Table$Prob_Orth / Bayes.Table$Prob_Paral) / 100

    densityApproach <-
      list(
        orthoDensity = fun.den.ortho,
        paralDensity = fun.den.paral,
        orthoBF = Bayes.curve.orth,
        paralBF = Bayes.curve.paral,
        orthoProbability = Prob.orth,
        paralProbability = Prob.paral,
        BayesTable = Bayes.Table,
        BayesTable.all = Bayes.Table.all
      )

    return(densityApproach)
  }


GetQueryFeature <- function(query.feature.file, selected.genes.DF) {
  Feature.Query <-
    read.table(
      file = query.feature.file,
      header = T,
      sep = "\t",
      comment.char = "",
      quote = "",
      stringsAsFactors = F
    )
  # PGP genes
  Feature.Query.PGP <-
    Feature.Query[which(Feature.Query$product_accession %in% selected.genes.DF$Acc_Ver),]
  Feature.Query.PGP <-
    merge.data.frame(
      x = Feature.Query.PGP[-15] ,
      y = selected.genes.DF[c("Acc_Ver", "Gen", "Pathway")],
      by.x = "product_accession",
      by.y = "Acc_Ver",
      all.x = T,
      all.y = F
    )
  names(Feature.Query.PGP)[names(Feature.Query.PGP) == "Gen"] <-
    "symbol"
  return(Feature.Query.PGP)
}

GetSubjectFeatureBayes <- function(feature.subject, bayes.table) {
  bayes.table <-
    bayes.table[order(bayes.table$ProteinID_Target,-bayes.table$ID_perc), ]
  bayesTableUniq <-
    bayes.table[!duplicated(bayes.table$ProteinID_Target), ]

  BayesTableOrth <- bayesTableUniq[bayesTableUniq$Ortolog.by.Prob, ]

  strainList <- names(feature.subject)

  feature.subject.PGP <- sapply(feature.subject, function(df){
    f.pgp <- merge.data.frame(x = df, y = BayesTableOrth,
                              by.x = "product_accession",
                              by.y = "ProteinID_Target")
    f.pgp <- f.pgp[order(f.pgp$genomic_accession, f.pgp$start),]
    return(f.pgp)

  },simplify = F)
  return(feature.subject.PGP)
}

SyntenyTableMaker <- function( feature.Subject,
                               feature.Subject.PGP,
                               feature.Query.PGP){
  Paths <- unique(feature.Query.PGP$Pathway)
  Ngens_S <- sapply(Paths, function(path){
    sum(feature.Subject.PGP$Pathway == path)
  })
  Ngens_Q <- sapply(Paths, function(path){
    sum(feature.Query.PGP$Pathway == path)
  })
  Complete <- Ngens_S /Ngens_Q
  Contig <- sapply(Paths, function(path){
    unique(feature.Subject.PGP$genomic_accession[feature.Subject.PGP$Pathway == path])
  },simplify = F, USE.NAMES = T)

  Order <- sapply(Paths, function(path){
    # NOTE: Here the order is given when the genes are in the same contig
    # If they are in different contig and the contigs aren't ordered the Order = FALSE
    if (Complete[[path]]>= 0.75){
      genes.query  = feature.Query.PGP$symbol[feature.Query.PGP$Pathway == path]
      genes.subject = feature.Subject.PGP$symbol[feature.Subject.PGP$Pathway == path]
      if(identical(genes.query, genes.subject) |
         identical(genes.query, genes.subject[length(genes.subject):1])){
        return(TRUE)
      }else {
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
  })

  Border <- sapply(Paths, function(path){
    if (Complete[[path]]>= 0.75){
      sapply(Contig[[path]], function(contig){
        log.contig <- feature.Subject$genomic_accession == contig
        log.contig.pgp <- feature.Subject.PGP$genomic_accession == contig
        log.path <- feature.Subject.PGP$Pathway == path
        contig_max = max(feature.Subject$end[log.contig])
        contig_pgp_max = max(feature.Subject.PGP$end[log.contig.pgp & log.path])
        contig_min = min(feature.Subject$start[log.contig])
        contig_pgp_min = min(feature.Subject.PGP$start[log.contig.pgp & log.path])

        return((contig_max == contig_pgp_max) | (contig_min == contig_pgp_min))

      },USE.NAMES = F)
    }else{
      return(FALSE)
    }
  })

  length.query <- sapply(Paths, function(path){
    end = max(feature.Query.PGP$end[feature.Query.PGP$Pathway == path])
    start = min(feature.Query.PGP$start[feature.Query.PGP$Pathway == path])
    return (end - start)
  })
  Ratio <- sapply(Paths, function(path){
    if (Complete[[path]]>= 0.75){
      genes.query  = feature.Query.PGP$symbol[feature.Query.PGP$Pathway == path]
      genes.subject = feature.Subject.PGP$symbol[feature.Subject.PGP$Pathway == path]
      missing_genes = setdiff(genes.query, genes.subject)
      if(length(missing_genes)!=0){
        if(length(Contig[[path]]) == 1){
          ends = feature.Subject.PGP$end[feature.Subject.PGP$Pathway ==path]
          starts = feature.Subject.PGP$start[feature.Subject.PGP$Pathway == path]
          length.subject <- sum(ends-starts)
        } else if (length(Contig[[path]]) > 1){
          length.subject <- sum(sapply(Contig[[path]], function(contig){
            ends = feature.Subject.PGP$end[
              feature.Subject.PGP$Pathway == path &
                feature.Subject.PGP$genomic_accession == contig]
            starts = feature.Subject.PGP$start[
              feature.Subject.PGP$Pathway == path &
                feature.Subject.PGP$genomic_accession == contig]
            return(sum(ends-starts))
          }))
        }else {
          length.subject = 0
        }
        delete_gen_lenght <- sum(sapply(missing_genes, function(gene){
          log.gen <- feature.Query.PGP$symbol == gene
          return(feature.Query.PGP$end[log.gen] - feature.Query.PGP$start[log.gen])
        },USE.NAMES = F))
        return((abs(length.query[[path]] - delete_gen_lenght) - length.subject)/
                 (length.query[[path]] - delete_gen_lenght))
      }else{
        if(length(Contig[[path]]) == 1){
          end = max(feature.Subject.PGP$end[feature.Subject.PGP$Pathway == path])
          start = min(feature.Subject.PGP$start[feature.Subject.PGP$Pathway == path])
          length.subject <- end - start
        }else if (length(Contig[[path]]) >1){
          length.subject <- sum(sapply(Contig[[path]], function(contig){
            end = max(feature.Subject.PGP$end[
              feature.Subject.PGP$Pathway == path &
                feature.Subject.PGP$genomic_accession == contig])
            start = min(feature.Subject.PGP$start[
              feature.Subject.PGP$Pathway == path &
                feature.Subject.PGP$genomic_accession == contig])
            return(sum(end-start))
          }))
        }else {
          length.subject <- 0
        }
        return(abs(length.query[[path]] - length.subject)/ length.query[[path]])
      }
    }else{
      return(NA)
    }
  })

  Score.P <- sapply(Paths, function(path){
    if (Complete[[path]] >=0.75){
      sum(feature.Subject.PGP$Score[feature.Subject.PGP$Pathway == path])
    }else{
      return(NA)
    }
  })

  Presence_complete <- Complete == 1.0
  Presence_contig_missing_genes <- sapply(Paths, function(path){
    (Ngens_S[[path]] - 1 + length(Contig[[path]])) / Ngens_Q[[path]] >= 1
  })

  UniContig_Border <- sapply(Paths, function(path){
    length(Contig[[path]]) ==1 & any(Border[[path]])
  })

  Presence_score <- mean(feature.Subject.PGP$Score) * Ngens_S <= Score.P

  Presence_ratio <- Ratio <0.5

  Presence <- sapply(Paths, function(path){

    if(Presence_complete[[path]] & Presence_ratio[[path]] & Presence_score[[path]] & Order[[path]]){
      return(3)
    }else if(Presence_complete[[path]] & Presence_ratio[[path]]){
      return(2)
    }else if(Presence_complete[[path]] | Presence_contig_missing_genes[[path]] & any(Border[[path]]) | UniContig_Border[[path]] | Complete[[path]] >=0.75 & Presence_ratio[[path]]){
      return(1)
    }else {
      return(0)
    }
  })

  Synteny_table <- data.frame(
    Paths = Paths,
    Contigs = sapply(Contig, paste0, collapse = "; "),
    Border = sapply(Border, paste0, collapse ="; "),
    Complete = round(Complete * 100, 0),
    Arrangement = Order,
    Ratio = round(Ratio, digits = 4),
    Score = round(Score.P, 2),
    Presence = Presence,
    stringsAsFactors = F)

  rownames(Synteny_table) <- seq_along(Paths)
  return(Synteny_table)
}

GetSyntenyTable <-
  function( feature.subject,
            feature.subject.pgp,
            feature.query.pgp ) {
    strainList <- names( feature.subject.pgp )
    # Synteny_table_list <- list()
    Synteny_table_list <- mapply( SyntenyTableMaker,
                                  feature.subject,
                                  feature.subject.pgp,
                                  MoreArgs = list( feature.Query.PGP = feature.query.pgp ),
                                  SIMPLIFY = F, USE.NAMES = T )


    Synteny_table_list_present <- lapply(Synteny_table_list, function(df){
      subset(df, Presence != 0)
    })
    return(list(Syn_present = Synteny_table_list_present, Syn_all = Synteny_table_list))
  }

GetPangenome <- function(synteny.table.list, genes.pgp) {
  synteny.list <-  lapply(synteny.table.list, "[[", 1)
  synteny.table <- unlist(unname(synteny.list), recursive = F)
  paths <- unique(genes.pgp$Pathway)

  pangenome <- mapply(function(st, str){
    pg <- data.frame(row.names= paths)
    pg[,str] <- 0
    pg[st$Paths,] <- st$Presence[st$Paths %in% paths]
    return(pg)
  }, synteny.table, names(synteny.table),SIMPLIFY = F)

  pangenomedf <- data.frame(pangenome,check.names = F)
  return(pangenomedf)
}


#### ----------------------------------------RSCRIPT ----------------------------------------
args <- commandArgs(trailingOnly = T)
pdf(NULL)

if (!require("data.table")) install.packages("data.table")
main.strain <- args[1]
phylo.group.df <- read.csv(file = args[2], header = T, col.names = c("Strains", "Groups"), stringsAsFactors = F)

if(main.strain %in% phylo.group.df$Strains){
  with.ref = T
  phylo.group.df <- phylo.group.df[!phylo.group.df$Strains %in% main.strain,]
}else {
  with.ref = F
}
# raw.blast <- readRDS(args[3])
blast.dir <- args[3]
pgp.genes <- read.csv(args[4], header = T, stringsAsFactors = F)
subject.feature <- readRDS(args[5])
query.feature.file <- args[6]
out.folder <- args[7]

phylo.group.list <- sapply(unique(phylo.group.df$Groups), function(g) phylo.group.df$Strains[phylo.group.df$Groups==g],simplify = F, USE.NAMES = T)

strains <- phylo.group.df$Strains

raw.blast <- ReadBlastTables(main.strain, strains, blast.dir)

query.feature.pgp <- GetQueryFeature(query.feature.file, pgp.genes)

raw.blast.list <- list()
subject.feature.list <- list()
for (groups in names(phylo.group.list)){
  raw.blast.list[[groups]] <- raw.blast[phylo.group.list[[groups]][!phylo.group.list[[groups]] == main.strain]]
  subject.feature.list[[groups]] <- subject.feature[phylo.group.list[[groups]][!phylo.group.list[[groups]] == main.strain]]
}
rm(raw.blast, subject.feature)

ortho.tables.list <- lapply(X = raw.blast.list, FUN = OrthoTables)


GeMPro_DB <- mapply(function(ph, ort){
  Get.GeMPro.DB(ph, ortho.tables.list = ort, selected.genes = pgp.genes)
},
names(phylo.group.list),
ortho.tables.list,
SIMPLIFY = F,
USE.NAMES = T
)

Bayes.table <- sapply(GeMPro_DB, "[", 7)

subject.feature.pgp.list <- mapply(GetSubjectFeatureBayes,
                                   subject.feature.list,
                                   Bayes.table,
                                   SIMPLIFY = F)

synteny.table.list <- mapply(GetSyntenyTable,
                             subject.feature.list,
                             subject.feature.pgp.list,
                             MoreArgs = list(
                               feature.query.pgp = query.feature.pgp
                             ),
                             SIMPLIFY = F)

dir.create(out.folder, showWarnings = F, recursive = T)
saveRDS(ortho.tables.list, file = paste0(out.folder, "IWNO.RDS"))
saveRDS(GeMPro_DB, file = paste0(out.folder, "GeMPro_DB.RDS"))
# saveRDS(subject.feature.pgp.list, file = paste0(out.folder, "subject.feature.pgp.RDS"))
saveRDS(synteny.table.list, file = paste0(out.folder, "synteny.list.RDS"))
pangenome <- GetPangenome(synteny.table.list, pgp.genes)
if(with.ref){
  pangenome[,main.strain] <- 3
}

saveRDS(pangenome, file = paste0(out.folder, "pathway.profiles.RDS"))

Colors <- c("#2dc72d","#2cb9d9","#d94c2c","#95362b","#64a871","#6f6398","#cf041c","#ff49ea","#036d5a","#9b014b","#0552d7","#daa520")

itol.data.set <- c(
  "DATASET_EXTERNALSHAPE",
  "SEPARATOR COMMA",
  "DATASET_LABEL, Pathways",
  "COLOR,#00ff00",
  paste0("FIELD_COLORS,", paste0(Colors[1:nrow(pangenome)], collapse = ",")),
  paste0("FIELD_LABELS,",  paste0(rownames(pangenome), collapse = ",")),
  "LEGEND_TITLE,Pathways",
  paste0("LEGEND_SHAPES,", paste0(rep("2", times= nrow(pangenome)), collapse = ",")),
  paste0("LEGEND_COLORS,", paste0(Colors[1:nrow(pangenome)], collapse = ",")),
  paste0("LEGEND_LABELS,",  paste0(rownames(pangenome), collapse = ",")),
  "HEIGHT_FACTOR,5",
  "HORIZONTAL_GRID,0",
  "VERTICAL_GRID,0",
  "SHAPE_SPACING,0",
  "SHAPE_TYPE,2",
  "SHOW_VALUES,0",
  "DATA",
  ""
)
writeLines(itol.data.set, paste0(out.path, "pathway.profiles.txt"))
write.table(t(pangenome), file = paste0(out.path, "pathway.profiles.txt"), quote = F, row.names = T, sep = ",", append = T, col.names = F)
