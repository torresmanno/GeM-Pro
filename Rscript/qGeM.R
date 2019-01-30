# FUNCTIONS ---------------------------------------------------------------
# Read Blast ----
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


# Ortho Tables ----
ort.table <- function(df, cov.cutoff=50){

  column.names <- c("ProteinID_Query",
                    "ProteinID_Target",
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
  df <- setNames(df, column.names)

  # Orther the tables by Gen and BitScore
  df <- df[order(df$ProteinID_Query, -df$BitScore), ]

  df <- df[df$Coverage >= cov.cutoff, ]
  # Select the uniq values
  df.II <- df[!duplicated(df$ProteinID_Query), ]
}



# Prediction ----
Predictor <- function(test,classificator) {
  library(MASS)
  test.X <- test[, "P.Ident", drop = F]
  # Preddiction =========
  if (classificator$best_test == "LDA") {
    gempro = predict(classificator$fit, test.X)$posterior[, 2]
    pred = ifelse(gempro > classificator$Threshold, 1, 0)
    test$Class = pred
  } else if (classificator$best_test == "LR") {
    gempro = predict(classificator$fit, test, type = "response")
    pred = ifelse(gempro > classificator$Threshold, 1, 0)
    test$Class = pred
  } else if (classificator$best_test == "KNN") {
    test.X <- test[, "P.Ident", drop = F]
    pred = class::knn(classificator$train.X, test.X, classificator$train.Y, k = classificator$bestK)
    test$Class = pred
  }
  orth.test <- test[test$Class == 1, ]
  north.test <- test[test$Class == 0, ]
  return(orth.test)
}


# Features ----

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
  # Lista de genes PGP
  Feature.Query.PGP <-
    Feature.Query[which(Feature.Query$product_accession %in% selected.genes.DF$Acc_Ver),]
  Feature.Query.PGP <-
    merge.data.frame(
      x = Feature.Query.PGP[-15] ,
      y = selected.genes.DF,
      by.x = "product_accession",
      by.y = "Acc_Ver",
      all.x = T,
      all.y = F
    )
  names(Feature.Query.PGP)[names(Feature.Query.PGP) == "Gen"] <-
    "symbol"
  return(Feature.Query.PGP)
}

GetSubjectFeature.orth <- function(feature.subject, orth.table, pgp.df) {

  names(pgp.df)[names(pgp.df) == "Gen"] <- "symbol"
  feature.subject.PGP <- list()
  for (strain in names(orth.table)){
    df <- merge(x = feature.subject[[strain]], y = orth.table[[strain]],
                                           by.x = "product_accession", by.y = "ProteinID_Target")
    df2 <- merge(x = df, y = pgp.df, by.x = "ProteinID_Query", by.y = "Acc_Ver")
    feature.subject.PGP[[strain]] <- df2[order(df2$genomic_accession, df2$start),]
  }
  return(feature.subject.PGP)
}

# Sinteny -----

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
            return(sum(end-start))
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

  # Score.P <- sapply(Paths, function(path){
  #   if (Complete[[path]] >=0.75){
  #     sum(feature.Subject.PGP$Score[feature.Subject.PGP$Pathway == path])
  #   }else{
  #     return(NA)
  #   }
  # })

  Presence_complete <- Complete == 1.0
  Presence_contig_missing_genes <- sapply(Paths, function(path){
    (Ngens_S[[path]] - 1 + length(Contig[[path]])) / Ngens_Q[[path]] >= 1
  })

  UniContig_Border <- sapply(Paths, function(path){
    length(Contig[[path]]) ==1 & any(Border[[path]])
  })

  # Presence_score <- mean(feature.Subject.PGP$Score) * Ngens_S <= Score.P

  Presence_ratio <- Ratio <0.5

  Presence <- sapply(Paths, function(path){

    if(Presence_complete[[path]] & Presence_ratio[[path]]  & Order[[path]]){
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
    # Score = round(Score.P, 2),
    Presence = Presence,
    stringsAsFactors = F)

  rownames(Synteny_table) <- seq_along(Paths)
  return(Synteny_table)
}

GetSyntenyTable <-
  function( feature.subject,
            feature.subject.pgp,
            feature.query.pgp ) {

    feature.subject.pgp <- lapply(feature.subject.pgp, function(df){
      names(df)[names(df) == "BitScore"] <- "Score"
      return(df)
    })
    Synteny_table_list <- list()
    for(strain in names(feature.subject.pgp)){
      Synteny_table_list[[strain]] <-
        SyntenyTableMaker(
          feature.Subject = feature.subject[[strain]],
          feature.Subject.PGP = feature.subject.pgp[[strain]],
          feature.Query.PGP = feature.query.pgp
        )
    }

    Synteny_table_list_present <- lapply(Synteny_table_list, function(df){
      subset(df, Presence != 0)
    })
    return(list(Syn_present = Synteny_table_list_present, Syn_all = Synteny_table_list))
  }

# Pangenome ---------------------------------------------------------------

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

#### RSCRIPT ----
# IMPORTANT: needed the library mass for predict.lda
args <- commandArgs(TRUE)

main.strain <- args[1]
phylo.group.df <- read.csv(file = args[2], header = T, col.names = c("Strains", "Groups"), stringsAsFactors = F)
blast.path <- args[3]
pgp.genes <- read.csv(args[4], header = T, stringsAsFactors = F)
best_class <- args[5]
# group_name <- gsub("_Best_class.RDS","", basename(args[5]))
subject.feature <- readRDS(args[6])
query.feature.file <- args[7]
out.path <- args[8]

# check if reference strain is in the strains to compare.
if(main.strain %in% phylo.group.df$Strains){
  with.ref = T
  phylo.group.df <- phylo.group.df[!phylo.group.df$Strains %in% main.strain,]
}else {
  with.ref = F
}

groups <- unique(phylo.group.df$Groups)
strains <-
  sapply(groups, function(g)
    phylo.group.df$Strains[phylo.group.df$Groups == g], simplify = F, USE.NAMES = T)

best_class_list <-
  sapply(groups, function(g)
    readRDS(paste0(best_class, g, "_Best_class.RDS")), simplify = F, USE.NAMES = T)

raw.blast <-
  sapply(groups, function(g)
    ReadBlastTables(main.strain, strains[[g]], blast.path), simplify = F, USE.NAMES = T)

test.list <-
  sapply(groups, function(g)
    lapply(raw.blast[[g]], ort.table), simplify = F, USE.NAMES = T)

Prediction_list <-
  sapply(groups, function(g)
    lapply(test.list[[g]], Predictor, classificator = best_class_list[[g]]), simplify = F, USE.NAMES = T)

query.feature.pgp <- GetQueryFeature(query.feature.file, pgp.genes)

subject.feature.pgp.list <-
  sapply(groups, function(g)
    GetSubjectFeature.orth(subject.feature, Prediction_list[[g]], pgp.genes), simplify = F, USE.NAMES = T)

subject.feature.list <-
  sapply(groups, function(g)
    subject.feature[strains[[g]]])
synteny.table.list <-
  sapply(groups, function(g)
    GetSyntenyTable(
      subject.feature.list[[g]],
      subject.feature.pgp.list[[g]],
      feature.query.pgp = query.feature.pgp
    ), simplify = F, USE.NAMES = T)

dir.create(out.path, showWarnings = F,recursive = T)

saveRDS(Prediction_list, paste0(out.path, "/orthologous.RDS"))
saveRDS(synteny.table.list, paste0(out.path, "/synteny.RDS"))
pangenome <- GetPangenome(synteny.table.list, pgp.genes)

# If the ref strain is in the list of genomes to compare it add the main strain with all the pathways.
if(with.ref){
  pangenome[,main.strain] <- 3
}

# Create ITOL file
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

saveRDS(pangenome, file =paste0(out.path, "/pathway.profiles.RDS"))

writeLines(itol.data.set, paste0(out.path, "pathway.profiles.txt"))
write.table(t(pangenome), file = paste0(out.path, "pathway.profiles.txt"), quote = F, row.names = T, sep = ",", append = T, col.names = F)
