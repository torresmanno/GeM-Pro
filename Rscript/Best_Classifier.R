# Functions
Best_Classificator <- function(Bayes.table.all){

  # Classiffication
  BestThreshold <- function(trait, classif, thresholds, cl){
    set.seed(1234)
    folds = sample(rep(1:10, length = nrow(training)))
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    if (classif =="LDA"){
      # error <- parallel::parSapply(cl, 1:10, FUN = function(i){
      errors <- foreach::foreach(i = 1:10, .combine = "rbind") %dopar%{
        fit = MASS::lda(Class ~ P.Ident, data = trait[folds !=i, ])
        probs = predict(fit, trait[folds == i, ])$posterior[, 2]
        sapply(1:length(thresholds), function(j){
          pred = ifelse(probs > thresholds[j], 1, 0)
          mean(pred != trait[folds == i, "Class", drop = F])
        })
      }

    } else if (classif== "LR") {
      # error <- parallel::parSapply(cl, 1:10, FUN = function(i){
      errors <- foreach::foreach(i = 1:10, .combine = "rbind") %dopar% {
        fit = glm(Class ~ P.Ident , data = trait[folds != i, ], family = "binomial")
        probs = predict(fit, trait[folds == i, ],type = "response")
        sapply(1:length(thresholds), function(j){
          pred = ifelse(probs > thresholds[j], 1, 0)
          mean(pred != trait[folds == i, "Class", drop = F])
        })
      }
    }
    # parallel::stopCluster(cl)

    CVerrors = colMeans(errors)
    names(CVerrors) = thresholds

    selectedThres = thresholds[which.min(CVerrors)]

    return(data.frame(Threshold = selectedThres, error = min(CVerrors)))
  }

  getBestThreshold <- function(training, classifficator){
    no_cores <- parallel::detectCores() - 1

    # Initiate cluster
    cl <- parallel::makeCluster(no_cores)


    thresholds = seq(0.1, 1, 0.1)
    bestThres <- BestThreshold(training, classifficator, thresholds, cl)
    selectedThres <- bestThres$Threshold

    threshold01 <- seq(selectedThres - 0.09, selectedThres + 0.09, 0.01)
    bestThres <- BestThreshold(training, classifficator, threshold01, cl)
    selectedThres <- bestThres$Threshold

    threshold001 = seq(selectedThres - 0.009, selectedThres + 0.009, 0.001)

    bestThres <- BestThreshold(training, classifficator, threshold001, cl)
    selectedThres <- bestThres$Threshold

    threshold0001 = seq(selectedThres - 0.0009, selectedThres + 0.0009, 0.0001)
    bestThresFinal <- BestThreshold(training, classifficator, threshold0001, cl)

    parallel::stopCluster(cl)
    return(bestThresFinal)
  }

  Bayes.table.all$Class <- ifelse(Bayes.table.all$Ortolog.by.Prob, 1, 0)
  Bayes.sub <- Bayes.table.all[,c(1, 2, 3, 9)]
  Bayes.o <- Bayes.sub[order(Bayes.sub$Subject.ID, -Bayes.sub$P.Ident),]
  Bayes.o$Class[duplicated(Bayes.o$Subject.ID)] <- 0
  data <- Bayes.o[,c(3, 4)]

  set.seed(1)
  Train <- sample(c(T, F), nrow(data), T, c(0.9, 0.1))
  training <- data[Train, , drop=T]

  Test <- !Train
  test <- data[Test, , drop=F]

  # LDA ===========
  bestThresholdLDA <- getBestThreshold(training, "LDA")
  selectedThresLDA <- bestThresholdLDA$Threshold

  fitLDA = MASS::lda(Class ~ P.Ident, data = training)
  probsLDA = predict(fitLDA, training)$posterior[, 2]
  predLDA = ifelse(probsLDA > selectedThresLDA, 1, 0)
  LDAerror = mean(predLDA != training[, "Class", drop = F])

  # Regression ============
  bestThresholdLR <- getBestThreshold(training, "LR")
  selectedThresLR <- bestThresholdLR$Threshold

  fitLR = glm(Class ~ P.Ident, data = training, family = "binomial")
  probsLR = predict(fitLR, training, type = "response")
  predLR = ifelse(probsLR > selectedThresLR, 1, 0)
  LRerror = mean(predLR != training[, "Class", drop = F])

  # Best Class ----

  if (LDAerror <= LRerror) {
    test.X <- test[, "P.Ident", drop = F]
    probs = predict(fitLDA, test.X)$posterior[, 2]
    pred = ifelse(probs > selectedThresLDA, 1, 0)
    LDA_Test_Error = mean(pred != test[, "Class", drop = F])
    Best_class <- list(
      LDAerror = LDA_Test_Error,
      Threshold = selectedThresLDA,
      fit = fitLDA,
      best_test = "LDA",
      CM = caret::confusionMatrix(as.factor(pred), as.factor(test$Class), positive = "1")
    )
  } else if (LRerror <= LDAerror){
    probs = predict(fitLR, test, type = "response")
    pred = ifelse(probs > selectedThresLR, 1, 0)
    LR_Test_Error = mean(pred != test[, "Class", drop = F])
    Best_class <- list(
      LRerror = LR_Test_Error,
      Threshold = selectedThresLR,
      fit = fitLR,
      best_test = "LR",
      CM = caret::confusionMatrix(as.factor(pred), as.factor(test$Class), positive = "1")
    )
  }
  return(Best_class)
}

# Script

args <- commandArgs(trailingOnly = T)
dl <- readRDS(args[1])
out.path <- args[2]


if (!require("caret")) install.packages("caret")
if (!require("MASS")) install.packages("MASS")
if (!require("foreach")) install.packages("foreach")
if (!require("doParallel")) install.packages("doParallel")
if (!require("parallel")) install.packages("parallel")
if (!require("methods")) install.packages("methods")
if (!require("e1071")) install.packages("e1071")

bt.all <- lapply(dl, "[[", 8)
library(methods)
library(e1071)

dir.create(out.path, showWarnings = F, recursive = T)
mapply(function(bt, name){
  saveRDS(Best_Classificator(bt), paste0(out.path, "/", name, "_Best_class.RDS"))
}, bt.all, names(bt.all))
# Best_Class <- lapply(bt.all, Best_Classificator)

# saveRDS(Best_Class, "Best_Class.RDS")
