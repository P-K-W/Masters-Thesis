## Load packages
library(readxl)
library(caret)
library(e1071)
library(corrplot)

## Read data
data <- read_excel("indecies_data.xlsx")
data$Group <- factor(data$Group, levels = c("HC","OCD"))

## Compute correlation matrix & heatmap
correlation_matrix <- cor(data[, 4:10])
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200))
png("correlation_heatmap_corrplot.png")
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200))
dev.off()

## Define predictors
all_preds    <- setdiff(names(data), c("Group","Dx"))
numeric_preds<- all_preds[sapply(data[all_preds], is.numeric)]
other_preds  <- setdiff(all_preds, numeric_preds)

## Auto‐select GOF & anticorrelation measures
gof_measures  <- grep("GOF", numeric_preds, value=TRUE)
anti_measures <- grep("isoQ|FractalDim", numeric_preds, value=TRUE)
other_preds   <- setdiff(all_preds, c(gof_measures, anti_measures))

## Build feature‐set list (3 × 2 = 6 combinations)
feature_sets <- list()
for(g in gof_measures){
  for(a in anti_measures){
    nm <- paste(g, "vs", a, sep = "_")
    feature_sets[[nm]] <- c(g, a, other_preds)
  }
}

## 5-fold cv 10 repeat nested
n_repeats <- 10      # repeat the entire 5-fold procedure 10 times
n_folds   <- 5       # 5 outer folds, and 5 inner folds for tuning

## Define inner-loop tuning grids as Jonsson et al. (2018)
grid_rbf <- expand.grid(
  C     = c(0.001, 0.01, 0.1, 1, 10),
  sigma = c(0.001, 0.01, 0.1)
)

grid_poly <- expand.grid(
  degree = c(2, 3),
  scale  = c(0.01, 0.1, 1),
  C      = c(0.001, 0.01, 0.1, 1, 10)
)

## Prepare objects to hold results for RBF and Poly
results_rbf <- data.frame(
  GOF      = character(),
  Anti     = character(),
  MeanAUC  = numeric(),
  SD_AUC   = numeric(),
  BestParams = I(list()),
  stringsAsFactors = FALSE
)

results_poly <- data.frame(
  GOF      = character(),
  Anti     = character(),
  MeanAUC  = numeric(),
  SD_AUC   = numeric(),
  BestParams = I(list()),
  stringsAsFactors = FALSE
)

## Begin nested‐CV (outer loop = 5 folds; repeat the 5-fold outer loop 10 times)
set.seed(123)  

for (nm in names(feature_sets)) {
  preds <- feature_sets[[nm]]
  parts <- strsplit(nm, "_vs_")[[1]]
  gof   <- parts[1]
  anti  <- parts[2]
  
  X_full <- data[, preds]
  y_full <- data$Group
  
  auc_values_rbf  <- c()
  best_tunes_rbf  <- list()
  
  auc_values_poly <- c()
  best_tunes_poly <- list()
  
  for (rep in 1:n_repeats) {
    seed_rep <- 123 + rep
    
    set.seed(seed_rep)
    
    folds <- createFolds(y_full, k = n_folds, list = TRUE, returnTrain = FALSE)

    for (fold_idx in seq_along(folds)) {
      test_idx  <- folds[[fold_idx]]
      train_idx <- setdiff(seq_len(nrow(data)), test_idx)
      
      X_train <- X_full[train_idx, , drop = FALSE]
      y_train <- y_full[train_idx]
      X_test  <- X_full[test_idx,  , drop = FALSE]
      y_test  <- y_full[test_idx]
      
      #RBF tuning
      ctrl_inner_rbf <- trainControl(
        method          = "cv",
        number          = n_folds,
        classProbs      = TRUE,
        summaryFunction = twoClassSummary,
        savePredictions = "none"  
      )
      fit_rbf_inner <- train(
        x          = X_train,
        y          = y_train,
        method     = "svmRadial",
        metric     = "ROC",
        trControl  = ctrl_inner_rbf,
        tuneGrid   = grid_rbf,
        preProcess = c("center", "scale")
      )
      # Record the best‐tuned parameters for this outer fold
      best_tunes_rbf <- append(best_tunes_rbf, list(fit_rbf_inner$bestTune))
      
      # Evaluate RBF on the OUTER test set
      probs_rbf <- predict(fit_rbf_inner, newdata = X_test, type = "prob")[, "OCD"]
      roc_rbf   <- roc(response = y_test, predictor = probs_rbf, levels = c("HC", "OCD"), direction = "<")
      auc_values_rbf <- c(auc_values_rbf, as.numeric(auc(roc_rbf)))
      
      # POLY tuning
      ctrl_inner_poly <- trainControl(
        method          = "cv",
        number          = n_folds,
        classProbs      = TRUE,
        summaryFunction = twoClassSummary,
        savePredictions = "none"
      )
      fit_poly_inner <- train(
        x          = X_train,
        y          = y_train,
        method     = "svmPoly",
        metric     = "ROC",
        trControl  = ctrl_inner_poly,
        tuneGrid   = grid_poly,
        preProcess = c("center", "scale")
      )
      # Record best poly parameters for this outer fold
      best_tunes_poly <- append(best_tunes_poly, list(fit_poly_inner$bestTune))
      
      # Evaluate Poly on the OUTER test set
      probs_poly <- predict(fit_poly_inner, newdata = X_test, type = "prob")[, "OCD"]
      roc_poly   <- roc(response = y_test, predictor = probs_poly, levels = c("HC", "OCD"), direction = "<")
      auc_values_poly <- c(auc_values_poly, as.numeric(auc(roc_poly)))
    }
  }
  
  # Compute mean & SD of the 50 AUCs
  mean_auc_rbf <- mean(auc_values_rbf)
  sd_auc_rbf   <- sd(auc_values_rbf)
  
  mean_auc_poly <- mean(auc_values_poly)
  sd_auc_poly   <- sd(auc_values_poly)
  
  # Append to results data.frames
  results_rbf <- rbind(
    results_rbf,
    data.frame(
      GOF        = gof,
      Anti       = anti,
      MeanAUC    = mean_auc_rbf,
      SD_AUC     = sd_auc_rbf,
      BestParams = I(list(best_tunes_rbf)),
      stringsAsFactors = FALSE
    )
  )
  
  results_poly <- rbind(
    results_poly,
    data.frame(
      GOF        = gof,
      Anti       = anti,
      MeanAUC    = mean_auc_poly,
      SD_AUC     = sd_auc_poly,
      BestParams = I(list(best_tunes_poly)),
      stringsAsFactors = FALSE
    )
  )
}

# Print final results
cat("=== RBF SVM Results (Nested 5‐fold × 10 Repeats) ===\n")
print(results_rbf)
cat("\n=== Polynomial SVM Results (Nested 5‐fold × 10 Repeats) ===\n")
print(results_poly)








