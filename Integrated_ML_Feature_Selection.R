library(glmnet)
library(randomForest)
library(e1071)
library(caret)
library(sigFeature)
library(limma)
library(ggpubr)
library(VennDiagram)

source('msvmRFE.R') 

input_file <- "merge.txt"
disease_file <- "disease.txt"
sample_file <- "sample.txt"
output_dir <- "./results/"
if(!dir.exists(output_dir)) dir.create(output_dir)

rt <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
data_mat <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = list(rownames(exp), colnames(exp)))
data_mat <- avereps(data_mat)
data_mat <- t(data_mat)

disease_genes <- read.table(disease_file, header = FALSE, sep = "\t")[, 1]
data_mat <- data_mat[, disease_genes]
sample_info <- read.table(sample_file, sep = "\t", header = FALSE, row.names = 1)
data_mat <- data_mat[rownames(sample_info), ]

X_matrix <- as.matrix(data_mat)
num_con <- as.vector(as.matrix(table(sample_info[, 1]))[1, 1])
group_labels <- c(rep("0", num_con), rep("1", nrow(data_mat) - num_con))
Y_vector <- as.matrix(group_labels)
Y_factor <- as.factor(group_labels)

# lasso
set.seed(123)
cv_fit <- cv.glmnet(X_matrix, Y_vector, family = "binomial", nfolds = 10, alpha = 1)
lasso_fit <- glmnet(X_matrix, Y_vector, family = "binomial")
coef_mat <- coef(lasso_fit, s = cv_fit$lambda.min)
lasso_genes <- row.names(coef_mat)[which(coef_mat != 0)]
lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"]

write.table(lasso_genes, file = paste0(output_dir, "LASSO_selected_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(paste0(output_dir, "Fig6A_LASSO.pdf"), width = 7, height = 5)
layout(matrix(c(1,1,2,2), 2, 2, byrow = FALSE))
plot(lasso_fit, xvar = "lambda")
plot(cv_fit)
dev.off()

# randomforest
set.seed(123)
rf_init <- randomForest(x = X_matrix, y = Y_factor, ntree = 500)
opt_trees <- which.min(rf_init$err.rate[, 1])
set.seed(123)
rf_final <- randomForest(x = X_matrix, y = Y_factor, ntree = opt_trees)
importance_df <- as.data.frame(importance(rf_final))
rf_genes <- rownames(importance_df)[importance_df[,1] > 2]

write.table(rf_genes, file = paste0(output_dir, "RF_selected_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(paste0(output_dir, "Fig6B_RF.pdf"), width = 6, height = 6)
varImpPlot(rf_final, main = "Random Forest Variable Importance")
dev.off()

# sum
set.seed(123)
svm_input <- as.data.frame(cbind(Type = group_labels, X_matrix))
svm_input$Type <- as.factor(svm_input$Type)
nfold <- 10
folds <- rep(1:nfold, len = nrow(svm_input))[sample(nrow(svm_input))]
folds_list <- lapply(1:nfold, function(x) which(folds == x))
svm_results <- lapply(folds_list, svmRFE.wrap, svm_input, k = 10, halve.above = 100)
top_features <- WriteFeatures(svm_results, svm_input, save = FALSE)
feat_sweep <- lapply(1:min(nrow(top_features), 20), FeatSweep.wrap, svm_results, svm_input)
errors <- sapply(feat_sweep, function(x) ifelse(is.null(x), NA, x$error))
opt_num <- which.min(errors)
svm_genes <- top_features$FeatureName[1:opt_num]

write.table(svm_genes, file = paste0(output_dir, "SVM_selected_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(paste0(output_dir, "Fig6C_SVM.pdf"), width = 5, height = 5)
PlotErrors(errors, no.info = min(prop.table(table(svm_input$Type))))
dev.off()

consensus_genes <- intersect(intersect(lasso_genes, rf_genes), svm_genes)
write.table(consensus_genes, file = paste0(output_dir, "Consensus_Genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

venn_obj <- venn.diagram(
  x = list(LASSO = lasso_genes, RF = rf_genes, SVM = svm_genes),
  filename = NULL,
  fill = c("#E64B35FF", "#4DBBD5FF", "#00A087FF"),
  alpha = 0.6, cex = 1.5, fontface = "bold", cat.cex = 1.2
)
pdf(paste0(output_dir, "Fig6D_Venn.pdf"), width = 6, height = 6)
grid.draw(venn_obj)
dev.off()

writeLines(capture.output(sessionInfo()), paste0(output_dir, "SessionInfo.txt"))
