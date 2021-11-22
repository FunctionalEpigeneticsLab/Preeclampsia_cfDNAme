suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(boot))
suppressMessages(library(limma))
suppressMessages(library(glmnet))
suppressMessages(library(ROCR))
suppressMessages(library(gplots))
suppressMessages(library(ggfortify))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grDevices))

apply_model <- function(input, train_ratio, top_n, alpha){
  train=sample(1:nrow(input), nrow(input)*train_ratio)
  test=(-train)
  train_data = input[train,]
  while(is.na(table(input[-train,'condition1'])[1]) | is.na(table(input[-train,'condition1'])[2]) | table(input[-train,'condition1'])[2] <= 1 | table(input[-train,'condition1'])[1] <= 1){
    train=sample(1:nrow(input), nrow(input)*train_ratio)
    test=(-train)
    train_data = input[train,]
  }
  
  # limma test
  # design <- model.matrix(~ train_data$condition1 + train_data$MeanMethylation+ train_data$GA_NIPT)
  # design <- model.matrix(~ train_data$condition1 +  train_data$FF+ train_data$GA_NIPT+ train_data$exp)
  design <- model.matrix(~ train_data$condition1 +  train_data$FF + train_data$GA_NIPT)
  new_train_data <- t(train_data[,!(colnames(train_data) %in% c("exp","sample","condition1", "GAD.days", "GA_NIPT", "MeanMethylation", "coverageRate", "FF"))])
  fit <- lmFit(new_train_data, design)  # Column 1 contains row-namesb
  fit2 <- eBayes(fit)
  top <- rownames(head(fit2$p.value[order(fit2$p.value[,2], decreasing = F),],top_n))
  
  x = model.matrix(condition1~. , input[,c(top,'condition1')])
  y = droplevels(as.factor(input$condition1))
  y_test = y[test]
  cv.lasso = cv.glmnet(x[train,],y[train], alpha = alpha, family = 'binomial') # number of folds - default is 10.

  # Apply model to testing dataset
  lasso.prob <- predict(cv.lasso,type="response", newx = x[test,], s = 'lambda.min')
  lasso.coef <- predict(cv.lasso,type="coefficients", s = 'lambda.min')[1:ncol(x)+1,]
  nonzero <- list(lasso.coef[lasso.coef!=0])
  pred <- prediction(lasso.prob, y_test)

  # ROC curve
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")@y.values # shows calculated AUC for model
  PE_score<-list(lasso.prob[,1])
  return(c(auc,nonzero,perf, PE_score))
}

