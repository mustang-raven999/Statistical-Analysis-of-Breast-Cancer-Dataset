# install.packages("mlbench")
library(mlbench)
data("BreastCancer")
dim(BreastCancer)
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("GGally")

library(ggplot2)



# 699 rows and 11 columns
# n = 699
# p = 10
# 1 cell id column == samples

head(BreastCancer)

# exporting data to excel
# write.csv(BreastCancer , file = "BreastCancer.csv")
# df.BreastCancer <- na.omit(BreastCancer)

# renamed
nbc <- BreastCancer
nbc


# removing na values

nbc <- na.omit(nbc)

# Removing ID column as it is casuing issue in Regression
#df = subset(mydata, select = -c(x,z) )

nbc <- subset(nbc, select = -c(Id))
head(nbc)


# changing ordinal values to numeric 
str(BreastCancer)

nbc$Cl.thickness <- as.numeric(nbc$Cl.thickness)
nbc$Cell.size <- as.numeric(nbc$Cell.size)
nbc$Cell.shape <- as.numeric(nbc$Cell.shape)
nbc$Marg.adhesion <- as.numeric(nbc$Marg.adhesion)
nbc$Epith.c.size <- as.numeric(nbc$Epith.c.size)
nbc$Bare.nuclei <- as.numeric(nbc$Bare.nuclei)
nbc$Bl.cromatin <- as.numeric(nbc$Bl.cromatin)
nbc$Normal.nucleoli <- as.numeric(nbc$Normal.nucleoli)
nbc$Mitoses <- as.numeric(nbc$Mitoses)

str(nbc)
str(BreastCancer)






#summaries

summary(nbc)


# Understanding relationship between different predictor variables and response
# variable

# cell thickness v Class

# cell_th <- ggplot(nbc , aes(x=Cl.thickness, fill=Class, color = Class)) +
#   
#   geom_histogram(aes(y=density), position="identity", alpha=0.3, bins=10 )+
#   geom_density(alpha=0.8)+scale_color_manual(values=c("darkorange", "darkorchid"))+
# labs(title="Cell thickness v Class",x="Cl.thickness", y = "Density")
# 
# cell_th
library(ggplot2)

library(GGally)
install.packages("GGally")
data(nbc)
ggpairs(nbc, columns = 1:9, ggplot2::aes(colour=Class))

# #Scatterplots between two predictor variables 
# ggplot(nbc, aes(x = Cl.thickness, y = Bl.cromatin, colour = Class)) +
#   geom_point() +
#   scale_colour_manual(values = c("darkorange","purple")) 
# Scatterplots are not a good measure to find any correlation between the variables

#####################################################
####################################################

#  Logistic regression

logistic <- glm(Class~ ., data = nbc, family = "binomial")
summary(logistic)

# from the summary of the logistic regression we get that Cl.thickness, Marg.adhesion
# Bare.nuclei, Bl.cromatin are the 4 predictor variables that have p value<0.05
# and z value>abs[2] ie; more than 2 standard deviation away.
# This means that they have a statistically significant role in predicting the class
# ie; wheather a cell is benign or malignent

###########################################################

### Best Subset selection

# bestglm package
install.packages("bestglm")
library(bestglm)

bst_fit_aic <- bestglm(nbc , family = binomial , IC = "AIC")
bst_fit_bic <- bestglm(nbc, family = binomial , IC= "BIC")


bst_fit_aic
bst_fit_bic

bst_fit_aic$Subsets
bst_aic_mod = bst_fit_aic$ModelReport$Bestk
bst_aic_mod

# model with loweat AIC score is selected

# bst_fit_aic$Subsets shows us that M7 is the best model that can be used to predict
# the response variable. This means that a cluster of 7 variable model is used to
# predict the class of the cell ie: benign or Malignent
#
# We reached this conclusion after seeing True values for all variables except
# Cell.size and Epith.c.size. and also m7 having the lowest aic score.

bst_fit_bic$Subsets

bst_fit_modbic = bst_fit_bic$ModelReport$Bestk
bst_fit_modbic

# model with loweat BIC score is selected

# bst_fit_$Subsets shows us that M5 is the best model that can be used to predict
# the response variable. This means that a cluster of 5 variable model is used to
# predict the class of the cell ie: benign or Malignent
#
# We reached this conclusion after seeing True values for all variables except
# Cell.size , cell.shape ,Mitosis and Epith.c.size. and also m7 having
# the lowest aic score.

## Create multi-panel plotting device
p = ncol(nbc[1:ncol(nbc)-1])

par(mfrow=c(1,2))
## Produce plots, highlighting optimal value of k
plot(0:p, bst_fit_aic$Subsets$AIC, xlab="Number of predictors", ylab="AIC", type="b")
points(bst_aic_mod, bst_fit_aic$Subsets$AIC[bst_aic_mod +1], col="red", pch=16)
plot(0:p, bst_fit_bic$Subsets$BIC, xlab="Number of predictors", ylab="BIC", type="b")
points(bst_fit_modbic, bst_fit_bic$Subsets$BIC[bst_fit_modbic +1], col="red", pch=16)

##### Now creating a logistic regression for a model with our 5 predictors

pstar = 5
bst_fit_bic$Subsets[pstar+1 , ]

indices = as.logical(bst_fit_bic$Subsets[pstar+1, 2:(p+1)])
indices

new_nbc = data.frame(nbc[,indices])
head(new_nbc)

logfit = glm(Class~., data = new_nbc, family = "binomial")
summary(logfit)



#######################################################################
#######################################################################
### Ridge and Lasso


library(glmnet)

# Choose grid of values for the tuning parameter
grid = 10^seq(-4, 1, length.out = 100)
# Fit a model with LASSO penalty for each value of the tuning parameter
lasso_fit = glmnet(nbc[,1:9], nbc$Class, family="binomial", alpha=1, standardize=FALSE, lambda=grid)
# Examine the effect of the tuning parameter on the parameter estimates
plot(lasso_fit, xvar = "lambda", col = rainbow(p), label = TRUE)

# chossing the optimal lambda value

lasso_cv_fit = cv.glmnet(as.matrix(nbc[,1:9]), nbc$Class, family="binomial", alpha=1, standardize=FALSE, lambda=grid , type.measure= "class")

plot(lasso_cv_fit)

min_lambda = lasso_cv_fit$lambda.min
min_lambda


grid = 10^seq(-4,1, length.out=100)
lasso_cv = cv.glmnet(as.matrix(nbc[,1:9]), nbc$Class, alpha = 1, family = "binomial", lambda = grid, type.measure = "class")


a <- lasso_cv$lambda.min
a
plot(lasso_cv)


coef(lasso_cv_fit, s = min_lambda)
coef(lasso_cv, s = a)

###########################################################
##############################################################

# LDA
# install.packages("topicmodels")
# library(topicmodels)

install.packages("MASS")
library(MASS)

lda_fit = lda(Class ~ ., data = nbc)
lda_fit
plot(lda_fit)

#################################################################


#Cross validation

# n = nrow(nbc) - 1
#
# # 10-fold cross validation...
# nfolds = 10
#
# # Sample fold-assignment index
# fold_index = sample(nfolds, n, replace=TRUE)
#
# head(fold_index)

#
# install.packages("caret")
# library(caret)
# set.seed(2000)
# new_nbc$Class = as.factor(new_nbc$Class)
# str(nbc$Class)
#
# str(new_nbc)
#
# index = createDataPartition(new_nbc$Class, p=.8, list = FALSE, times = 1)
#
# nbc_train = new_nbc[index,]
# nbc_test = new_nbc[-index,]
#
#
# ### Test error for training and testing
#
# # Logistic BIC subset Test error
#
# logfit_train = glm(Class~ ., data = nbc_train, family = "binomial")
# summary(logfit_train)
#
# phat_test = predict(logfit_train , nbc_test , type = "response")
# yhat_test = ifelse(phat_test > 0.5, 'malignant', 'benign')
#
# ##Compute test Error
#
# error_logistic = 1 - mean(nbc_test$Class == yhat_test)
# error_logistic
#
#
# # Lasso Penalty test error
#
# head(new_nbc)
# lasso_train = glmnet()
#
# lasso_fit = glmnet(nbc[,1:9], nbc$Class, family="binomial", alpha=1, standardize=FALSE, lambda=grid)
# d

###########################################
###########################################

# Cross validation using test error and 10 fold method


# number of rows..
n = nrow(df) - 1

# 10-fold cross validation...
nfolds = 10

# Sample fold-assignment index
fold_index = sample(nfolds, n, replace=TRUE)

# Print first few fold-assignments
head(fold_index)


class_cv = function(nbc, new_nbc, fold_ind, model_type, optim_lambda){
  nfolds = max(fold_ind) # number of folds...
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.") # validating partitions...
  cv_errors=numeric(nfolds) # empty list of size equal to number of folds...
  for(fold in 1:nfolds) # fitting different models based on model_type for each fold and making prediction then finding yhat(predicted class)...
  {
    if(model_type == 'IC'){ # Logistic regression using features selected from subset selection AIC and BIC...
      tmp_fit = glm(Class ~ ., data = new_nbc[fold_ind != fold,], family = "binomial")
      pred = predict(tmp_fit, new_nbc[fold_ind == fold, 1:ncol(new_nbc)-1])
      yobs = new_nbc$Class[fold_ind == fold]
      yhat = ifelse(pred > 0.5, 'malignant', 'benign')
    }
    else if(model_type == 'LASSO'){ # LASSO Regularized Logistic Regression...
      tmp_fit = glmnet(nbc[fold_ind != fold, 1:9], nbc[fold_ind != fold, 10], alpha = 1, family = "binomial", lambda = optim_lambda)
      pred = predict(tmp_fit, as.matrix(nbc[fold_ind == fold, 1:9]))
      yobs = nbc$Class[fold_ind == fold]
      yhat = ifelse(pred > 0.5, 'malignant', 'benign')
    }
    else if(model_type == 'LDA'){ # Linear discriminant analysis...
      tmp_fit = lda(Class ~ ., data = nbc[fold_ind != fold,])
      pred = predict(tmp_fit, nbc[fold_ind == fold, 1:9])
      yobs = nbc$Class[fold_ind==fold]
      yhat = ifelse(pred$x > 0.5, 'malignant', 'benign')
    }
    cv_errors[fold] = 1 - mean(yobs == yhat) # Calculate the training error
  }
  # Calculating mean error of all the folds...
  fold_sizes = numeric(nfolds) # placeholder for length of each fold...
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold)) # calculating length of each fold...
  test_error = weighted.mean(cv_errors, w = fold_sizes) # calculating weighted mean of the cross validation errors...
  return(test_error)
}

# library(glmnet)
# library(MASS)
# model_type = IC for AIC, BIC logistic regression...
# model_type = LASSO for logistic regression with lasso regularization...
# model_type = LDA for Linear Discriminant Analysis...
model_type = c('IC', 'LASSO', 'LDA')
error_list = c()
for(i in model_type){
  error = class_cv(nbc, new_nbc, fold_index, i, lasso_cv$lambda.min) # performing cross validation on each type of models
  error_list = append(error_list, error)
}


# Creating dataframe of model used and thier error rate values...
model_name = c('Logistic Regression with BIC', 'Logistic Regression with LASSO regularization', 'Linear Discriminant Analysis')
nbc_error = data.frame(model_name, error_list)
colnames(nbc_error) = c("Model", "Error Rate")
nbc_error

