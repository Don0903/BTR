##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##----------------------------- BACHELOR THESIS---------------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ PACKAGES AND LIBRARIES  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("ARTofR")
install.packages("rlang")
install.packages("caret")
install.packages("randomForest")
install.packages("glmnet")
install.packages("stringr")
install.packages('caTools')
install.packages("devtools")
install.packages("tictoc")
install.packages('varImp')
install.packages("tidyverse")
install.packages("shadowtext")
install.packages("scales")
library(scales)
#varimp
install.packages("vip")
library(vip)

library(ggplot2)
library(caTools)
library("stringr")
library(caret)
library(tidyverse)
library("readxl")
library(tictoc)
library(glmnet)
library(dplyr)
#for changing rownames to a column
library(tibble)
library(ranger)
library(gbm)
library(grid)
library(shadowtext)
library(forcats)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              DATA PREPARATION                            ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LOAD THE DATA 
lipgene_data <- read.csv("C:\\Users\\hoken\\surfdrive\\Shared\\Akiya_Hoken_2022\\LIPGENE data\\LIPGENE_glucmap_semi.csv",
                         header = T,
                         dec= ".",
                         sep = ",",
                         stringsAsFactors = F,
                         na.strings ="#LEEG!"
)

#Delete any duplicates, or rows containing no information (all and almost all NA's, which is dropout of 2,4,5,6)
# duplicated(lipgene_data) #found no duplicates
lipgene_data  <- filter(lipgene_data, !(lipgene_data$dropout %in% c(2,4,5,6)))

#Eliminate useless data (we eliminate "diet" since we have "dietgroup" with same info) (SI_sqrt_ch is nzv, thus we eliminate)
#also eliminating si_sqrt_ch because its fucking up bagImpute
eliminate <- c("phase","partnernr","center","subjectnr","diet","MF_LF","dropout","SI_sqrt_ch")

lipgene_data <- lipgene_data[, !(names(lipgene_data) %in% eliminate)] 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                PREPROCESSING                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#attempt for pre processing the entire data set (we only found one that had nzv)
# nzv <- nearZeroVar(lipgene_data)
# filteredLipgene <- lipgene_data[,-nzv]
set.seed(5)
imp <- caret::preProcess(lipgene_data, method = "bagImpute", k=5)
lipgene_data <- predict(imp,lipgene_data)

#save point to not keep running same 
save <- lipgene_data
lipgene_data <- save

#Function for creating a column with their skeletal muscle mass 
#CHANGE THE TARGET VARIABLE HERE, BUT NOT THE NAME SO U DONT HAVE TO CHANGE EVERYTHING!!
lipgene_data$smm_delta <- with(lipgene_data, 
                            (lipgene_data$weight_ch*0.065 + ((((lipgene_data$height*100)^2/(lipgene_data$impedance_post))-((lipgene_data$height*100)^2/(lipgene_data$impedance_post)))*0.2394))
                              )
#get rid of inf numbers from divinding 0
lipgene_data[sapply(lipgene_data, is.infinite)] <- 0

#get rid of weight related and impedance and HEIGHTTTTTT
eliminate_weights <- c("weight_pre","weight_mid","weight_ch","weightch_mid",
                       "weight_post","BMI_pre","BMI_post","BMI_ch","BMI_ln_ch","impedance_pre","impedance_post","height")
lipgene_data <- lipgene_data[, !(names(lipgene_data) %in% eliminate_weights)] 

#now for checking for correlated predictors
lipgeneCor <- cor(lipgene_data)
summary(lipgeneCor[upper.tri(lipgeneCor)])  #highest has 0.99755 correlation, which is bad!! thus we get rid of them
highcor <- findCorrelation(lipgeneCor, names=TRUE, cutoff = .85)
print(highcor)
# lipgene_data <- lipgene_data[,-highcor]
# idk <- cor(lipgene_data)
# summary(idk[upper.tri(idk)]) #now down to 160 variables!!

#attempt to delete pre/post as reduction in dimension?
lipgene_data <- lipgene_data %>% select(-(contains(
  c("_ch","_ln_","_kg_","_ltr_","_mean","_sqrt_","WHR","HOMA","LBM","PUFA_PCT_P",
    "_post","TRL_TG_pre","_sin_",
    "Total_physical_activity","apoB48_pre"))))

lipgeneCor <- cor(lipgene_data)
summary(lipgeneCor[upper.tri(lipgeneCor)]) #now down to 0.831209
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               DATA SPLITTING                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create empty dataframes (for each type of diet)
diet_A <- data.frame()
diet_A_M <- data.frame()
diet_A_F <- data.frame()
diet_B <- data.frame()
diet_B_M <- data.frame()
diet_B_F <- data.frame()
diet_C <- data.frame()
diet_C_M <- data.frame()
diet_C_F <- data.frame()
diet_D <- data.frame()
diet_D_M <- data.frame()
diet_D_F <- data.frame()

#split each diet groups and their sex (categorical)
#VERY INEFFICIENT CODE FOR NOW!!!!
for (i in 1:nrow(lipgene_data)) {         #for every observations,
  if (lipgene_data[i,]$dietgroup == 1){   #if the dietgroup is 1, then store it in diet_A
    diet_A = rbind(diet_A,lipgene_data[i,])
    if (lipgene_data[i,]$gender==0){
      diet_A_M = rbind(diet_A_M,lipgene_data[i,])
    }
    else{
      diet_A_F = rbind(diet_A_F,lipgene_data[i,])
    }
  }
  else if (lipgene_data[i,]$dietgroup == 2){
    diet_B = rbind(diet_B,lipgene_data[i,])
    if (lipgene_data[i,]$gender==0){
      diet_B_M = rbind(diet_B_M,lipgene_data[i,])
    }
    else{
      diet_B_F = rbind(diet_B_F,lipgene_data[i,])
    }
  }
  else if (lipgene_data[i,]$dietgroup == 3){
    diet_C = rbind(diet_C,lipgene_data[i,])
    if (lipgene_data[i,]$gender==0){
      diet_C_M = rbind(diet_C_M,lipgene_data[i,])
    }
    else{
      diet_C_F = rbind(diet_C_F,lipgene_data[i,])
    }
  }
  else if (lipgene_data[i,]$dietgroup == 4){
    diet_D = rbind(diet_D,lipgene_data[i,])
    if (lipgene_data[i,]$gender==0){
      diet_D_M = rbind(diet_D_M,lipgene_data[i,])
    }
    else{
      diet_D_F = rbind(diet_D_F,lipgene_data[i,])
    }
  }
}

#now remove info on diet groups and sex
eliminate_ds <- c("gender","dietgroup")
# AGAGIN, ABSOLUTELY HORRIFIC LOOKING CODE BUT ITS THERE FOR NOW
diet_A <- diet_A[,!(names(diet_A) %in% eliminate_ds)]
diet_A_F <- diet_A_F[, !(names(diet_A_F) %in% eliminate_ds)] 
diet_A_M <- diet_A_M[, !(names(diet_A_M) %in% eliminate_ds)]
diet_B <- diet_B[,!(names(diet_B) %in% eliminate_ds)]
diet_B_F <- diet_B_F[, !(names(diet_B_F) %in% eliminate_ds)] 
diet_B_M <- diet_B_M[, !(names(diet_B_M) %in% eliminate_ds)] 
diet_C <- diet_C[,!(names(diet_C) %in% eliminate_ds)]
diet_C_F <- diet_C_F[, !(names(diet_C_F) %in% eliminate_ds)] 
diet_C_M <- diet_C_M[, !(names(diet_C_M) %in% eliminate_ds)]
diet_D <- diet_D[,!(names(diet_D) %in% eliminate_ds)]
diet_D_F <- diet_D_F[, !(names(diet_D_F) %in% eliminate_ds)] 
diet_D_M <- diet_D_M[, !(names(diet_D_M) %in% eliminate_ds)] 

# apparently this doesnt work neither :(
# diets = list(diet_A,diet_A_F,diet_A_M,diet_B,diet_B_F,diet_B_M,
#              diet_C,diet_C_F,diet_C_M,diet_D,diet_D_F,diet_D_M)
# 
# eliminateDS <- function(diets,variables){
#   diets <- diets[,!(names(diets) %in% variables)]
# }
# 
# for(i in 1:length(diets)){
#   diets[[i]] = eliminateDS(diets[[i]],eliminate_ds)
# }

#SPLIT IS ACTUALLY SO SHIT 
# diets = split(lipgene_data, f = lipgene_data$dietgroup)
# diet_A<- as.data.frame(diets[1])
# diet_A_g <- split(diet_A,f = diet_A$X1.gender)
# diet_A_M <- as.data.frame(diet_A_g[1])
#

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##------------------------- MACHINE LEARNING MODELS-----------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MODELRUN, no longer in use!!!
# 
# #function for the models
# model_predict <- function(diet,model_name){
#   tic()
#   #reproducibility 
#   set.seed(5)
#  
#   # Get the number of observations
#   n_obs <- nrow(diet)
#   
#   # Identify row to split on: split
#   split <- round(n_obs * 0.75)
#   
#   # Create train
#   train <- diet[1:split, ]
#   
#   # Create test
#   test <- diet[(split+1):nrow(diet),]
#   
#   #split predictor variable to target variable, favourable with caret packages
#   x <- train[,colnames(train) != "smm_delta"]
#   y <- train$smm_delta
#   
#   #model attempt
#   model <- train(
#     x = x,
#     y = y,
#     preProcess = c("center","scale"),
#     method = model_name,
#     
#     # #for ranger varImp
#     # importance = 'impurity',
#     
#     trControl = trainControl(
#       method =  "cv", 
#       number = 10,
#       verboseIter = TRUE
#     )
#   )
#   predicted <- predict(model, test)
#   print(predicted)
#   plot(predicted, test$smm_delta, main = model_name)
#   abline(reg=lm(test$smm_delta ~ predicted))
#   toc()
#   return(model)
# }
# #create model objects that return model
# glmnet <- model_predict(diet_A,"glmnet")
# ranger <- model_predict(diet_D,"ranger")
# lm <- model_predict(diet_D,"lm")
# gbm <- model_predict(diet_D,"gbm")
# brnn <- model_predict(diet_D,"brnn")
# rf <- model_predict(diet_D,"rf")
# #gradient decent --> doesnt exist anymore, great!!
# # mlpSGD <- model_predict(diet_D,"mlpSGD")
# xgbLinear <- model_predict(diet_D,"xgbLinear")
# xgbTree<- model_predict(diet_D,"xgbTree")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ CONTROLLED, TEST/TRAIN DATASETS  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(5)
  
# Get the number of observations
n_obs <- nrow(diet_D)
  
# Identify row to split on: split
split <- round(n_obs * 0.75)
  
# Create train
train <- diet_D[1:split, ]
  
# Create test
test = diet_D[(split+1):nrow(diet_D),]
  
#split predictor variable to target variable, favourable with caret packages
x <- train[,colnames(train) != "smm_delta"]
y <- train$smm_delta
  
#my trainControl setting to make it a fair comparison
myControl <- trainControl(
  method =  "cv", 
  number = 10,
  verboseIter = TRUE
)

#adaptive resampling exapmle --> not going to use, since manual tuning is 
#more fun + less computational for models with large number of hyperparameters
# fitControl <- trainControl(
#   method = "adaptive_cv",
#   number = 3,
#   repeats =3,
#   adaptive = list(
#     min = 2,   #minimum number of resamples per hyperparmeter (default at 5)
#     alpha = 0.05,   #confidence level for removin g hyperparam (chanign this doenst influence that much)
#     method = "gls",   #either "gls" or "BT" - resampling method. gls for linear models, or BT for models with large number of hyperparameters
#     complete = TRUE   #when set to true, do we want to generate the resampling set?
#   ),
#   search = "random"
# )

#---------------------------Linear Regression Models----------------------------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                   GLMNET                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(5)
# custom tuning grid
myGrid <- expand.grid(
  #alpha =0,   #best fit
  alpha = seq(0,1, length =5),   #0 - pure ridge to 1- pure lasso.
  #lambda = 1 #bestfit
  lambda = seq(1,0.0001,length =100) #go from 0.0001 to 1, having 100 stops in between. (diet_C collapses at anything higher that 0.0875)
)

#Model object
tic() #time code execution
glmnet <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  method = "glmnet",
  tuneGrid = myGrid,
  trControl = myControl
)
#create predicted model object
p_glmnet <- predict(glmnet, test)
toc()
#plot the test data against the predicted, with best fit line
plot(p_glmnet, test$smm_delta, main = "glmnet")
abline(reg=lm(test$smm_delta ~ p_glmnet))
#out-ofsample rmse
rmse_glmnet <- postResample(pred = p_glmnet, obs = test$smm_delta)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    PLS                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(5)
# custom tuning grid
myGrid <- expand.grid(
  ncomp = c(1:45) #number of components to include in the model
)

#model
tic()
pls <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  mode="regression",
  method = "pls",
  tuneGrid = myGrid,
  trControl = myControl
)
#predict object
p_pls <- predict(pls, test)
toc()

#plot 
plot(p_pls, test$smm_delta, main = "pls")
abline(reg=lm(test$smm_delta ~ p_pls))

#out-ofsample rmse
rmse_pls <- postResample(pred = p_pls, obs = test$smm_delta)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    ENET                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(5)
#custom tuning grid
myGrid <- expand.grid(
  lambda = c(1e-5,1e-4,1e-3,1e-2,1e-1,0,1,10,100,1000),   #weight decay
  #lambda = 1e-5,  # best fit
  fraction = seq(0,1, length = 100)   #fraction of ratio between L1/L2 pentalty (diet_C, including 0 collapses for some reason, so set to 0.0001)
  #fraction = 0.01010101   #best fit
)

#model
tic()
enet <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  method = "enet",
  tuneGrid = myGrid,
  trControl = myControl
)

p_enet <- predict(enet, test)
toc()

plot(p_enet, test$smm_delta, main = "enet")
abline(reg=lm(test$smm_delta ~ p_enet))

rmse_enet <- postResample(pred = p_enet, obs = test$smm_delta)

#--------------------------Ensemble Regression Models---------------------------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    GBM                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(5)
#custom tuning grid
myGrid <- expand.grid(
  n.trees = seq(10,200, length = 50), #number of trees, from 10 to 200 having 50 stops
  interaction.depth = seq(1,10,length=10), #complexity of the trees
  shrinkage = c(0.1,0.01,0.001),  #learning rate
  n.minobsinnode = 10  #minimum number of training set samples in a node

  # #best tune 
  # n.trees = 10, 
  # interaction.depth = 5,
  # shrinkage = 0.1,  
  # n.minobsinnode = 10  
)

#model
tic()
gbm <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"), 
  method = "gbm",
  tuneGrid = myGrid,
  trControl = myControl
)

p_gbm <- predict(gbm, test)
toc()

plot(p_gbm, test$smm_delta, main = "gbm")
abline(reg=lm(test$smm_delta ~ p_gbm))

rmse_gbm <- postResample(pred = p_gbm, obs = test$smm_delta) #out of sample rmse


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                   RANGER                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(5)
#custom tuning grid
myGrid <- expand.grid(
  mtry = c(2,3,4,5,10,20),   #Number of variables to possibly split at in each node.
  min.node.size = c(10,20,30,40,50),    #Minimal node size
  splitrule = c("variance","maxstat","extratrees")
  # mtry = 10,   #best tune
  # min.node.size = 30,    #best tune
  # splitrule = "extratrees"   #best
)

#model
tic()
ranger <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  method = "ranger",
  importance = 'impurity',  #this is required for ranger model's varImp to work (why? idk)
  tuneGrid = myGrid,
  trControl = myControl
)

p_ranger <- predict(ranger, test)
toc()

plot(p_ranger, test$smm_delta, main = "ranger")
abline(reg=lm(test$smm_delta ~ p_ranger))

rmse_ranger <- postResample(pred = p_ranger, obs = test$smm_delta)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  TREEBAG                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#model 
set.seed(5)
tic()
treebag <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  method = "treebag",
  trControl = myControl
)

p_treebag <- predict(treebag, test)
toc()

plot(p_treebag, test$smm_delta, main = "treebag")
abline(reg=lm(test$smm_delta ~ p_treebag))

rmse_treebag <- postResample(pred = p_treebag, obs = test$smm_delta)

#-----------------------Neural Network Regression Models------------------------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    BRNN                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(5)
#custom tuning grid
myGrid <- expand.grid(
  neurons = c(1,2,3,4,5,7)   #Number of variables to possibly split at in each node.
)

#model
tic()
brnn <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  method = "brnn",
  tuneGrid = myGrid,
  trControl = myControl
)

p_brnn <- predict(brnn, test)
toc()

plot(p_brnn, test$smm_delta, main = "brnn")
abline(reg=lm(test$smm_delta ~ p_brnn))

rmse_brnn <- postResample(pred = p_brnn, obs = test$smm_delta)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    NNET                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(5)
#custom tuning grid
myGrid <- expand.grid(
  size = c(1,2,3,4,5,6,7,8,9,10),  #number of hidden layers
  decay = c(10,100,500,800)
  # size = 4,  #best tune
  # decay = 300
)

#model 
set.seed(5)
tic()
nnet <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  method = "nnet",
  linout = TRUE, #switch for linear output units. Default logistic output units.
  tuneGrid = myGrid,
  trControl = myControl
)

p_nnet <- predict(nnet, test)
toc()

plot(p_nnet, test$smm_delta, main = "nnet")
abline(reg=lm(test$smm_delta ~ p_nnet))

rmse_nnet <- postResample(pred = p_nnet, obs = test$smm_delta)

#--------------------------Stepwise Regression Models---------------------------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                 glmStepAIC                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(5)

#model
tic()
glmStepAIC <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  method = "glmStepAIC",
  trControl = myControl
)

p_glmStepAIC <- predict(glmStepAIC, test)
toc()

plot(p_glmStepAIC, test$smm_delta, main = "glmStepAIC")
abline(reg=lm(test$smm_delta ~ p_glmStepAIC))

rmse_glmStepAIC <- postResample(pred = p_glmStepAIC, obs = test$smm_delta)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  leapSeq                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#model
set.seed(5)

myGrid <- expand.grid(
  nvmax = c(2,3,4,5,10,15,20,30)
  #nvmax = 2
)

tic()
leapSeq <- train(
  x = x,
  y = y,
  preProcess = c("center","scale"),
  tuneGrid = myGrid,
  method = "leapSeq",
  trControl = myControl
)

p_leapSeq <- predict(leapSeq, test)
toc()

plot(p_leapSeq, test$smm_delta, main = "leapSeq")
abline(reg=lm(test$smm_delta ~ p_leapSeq))

rmse_leapSeq <- postResample(pred = p_leapSeq, obs = test$smm_delta)


#----------------------------- MODEL ATTEMPTS END-------------------------------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##--------------------------------- RESULTS-------------------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ IN-SAMPLE RMSE  ----
##~~~~~~~~~~~~~~~~~~~~~~~~
#get summary of comparison between each models, IN-SAMPLE RMSE'S 
model_list <- list(glmnet,pls,enet,gbm,ranger,treebag,brnn,nnet,glmStepAIC,leapSeq)
model_names <- c("glmnet","pls","enet","gbm","ranger","treebag","brnn","nnet","glmStepAIC","leapSeq")
results <- summary(resamples(model_list))
#store the in-sample-rmse results, 
is_rmse_results <- as.data.frame(results$statistics[2])
#and rename the row titles to its respective model names
rownames(is_rmse_results) <- model_names
is_rmse_results

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ OUT-OF-SAMPLE RMSE  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#store all out of sample rsme's as vector
os_rmse <- c(rmse_glmnet[1],rmse_pls[1],rmse_enet[1],rmse_gbm[1],rmse_ranger[1],
             rmse_treebag[1],rmse_brnn[1],rmse_nnet[1],rmse_glmStepAIC[1],rmse_leapSeq[1]) 
#make empty dataframe
os <- data.frame(row.names = model_names, os_rmse)
#print the lowest rmse and its name
os
print(os[which(os$os_rmse == min(os)),,drop=FALSE])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Did the model overfit?  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#If the training set performed better than the test set, then the model is overfitted
overfit <- data.frame(row.names = model_names,  os$os_rmse-is_rmse_results$RMSE.Mean) #test error - train error = overfit (positive value means overfitting)
#rename column
colnames(overfit) <- "overfit"
overfit

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ VARIABLE IMPORTANCE  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#apply varimp function to all elements of the model list
list_varimp <- lapply(model_list,varImp)
#store all this in one dataframe
df_varimp <- data.frame(rows=43) #hard-coded 43 just because its convenient

for(i in 1:length(list_varimp)){ #for the length of the list, get the varimp list, then bind onto the dataframe
  x <- list_varimp[[i]][1]
  df_varimp <- cbind(df_varimp,x)
}

#remove weird column named "rows"
df_varimp = df_varimp[,!names(df_varimp) %in% "rows"]
#manually replace the varimp for enet, brnn, and the two stepwise because varImp() is not supported on these models
df_varimp$Overall.2 = filterVarImp(test,p_enet)[-1,] #[-1,] because the overall title in the way 
df_varimp$Overall.6 = filterVarImp(test,p_brnn)[-1,]
df_varimp$Overall.8 = filterVarImp(test,p_glmStepAIC)[-1,]
df_varimp$Overall.9 = filterVarImp(test,p_leapSeq)[-1,]

#fit the names to its corresponding model
colnames(df_varimp) <- model_names

#add the row numbers up to see the total importance? (maybe irrelevant)
df_varimp = as.data.frame(scale(df_varimp))
df_varimp$sum_scaled <- rowSums(df_varimp)


# #EXPORT ALL RELEVANT DATA AS CSV
# #in sample rmse
# write.csv(is_rmse_results,"C:\\Users\\hoken\\surfdrive\\Shared\\Akiya_Hoken_2022\\2B_is_rmse.csv")
# #ou-of-sample rmse
# write.csv(os,"C:\\Users\\hoken\\surfdrive\\Shared\\Akiya_Hoken_2022\\2B_os_rmse.csv")
# #does it overfit
# write.csv(overfit,"C:\\Users\\hoken\\surfdrive\\Shared\\Akiya_Hoken_2022\\2B_overfit.csv")
# #variable importance
# write.csv(df_varimp,"C:\\Users\\hoken\\surfdrive\\Shared\\Akiya_Hoken_2022\\2B_varimp.csv")

#=B2*(1-(SQRT(1-((1.96*SQRT(2))/(SQRT(124))))))  for excel calculating positive error rmse
#=B2*((SQRT(1+((1.96*SQRT(2))/(SQRT(124)))))-1)  for negative error rmse


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##------------------------------- VISUAL GRAPHS---------------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~
##  ~ RMSE  ----
##~~~~~~~~~~~~~~

#calculate the errors of rmse
error_function_pos <- function(val){
    return(val*(1-(sqrt(1-((1.96*sqrt(2))/(sqrt(124)))))))
}

error_function_neg <- function(val){
  return(val*((sqrt(1+((1.96*sqrt(2))/(sqrt(124)))))-1))
}

#assign errors
pos_err <- unlist(lapply(os$os_rmse,error_function_pos))
neg_err <- unlist(lapply(os$os_rmse,error_function_neg))

#plot the horizontal graph with error bars
os %>%
  ggplot( aes(reorder(row.names(os), -os_rmse, sum), os_rmse)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  geom_errorbar( aes(x=row.names(os), ymin=os_rmse-neg_err, ymax=os_rmse+pos_err), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  geom_text(aes(y=0,label = round(os_rmse,5)), hjust = -0.1)+
  coord_flip() +
  xlab("") +
  theme_bw() +
  labs(title = "Diet C - RMSE per model") + 
  ylab("RMSE Value") + 
  xlab("Model name")+ 
  theme(
  plot.title = element_text(color="red", size=16, face="bold.italic"),
  axis.title.x = element_text(color="blue", size=12, face="bold"),
  axis.title.y = element_text(color="blue", size=12, face="bold")
)


##~~~~~~~~~~~~~~~~
##  ~ VarImp  ----
##~~~~~~~~~~~~~~~~

head(df_varimp%>% arrange(desc(sum_scaled)),10) %>%
  ggplot( aes(reorder(row.names(head(df_varimp%>% arrange(desc(sum_scaled)),10)), -sum_scaled, sum), sum_scaled)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  geom_text(aes(label = round(sum_scaled,3)), vjust = -0.2) +
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  labs(title = "Diet C - Scaled Variable Importance Values ") + 
  ylab("Sum of Scaled Variable Importance") + 
  xlab("Variable Name") + theme(
    plot.title = element_text(color="red", size=16, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=12, face="bold"),
    axis.title.y = element_text(color="blue", size=12, face="bold")
  )


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##------------------------- COMBINED RESULT SECTION-----------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#...............................................................................
#                                                                              .
#  DO NOT RUN THIS SECTION ALL AT ONCE, this section should only be ran per    .
#  diet as commented - it won't work otherwise!! First run the models with     .
#  one diet, THEN come here to run ONLY the corresponding diet codes. Do this  .
#  with all the diets, then run the final "combined" code to get combined      .
#  results                                                                     .
#                                                                              .
#...............................................................................

##~~~~~~~~~~~~~~~~
##  ~ VARIMP  ----
##~~~~~~~~~~~~~~~~

#method to store different diet results
#RUN THIS PER DIET ONLY, do not run all at once

#DIET A
dietA_varimp <- head(df_varimp %>% arrange(desc(sum_scaled)),10)
dietA_varimp <- dietA_varimp[,11,drop=FALSE]

#DIET B
dietB_varimp <- head(df_varimp %>% arrange(desc(sum_scaled)),10)
dietB_varimp <- dietB_varimp[,11,drop=FALSE]

#DIET C
dietC_varimp <- head(df_varimp %>% arrange(desc(sum_scaled)),10)
dietC_varimp <- dietC_varimp[,11,drop=FALSE]

#DIET D
dietD_varimp <- head(df_varimp %>% arrange(desc(sum_scaled)),10)
dietD_varimp <- dietD_varimp[,11,drop=FALSE]

#graphing function for convenience
varimp_plotter <- function(diet_varimp,diet){
  diet_varimp%>% arrange(desc(summ)) %>%
    ggplot(aes(reorder(row.names(diet_varimp%>% arrange(desc(summ))), -summ, FUN=sum), summ)) +
    geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
    geom_text(aes(label = round(summ,3)), vjust = -0.2) +
    scale_x_discrete(guide = guide_axis(n.dodge=2))+
    labs(title = paste(diet, " - Scaled Variable Importance Values")) + 
    ylab("Sum of Scaled Variable Importance") + 
    xlab("Variable Name") + theme(
      plot.title = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="black", size=12, face="bold"),
      axis.title.y = element_text(color="black", size=12, face="bold")
    )
} 

varimp_plotter(dietA_varimp, "Diet A")
varimp_plotter(dietB_varimp, "Diet B")
varimp_plotter(dietC_varimp, "Diet C")
varimp_plotter(dietD_varimp, "Diet D")

##~~~~~~~~~~~~~~~~~~~~~~~
##  ~ RMSE COMBINED  ----
##~~~~~~~~~~~~~~~~~~~~~~~
#ERRORS (I dont know how but accidentally capitalised the "N" in Neg_err)
#DIET A
pos_err_A <- pos_err
Neg_err_A <- neg_err
#DIET B
pos_err_B <- pos_err
Neg_err_B <- neg_err
#DIET C
pos_err_C <- pos_err
Neg_err_C <- neg_err
#DIET D
pos_err_D <- pos_err
Neg_err_D <- neg_err

#RMSE's
#DIET A
os_A <- os
#DIET B
os_B <- os
#DIET C
os_C <- os
#DIET D
os_D <- os

#### QUICK GRAPH GENERATING FUNCTION FOR EASE OF USE####
rmse_plotter <- function(os,pos,neg,model){
  os %>%
    ggplot( aes(reorder(row.names(os), -os_rmse, sum), os_rmse)) +
    geom_bar(stat="identity", fill="#56B4E9", alpha=.6, width=.4) +
    geom_errorbar( aes(x=row.names(os), ymin=os_rmse-neg, ymax=os_rmse+pos), width=0.3, colour="black", alpha=0.9, size=1) +
    geom_text(aes(y=0,label = round(os_rmse,5)), hjust = -0.1)+
    coord_flip() +
    xlab("") +
    theme_bw() +
    labs(title = paste(model," - RMSE per model")) + 
    ylab("RMSE Value") + 
    xlab("Model name")+ 
    theme(
      plot.title = element_text(color="black", size=16, face="bold.italic"),
      axis.title.x = element_text(color="black", size=12, face="bold"),
      axis.title.y = element_text(color="black", size=12, face="bold")
    )
}

rmse_plotter(os_A,pos_err_A,Neg_err_A,"Diet A")
rmse_plotter(os_B,pos_err_B,Neg_err_B,"Diet B")
rmse_plotter(os_C,pos_err_C,Neg_err_C,"Diet C")
rmse_plotter(os_D,pos_err_D,Neg_err_D,"Diet D")
rmse_plotter(os_all,av_pos_err,av_neg_err,"Average")




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##combined section, run this AFTER you have ran the above codes for each diet----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~
##  ~ rmse  ----
##~~~~~~~~~~~~~~
#average errors
av_pos_err <- (pos_err_A+pos_err_B+pos_err_C+pos_err_D)/4
av_neg_err <-(Neg_err_A+Neg_err_B+Neg_err_C+Neg_err_D)/4

#average rmse per model
os_all <- (os_A + os_B + os_C + os_D)/4

#plot this
os_all %>%
  ggplot( aes(reorder(row.names(os_all), -os_rmse, sum), os_rmse)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  geom_errorbar( aes(x=row.names(os_all), ymin=os_rmse-av_neg_err, ymax=os_rmse+av_pos_err), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  geom_text(aes(y=0,label = round(os_rmse,5)), hjust = -0.1)+
  coord_flip() +
  xlab("") +
  theme_bw() +
  labs(title = "Average RMSE of All Diet Types") + 
  ylab("RMSE Value") + 
  xlab("Model name")+ 
  theme(
    plot.title = element_text(color="red", size=16, face="bold.italic"),
    axis.title.x = element_text(color="blue", size=12, face="bold"),
    axis.title.y = element_text(color="blue", size=12, face="bold")
  )

##~~~~~~~~~~~~~~~~
##  ~ Varimp  ----
##~~~~~~~~~~~~~~~~
#merge all varimps
# DietAB_varimp <- merge(dietA_varimp,dietB_varimp,all = TRUE,by=0) #merge the two by names, if new variable is present then add new row
# rownames(DietAB_varimp) <- DietAB_varimp[,1] #reformatting because new column with just the names were made
# DietAB_varimp[,1] <- NULL
# DietCD_varimp <- merge(dietC_varimp,dietD_varimp,all = TRUE,by=0)
# rownames(DietCD_varimp) <- DietCD_varimp[,1]
# DietCD_varimp[,1] <- NULL
# Dietall_varimp <- merge(DietAB_varimp,DietCD_varimp,all = TRUE,by=0)
# rownames(Dietall_varimp) <- Dietall_varimp[,1]
# Dietall_varimp[,1] <- NULL
# colnames(Dietall_varimp) <- c("Diet A","Diet B","Diet C","Diet D") #rename the column names because it was weird
# 
# #now reorder the variable names according to the total sum of each varimp scores
# #all na is now 0
# Dietall_varimp[is.na(Dietall_varimp)] <- 0
# Dietall_varimp$sum <- rowSums(Dietall_varimp)
# 
# #reorder
# Dietall_varimp <- Dietall_varimp %>% arrange(desc(sum))
# Dietall_varimp <- add_column(Dietall_varimp, row.names(Dietall_varimp), .before = 1)
# colnames(Dietall_varimp) <- c("Variables","Diet A","Diet B","Diet C","Diet D","sum") #little redundant but we move for now
# 
# Dietall_varimp <- data.frame(t(Dietall_varimp[-1]))
# barplot(as.matrix(Dietall_varimp))

#attempt with ggplot!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#DO NOT TOUCH THE FOLLOWING LINES UNTIL NOTIFIED AFTER FIRST RUN
dietA_varimp$diet <- "Diet A"
dietB_varimp$diet <- "Diet B"
dietC_varimp$diet <- "Diet C"
dietD_varimp$diet <- "Diet D"
dietA_varimp <- add_column(dietA_varimp, row.names(dietA_varimp), .before = 1)
dietB_varimp <- add_column(dietB_varimp, row.names(dietB_varimp), .before = 1)
dietC_varimp <- add_column(dietC_varimp, row.names(dietC_varimp), .before = 1)
dietD_varimp <- add_column(dietD_varimp, row.names(dietD_varimp), .before = 1)
rownames(dietA_varimp) <- NULL
rownames(dietB_varimp) <- NULL
rownames(dietC_varimp) <- NULL
rownames(dietD_varimp) <- NULL
colnames(dietA_varimp) <- c("Variables","summ","diet") #little redundant but we move for now
colnames(dietB_varimp) <- c("Variables","summ","diet") 
colnames(dietC_varimp) <- c("Variables","summ","diet") 
colnames(dietD_varimp) <- c("Variables","summ","diet") 
#UNTIL HERE, DO NOT RUN. BELOW THIS U CAN RUN AGAIN


alldiets <- bind_rows(dietA_varimp,dietB_varimp,dietC_varimp,dietD_varimp)
# alldiets <- add_column(alldiets, row.names(alldiets), .before = 1)
# colnames(alldiets) <- c("Variables","sum","diet") #little redundant but we move for now
# rownames(alldiets) <- NULL

alldiets <- alldiets %>%
  group_by(Variables) %>%
  mutate(label_y = cumsum(sum))

alldiets <- alldiets %>% arrange(desc(label_y))
alldiets$Variables <- factor(alldiets$Variables)

#same code again to reorder the variable numbers too apparently
alldiets <- alldiets %>%
  group_by(Variables) %>%
  mutate(label_y = cumsum(sum))

#plot
ggplot(alldiets, aes(y = fct_rev(fct_inorder(Variables)), x = sum, fill = diet)) +
  geom_col() +
  geom_text(aes(x = label_y-0.5*sum, label = round(sum,2)),  colour = "white") +
  labs(title = "Combined variable importance values per diet type") + 
  ylab("Variable Names") + 
  xlab("Scaled Importance Value")+ 
  theme(
    plot.title = element_text(color="black", size=16, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
  
#### RESULTS END ####             



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##---------------------- SECOND ROUND OF MODEL TESTING?-------------------------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Second round ####
#now pic the 10 most influential variables from the previous result, then get rid of other variables to see the model fit!
imp_vars <- c("BF_pct_pre", "impedance_pre", "BW_pct_pre", 
              "hip_pre", "CHO_en_pre","smm_delta")
diet_A_2 <- diet_A[,names(diet_A) %in% imp_vars]
diet_A_M_2 <- diet_A_M[,names(diet_A_M) %in% imp_vars]
diet_A_F_2 <- diet_A_F[,names(diet_A_F) %in% imp_vars]
diet_B_2 <- diet_B[,names(diet_B) %in% imp_vars]
diet_B_M_2 <-diet_B_M[,names(diet_B_M) %in% imp_vars] 
diet_B_F_2 <- diet_B_F[,names(diet_B_F) %in% imp_vars]
diet_C <- diet_C[,names(diet_C) %in% imp_vars]
diet_C_M <- diet_C_M[,names(diet_C_M) %in% imp_vars]
diet_C_F <- diet_C_F[,names(diet_C_F) %in% imp_vars]
diet_D <-diet_D[,names(diet_D) %in% imp_vars]
diet_D_M <- diet_D_M[,names(diet_D_M) %in% imp_vars]
diet_D_F <- diet_D_F[,names(diet_D_F) %in% imp_vars]

#------------ MODEL ATTEMPTS ------------#

set.seed(5)

# Get the number of observations
n_obs2 <- nrow(diet_A_2)

# Identify row to split on: split
split2 <- round(n_obs2 * 0.75)

# Create train
train2 <- diet_A_2[1:split2, ]

# Create test
test2 = diet_A_2[(split2+1):nrow(diet_A_2),]

#split predictor variable to target variable, favourable with caret packages
x2 <- train2[,colnames(train2) != "smm_delta"]
y2 <- train2$smm_delta

#my trainControl setting to make it a fair comparison
myControl <- trainControl(
  method =  "cv", 
  number = 10,
  verboseIter = TRUE
)
#------------linear regression models------------#

####GLMNET2 ATTEMPT####
set.seed(5)
# custom tuning grid
myGrid <- expand.grid(
  alpha = seq(0,1, length =5),   #0 - pure ridge to 1- pure lasso.
  lambda = seq(1,0.0001,length =100) #go from 0.0001 to 1, having 100 stops in between
)

#model attempt
tic() #time code execution
glmnet2 <- train(
  x = x2,
  y = y2,
  preProcess = c("center","scale"),
  method = "glmnet",
  tuneGrid = myGrid,
  trControl = myControl
)
#create model object
p_glmnet2 <- predict(glmnet2, test2)
toc()
#plot the test data against the predicted, with best fit line
plot(p_glmnet2, test2$smm_delta, main = "glmnet2")
abline(reg=lm(test2$smm_delta ~ p_glmnet2))
#out-ofsample rmse
rmse_glmnet2 <- postResample(pred = p_glmnet2, obs = test2$smm_delta)
####GLMNET END####

####nnet2 ATTEMPT####
set.seed(5)
#custom tuning grid
myGrid <- expand.grid(
  size = c(1,2,3,4,5),  #number of hidden layers
  decay = c(1e-3, 1e-2, 1e-1, 0.0, 1.0,10,100)
)

#model attempt
set.seed(5)
tic()
nnet2 <- train(
  x = x2,
  y = y2,
  preProcess = c("center","scale"),
  method = "nnet",
  linout = TRUE, #switch to linear output units. Default logistic output units.
  tuneGrid = myGrid,
  trControl = myControl
)

p_nnet2 <- predict(nnet2, test2)
toc()

plot(p_nnet2, test2$smm_delta, main = "nnet2")
abline(reg=lm(test2$smm_delta ~ p_nnet2))

rmse_nnet2 <- postResample(pred = p_nnet2, obs = test2$smm_delta)
#print(min(nnet$results[3]))


####nnetEND####