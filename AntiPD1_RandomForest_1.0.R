#install.packages("randomForest")
library(randomForest)

#require(pROC)
library(pROC)
library(graphics)

#use aao3290_Matson_SM-Table-S2.xlsx (from Matson et al., 2018: 42 patients with 63 bacteria relative abundances)
#rawOTU <- read.csv(file.choose(), header = T)
rawOTU <- read.csv("~/Desktop/Amy_HighSchool/Amy_10th/Resarch/62OTU_RACount_table.csv", header = T)
name = rawOTU$TaxaOTU
#transpose, change original column + column name to row + row name
OTU=as.data.frame(t(rawOTU[,-1]))
colnames(OTU)=name
OTU$ID=row.names(OTU)

#read in metadata containing response/nonresponse info (from Matson et al., 2018: response status of 41 patients ###one patient didn't have metadata)
#Meta <- read.csv(file.choose(), header = T)
#change directory for yourself
Meta <- read.csv("~/Desktop/Amy_HighSchool/Amy_10th/Resarch/AntiPD1_Metadata.csv", header = T)
Metaresponses = Meta[-1,1:2]
colnames(Metaresponses)=c("ID", "Response")
Metaresponses = Metaresponses[-(42:44),]

#merge Metaresponses with OTU
#change directory for yourself
OTU_Meta = merge(Metaresponses,OTU, by = "ID")
#> dim(OTU_Meta)
#[1] 41 65   ######63 bacteria and 41 patients

####### Sort the features based on their averaged mean decreased Gini score, which is averaged from total 63 leave-one-out test.

#make global data frame for importance scores (gini scores)
OTUginiAccuracy = data.frame("bacteria" = colnames(OTU_Meta)[-(1:2)])

##run model on all training sets (leave one out)
## takes validation sets using all 41 patients (1 patient per set "leave one out"), iterate through all 41
## training sets composed of 40 patients (leftoever, used to get working model)
for (i in 1:nrow(OTU_Meta[,-1])) {
  
  ValidSet = OTU_Meta[i,-1]
  TrainSet = OTU_Meta[-i,-1]
  # factorize the response variable
  TrainSet$Response <- factor(TrainSet$Response)
  ValidSet$Response <- factor(ValidSet$Response)
  
  # Create a Random Forest model with default parameters
  OTU_Model1 <- randomForest(Response ~ ., data = TrainSet, importance = TRUE)
  
  #finds the gini scores using the importance function
  OTUginiAccuracyTemp = as.data.frame(importance(OTU_Model1, type=2))
  
  ##record the gini scores into the global data frame above
  OTUginiAccuracy = cbind(OTUginiAccuracy, OTUginiAccuracyTemp$MeanDecreaseGini)
}

#average gini scores of each bacteria (collected from 41 iterations of training sets)
meanAccuracy_OTU=rowMeans(OTUginiAccuracy[,-1], na.rm = FALSE, dims = 1)
meanAccuracy_OTUtable$mean_decreaseGigi = meanAccuracy_OTU

#sort gini scores and matching bacteria names, obtain list of ranked bacteria based on importance
meanAccuracy_OTUtableSort = meanAccuracy_OTUtable[order(-meanAccuracy_OTUtable$mean_decreaseGigi),]
#### feature list done


######### calculate accuracy of the prediction using bacteria sets from the ranked list
###each round, carry out leave-one-out test
#global list for each "rank"
##"rank" = each bacteria set used (e.g. first set, 1 top bacteria; second set, 2 top bacteria; 28th set; 28 top bacteria)
EachRankAccuracy_Results = list()
RankAccuracyRatios = list()

#iteration; recursively adds the next bacteria to the "set" in each rank, starting from the first (most important) bacteria
### finds the PREDICTION ACCURACY based on VALIDATION SETS
for (j in 1:nrow(meanAccuracy_OTUtableSort)) {
  ##parse bacteria name -- allows it to stay as a data frame even when only 1 starting bacteria
  Feature = c("Response", as.vector(meanAccuracy_OTUtableSort$`names(meanAccuracy_OTU)`[1:j]))
  
  #add metadata back (response status)
  NewOTU_Meta = as.data.frame(OTU_Meta[,colnames(OTU_Meta) %in% Feature])
  
  summary = list()
  
  # use cross validation to find the prediction accuracy
  for (i in 1:nrow(NewOTU_Meta)) {
    ### same leave one out sets as before
    ValidSet = NewOTU_Meta[i,]
    TrainSet = NewOTU_Meta[-i,]
    #reassign factors
    TrainSet$Response <- factor(TrainSet$Response)
    ValidSet$Response <- factor(ValidSet$Response)
    
    # Create a Random Forest model with default parameters
    NewOTU_Model <- randomForest(Response ~ ., data = TrainSet, importance = TRUE)
    
    predValid <- predict(NewOTU_Model, ValidSet, type = "prob")
    predValid_class <- predict(NewOTU_Model, ValidSet, type = "class")
    ## Checking classification accuracy
    
    #some crude if else statements: checks to see if the prediction was accurate (correct +1, not correct +0)
    if(as.character(ValidSet[,1]) == "NonResponder") {
      Nonresponder_label = 1
    }else{
      Nonresponder_label = 0
    }
    if(as.character(predValid_class) == as.character(ValidSet[,1])){
      isCorrect = 1
    }else{
      isCorrect = 0
    }
    
    #### Prediction rate
    summary[[i]] = as.matrix(t(c(i, j, NewOTU_Model$err.rate[nrow(NewOTU_Model$err.rate),1],predValid[1],predValid[2], Nonresponder_label, as.character(ValidSet[,1]), isCorrect)))
    
  } #end of leave-one-out
  EachRankAccuracy_Results[[j]] = data.frame(t(sapply(summary, "[")))
  
  ##sum of isCorrect (correct prediction total)
  correctTotal = sum(as.numeric(as.character(EachRankAccuracy_Results[[j]][,7])))
  #find fraction out of total patients possible
  correctRatio = correctTotal/41
  
  #error average
  errorAverage = mean(as.numeric(as.character(EachRankAccuracy_Results[[j]][,3])))
  RankAccuracyRatios[[j]] = as.matrix(t(c(j,correctRatio, errorAverage)))
  
}
###### a summary dataframe!
### showcases: sample#, rank (# of bacteria from top), OOB error rate, nonresponder prediction rate, responder prediction rate, nonresponder label (if guessed nonresponder), the true responder status, and whether it was correct (1 or 0)
totalsummary = data.frame(do.call(rbind,EachRankAccuracy_Results))
colnames(totalsummary) = c("sample", "rank", "OOB", "Predict_Nonresponder_rate","Predict_Responder_rate", "Nonresponder_label", "Truth", "Correctness")


######### Calculate AUC

HighestRank_allRow <- max(as.numeric(totalsummary$rank))
auc_res <- list()

#iteration for each set of bacteria -- finds AUC for each iteration
for (NumberOfFeature in 1:HighestRank_allRow) {
  tryCatch({
    RankData <- totalsummary[totalsummary$rank == NumberOfFeature,]
    
    AUC_Area= auc(roc(as.factor(RankData$Nonresponder_label), as.numeric(as.vector(RankData$Predict_Nonresponder_rate))))
    
    auc_res[[NumberOfFeature]] <- data.frame(RankData[1,1:2], NumberOfFeature, AUC_Area)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
####data frame of auc
##features: sample#, rank, NumberOfFeature (same thing as rank), AUC Area
auc = data.frame(do.call(rbind,auc_res))

##### plot an ROC curve
# install.packages("pROC")
# library(pROC)
#png("~/Desktop/Amy_HighSchool/Amy_10th/Resarch/PD1_melanoma_Rank27_auc_otu_curve.png",width = 455,height = 435)
ggroc(roc(as.factor(RankData$Nonresponder_label), as.numeric(as.vector(RankData$Predict_Nonresponder_rate))), alpha = 0.5, colour = "red", linetype = 1, size = 2)
#dev.off()

#####old code: creates tables on the ranks sorted by percent correct and error rates given
# #by isCorrect number
# totalRatiosummary_SortByisCorrect = totalRatiosummary[order(-totalRatiosummary[,2]),]
# 
# #by error rates
# totalRatiosummary_SortByError = totalRatiosummary[order(totalRatiosummary[,3]),]
# 
# colnames(totalRatiosummary_SortByisCorrect) = c("Rank", "CorrectRatio_sorted","Error")
# colnames(totalRatiosummary_SortByError) = c("Rank", "CorrectRatio","Error_sorted")
# colnames(totalsummary) = c("Leave_one_out_Position", "Rank","Error","Truth","PredictedResults","isCorrect")
# 

