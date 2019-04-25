library(msEmpiRe)

#for evaluation at the end
library(PerfMeas)
#to extract featureData from data, also used for evaluation
library(Biobase)


#########################
# DIFFERENTIAL ANALYSIS #
#########################

f <- system.file("extdata", "c1_c3.data", package = "msEmpiRe")
p <- system.file("extdata", "c1_c3.pdata", package = "msEmpiRe")

#loading data from installed data sets. The dataset is the yeast spike-in benchmarking data of O'Connell et al., as discussed in the paper
data <- msEmpiRe::read.standard(f, p,
                                prot.id.generator=function(pep) unlist(strsplit(pep, "\\."))[1],
                                signal_pattern="c.*rep.*")


#extract the first two conditions
conditions <- extract_conditions(data)
conditions <- conditions[, c(1,2)]


#removing peptides that are detected in less than 2 samples per condition
data <- msEmpiRe::filter_detection_rate(data, condition=conditions)

#normalize, system.time is just to measure the time of processing
system.time(data <- msEmpiRe::normalize(data))

#main analysis, system.time is just to measure the time of processing
system.time(result <- de.ana(data))


##############
# EVALUATION #
##############
#as the data we have analysed is a benchmarking set, with known fold changes applied to a subset of the proteins
#(see paper for more detail),
#we can now evaluate the differential analysis, i.e. extract measures like precision, recall etc.
#all proteins are labelled, wether they are actually changing ("true"), or not ("false") 

#extracting the labels,aggregate per protein
tmp <- fData(data)
tmp$is_true <- sapply(tmp$is_true, function(val)
  {if(val=="true"){return(1)}
  return(0)})
print(head(tmp))


tmp <- aggregate(is_true ~ prot.id, tmp, sum)
tmp$is_true <- sapply(tmp$is_true, function(val)
  {
  if(val > 0)
  {
    return(1)
  }
  return(0)
})

#labels must have same order as result
rownames(tmp) <- tmp$prot.id
tmp <- tmp[as.character(result$prot.id), ]

result$label <- tmp$is_true

min_fc <- 0.35
scores <- cbind(colnames(result)[grep("p.val", colnames(result))], colnames(result)[grep("p.adj", colnames(result))])

#print out the scores of the evaluation
apply(scores, 1, function(score)
{
  print(score)
  tp <- sum( (result[, score[2]] <= 0.05) & (result$label==1) & (abs(result$log2FC) >= min_fc))
  real <- sum( (result$label==1))
  called <- sum( (result[, score[2]] <= 0.05) & (abs(result$log2FC) >= min_fc))

  print(sprintf("PREC: %g", tp / called))
  print(sprintf("REC: %g", tp / real))

  auroc <- AUC.single(1.0 - result[, score[1]], result$label)
  preccc <- list(precision.at.all.recall.levels(1.0 - result[, score[1]], result$label))
  auprc <- AUPRC(preccc)
  print(sprintf("AUROC: %g", auroc))
  print(sprintf("AUPRC: %g", auprc))
  print("----------")
})
