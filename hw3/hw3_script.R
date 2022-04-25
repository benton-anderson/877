# Problem 1 
# 1. Two group comparison: Apply EBSeq to the PBMC data to determine the number of genes differentially expressed (DE) between the two conditions.
#        1. How many genes are DE at false discovery rate (FDR) 5% and 10%? (3 pts)
#                3001 at 5% and 3396 at 10%
#        2. What are the estimated library size factors? What information do they provide? (3 pts)
#        3. Does the model fit well as assessed via diagnostics? Please show diagnostic plots and discuss. (3 pts)

library(EBSeq)
library(dplyr)

pbmc <- read.csv(file="C:\\Users\\bjanderson23\\Desktop\\classes\\877-stats\\hw3\\PBMC_COVID_NormExp.csv")

common_genes <- intersect(rownames(a549), rownames(nhbe))

genemat <- data(GeneMat)
GeneMat

conditions <- as.factor(c('C1', 'C1', 'C1', 'C2', 'C2', 'C2'))

Sizes <- MedianNorm(Data = pbmc[,-1])
Sizes

# Matrix is like a numpy array, all values must be the same dtype. 
# So if you convert a df with row names into matrix, you have to 
pbmc_matrix <- as.matrix(pbmc[, -1])
rownames(pbmc_matrix) <- pbmc[, 1]
pbmc_matrix

# pbmc_matrix <- as.matrix(pbmc)

EBOut2group <- EBTest(Data=pbmc_matrix, Conditions = conditions, sizeFactors = Sizes, maxround = 5)

EBOut2group

EBDERes5 <- GetDEResults(EBOut2group, FDR=0.05)
EBDERes10 <- GetDEResults(EBOut2group, FDR=0.1)
str(EBDERes5$DEfound)

head(EBDERes5$PPMat, 10)
order(EBDERes5$PPMat)

head(EBDERes5, 10)
order(EBDERes5)

write.csv(x = EBDERes5$PPMat, file="C:\\Users\\bjanderson23\\Desktop\\classes\\877-stats\\hw4\\fdr5.csv")

str(EBDERes10$DEfound)
head(EBDERes5$Status)

EBDERes5$PPMat %>% arrange(PPDE)
order(EBDERes5$PPMat[, 2])
colnames(EBDERes5$PPMat[, 0])
# PPMat = Posterior Probabilities Matrix
# PPMat contains PPEE and PPDE 
#  EE = Equivalently expressed, DE = Differentially Expressed
head(EBDERes5$PPMat, 5)

str(EBDERes5$Status)

qn <- QuantileNorm(GeneMat, 0.99)


apply(GeneMat, 2, median) / 371.5
MedianNorm(GeneMat)

# 3. Show diagnostic plots

GeneFC <- PostFC(EBOut)
str(GeneFC)

PlotPostVsRawFC(EBOut, GeneFC)

# The y-axis is the fold change of the raw data. Posterior FC is the adjusted fold change 
# There are two lines that deviate from y=x, and they have relatively high ranks (given by green color).
# The tutorial vignette states that the posterior FC tends to shrink genes with low expressions (small rank).
# But we see that relatively high ranked genes are shifted. I don't use R, so I can't dig into why that is and which genes are deviating. 

par(mfrow=c(1,2))
QQP(EBOut)
# The data points lie near to the y=x line, which suggests that using a Beta prior is appropriate. 
# There is a mild deviation on Condition 1 plot near 0.8, but I'm unsure whether that suggests any problem. 

DenNHist(EBOut)
# The fitted density closely matches the empirical distribution, suggesting that our distribution
# The fitted line refers to the alpha and beta
# parameters of the Beta distribution, given at the bottom of the plot. 
# If there is a large deviation from the fitted line, what does that say about our experiment? Bad experimental design? Bad sequencing run?
# What does a poorly fitted empirical distribution look like? 


## Four group comparison

nhbe <- read.csv(file='C:\\Users\\bjanderson23\\Desktop\\classes\\877-stats\\hw3\\NHBE_COVID_NormExp.csv', row.names=1)
length(nhbe[, 1])
a549 <- read.csv(file = "C:\\Users\\bjanderson23\\Desktop\\classes\\877-stats\\hw3\\A549_COVID_NormExp.csv", row.names=1)
common_genes <- intersect(rownames(a549), rownames(nhbe))
length(common_genes)
a549 <- a549[common_genes,]
nhbe <- nhbe[common_genes,]
combined <- cbind(a549, nhbe)
combined <- as.matrix(combined)
str(combined)
head(combined, 1)

conditions_4group <- as.factor(c('C1', 'C1', 'C1', 
                                 'C2', 'C2', 'C2',
                                 'C3', 'C3', 'C3',
                                 'C4', 'C4', 'C4'))
patterns <- GetPatterns(conditions_4group)
Sizes_4group <- MedianNorm(combined)
Sizes_4group

EBOut <- EBMultiTest(Data=combined, 
                     Conditions = conditions_4group, 
                     sizeFactors = Sizes_4group,
                     maxround = 1)
EBOut

multipp <- GetMultiPP(EBOut)
names(multipp)
multipp$PP
pps <- multipp$PP

df <- as.data.frame(apply(pps[, -1], 1, sum))
length(df[df[1] > 0.95])
length(df[df[1] > 0.90])











