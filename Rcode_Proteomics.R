library("lme4")
library("plyr")
library("ggplot2")
library("MASS")
# read data

#dat <- read.csv("P:/RFI Protemics/twofileswithproteindata/Proteins of Interst for Tech Paper.csv")

#dat <- read.csv("/run/user/1000/gvfs/smb-share:server=smb.stat.iastate.edu,share=ntyet/RFI Protemics/twofileswithproteindata/Proteins of Interst for Tech Paper.csv")
dat <- read.csv("Proteins of Interst for Tech Paper.csv")

dim(dat)

x <- colnames(dat) # column name of data
x


y <- c( "Image","X", paste("X.", 1:31, sep = "")) # Column name unsued

# Column name used in analysis including the id of protein and 16 gels  
# of 8 animals, each gel runs twice

use_col <- setdiff(x,y)
use_col
dat_used <- as.matrix(dat[, use_col])
dim(dat_used)


# The first row is group classification, the second row is the name of Std. Abund, therefore the data actually in use is the dat_used except the first 2 rows
dat_final <- matrix(as.numeric(dat_used[-c(1:2),]), 
                    nrow = nrow(dat_used)-2,
                    ncol = ncol(dat_used),
                    byrow = F) 
dim(dat_final)

head(dat_final)

# Group of each sample

group <- dat_used[1,-1] 
group

# Obtain data for each group: Depelted and Not.Depleted (Whole), the first column of dat_final contains name of protein
deplete <- dat_final[,-1][,group=="Depleted"]
dim(deplete)


whole <- dat_final[,-1][,group=="Whole"]
dim(whole)


# Find out which Cy is used for each run

group[group=="Depleted"] # Cy for Depleted group: rep(c(5,5,3,3),4)

group[group=="Whole"] # Cy for Whole group: rep(c(3,3,5,5),4)

#metadata <- read.table("/run/user/1000/gvfs/smb-share:server=smb.stat.iastate.edu,share=ntyet/RFI Protemics/twofileswithproteindata/Meta.data2.txt")
#metadata <- read.table("Meta.data2.txt")
#str(metadata)
# Sample.Name <- c("1807", "4810", "1209", "1906", "3908", "3106", "2107", "2712")
# line <- rep(c(1,2), each = 8)
# x <- rnorm(16, 0, 1)
# depleted <- "TRUE"

out_model <- function(x, depleted){ # x is the row of data (i.e., data of each protein spot)
  if (depleted == "TRUE"){
    cy <- as.factor(rep(c(5,5,3,3),4))
  } else {
    cy <- as.factor(rep(c(3,3,5,5),4))
  }
  animal <- as.factor(rep(1:8, each = 2))
  line <- as.factor(rep(c(1,2), each = 8))
  # check if all obsetvations for one Cy are missing or not
  
  if ((sum(is.na(x[cy==3]))==8|sum(is.na(x[cy==5]))==8) & 
        (sum(is.na(x[line==1]))==8|sum(is.na(x[line==2]))==8)){
    model <- lmer(x~ (1|animal), na.action="na.omit")
    mean_est <- summary(model)$coeff[,1]
    sd_est <- as.vector(sqrt(summary(model)$vcov))
    #str(summary(model))
  } 
  if ((sum(is.na(x[cy==3]))==8|sum(is.na(x[cy==5]))==8) & 
        ((sum(is.na(x[line==1]))!=8)&(sum(is.na(x[line==2]))!=8))){
    model <- lmer(x~ line + (1|animal), na.action="na.omit")
    mean_est <- summary(model)$coeff[1,1] + summary(model)$coeff[2,1]/2
    sd_est <- as.vector(sqrt(t(c(1,1/2)) %*%summary(model)$vcov%*%c(1,1/2)))
  }
  
  if ((sum(is.na(x[cy==3]))!=8|sum(is.na(x[cy==5]))!=8) & 
        ((sum(is.na(x[line==1]))==8)&(sum(is.na(x[line==2]))==8))){
    model <- lmer(x~ cy + (1|animal), na.action="na.omit")
    mean_est <- summary(model)$coeff[1,1] + summary(model)$coeff[2,1]/2
    sd_est <- as.vector(sqrt(t(c(1,1/2)) %*%summary(model)$vcov%*%c(1,1/2)))
  }
  
  if ((sum(is.na(x[cy==3]))!=8& sum(is.na(x[cy==5]))!=8) & 
        ((sum(is.na(x[line==1]))!=8)&(sum(is.na(x[line==2]))!=8))){ 
    model <- lmer(x~ cy + line + (1|animal), na.action="na.omit")
    #str(summary(model))
    mean_est <- summary(model)$coeff[1,1] + 
      summary(model)$coeff[2,1]/2 + 
      summary(model)$coeff[3,1]/2
    
    sd_est <- as.vector(sqrt(t(c(1, 1/2, 1/2))%*% summary(model)$vcov %*%c(1,1/2, 1/2)))  }
  
  
  return(c(mean_est, sd_est))
}

depleted_out <- laply(1:236, function(i)out_model(deplete[i,], depleted = "TRUE"))

whole_out <- laply(1:236, function(i)out_model(whole[i,], depleted = "FALSE"))

colnames(depleted_out) <- colnames(whole_out) <- c("mean","sd")
sd_out <- data.frame(
  logsd = log(c(whole_out[, "sd"], depleted_out[,"sd"])), 
  sample = rep(c("whole", "depleted"), each = dim(whole_out)[1]))

p <- ggplot(sd_out, aes(sample, logsd))
p + geom_boxplot(aes(fill = sample)) + 
  ggtitle("Standard Error for Each Sample Type") + 
  scale_fill_discrete(name= "Sample Type") + 
  xlab("Sample Type") + 
  ylab("log(Standard Error)") + 
  theme(text = element_text(size=18))


wilcox.test(log(depleted_out[,2]), log(whole_out[, 2]), alternative = "greater", paired = TRUE)
wilcox.test(depleted_out[,2], whole_out[, 2], alternative = "two.sided", paired = TRUE)


wilcox.test(log(depleted_out[,2]/whole_out[, 2]), alternative = "two.sided")
wilcox.test(log(depleted_out[,2]/whole_out[, 2]), alternative = "greater")


