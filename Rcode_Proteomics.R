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

use.col <- setdiff(x,y)
use.col
dat.used <- as.matrix(dat[, use.col])
dim(dat.used)


# The first row is group classification, the second row is the name of Std. Abund, therefore the data actually in use is the dat.used except the first 2 rows
dat.final <- matrix(as.numeric(dat.used[-c(1:2),]), 
                    nrow = nrow(dat.used)-2,
                    ncol = ncol(dat.used),
                    byrow = F) 
dim(dat.final)

head(dat.final)

# Group of each sample

group <- dat.used[1,-1] 
group

# Obtain data for each group: Depelted and Not.Depleted (Whole), the first column of dat.final contains name of protein
depleted.dat <- dat.final[,-1][,group=="Depleted"]
dim(depleted.dat)
depleted.dat[1,]

not.depleted.dat <- dat.final[,-1][,group=="Whole"]
dim(not.depleted.dat)
not.depleted.dat[1,]

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

result <- function(x, depleted){ # x is the row of data (i.e., data of each protein spot)
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
    mean.estimate <- summary(model)$coeff[,1]
    sd.estimate <- as.vector(sqrt(summary(model)$vcov))
    #str(summary(model))
  } 
  if ((sum(is.na(x[cy==3]))==8|sum(is.na(x[cy==5]))==8) & 
        ((sum(is.na(x[line==1]))!=8)&(sum(is.na(x[line==2]))!=8))){
    model <- lmer(x~ line + (1|animal), na.action="na.omit")
    mean.estimate <- summary(model)$coeff[1,1] + summary(model)$coeff[2,1]/2
    sd.estimate <- as.vector(sqrt(t(c(1,1/2)) %*%summary(model)$vcov%*%c(1,1/2)))
  }
  
  if ((sum(is.na(x[cy==3]))!=8|sum(is.na(x[cy==5]))!=8) & 
        ((sum(is.na(x[line==1]))==8)&(sum(is.na(x[line==2]))==8))){
    model <- lmer(x~ cy + (1|animal), na.action="na.omit")
    mean.estimate <- summary(model)$coeff[1,1] + summary(model)$coeff[2,1]/2
    sd.estimate <- as.vector(sqrt(t(c(1,1/2)) %*%summary(model)$vcov%*%c(1,1/2)))
  }
  
  if ((sum(is.na(x[cy==3]))!=8& sum(is.na(x[cy==5]))!=8) & 
        ((sum(is.na(x[line==1]))!=8)&(sum(is.na(x[line==2]))!=8))){ 
    model <- lmer(x~ cy + line + (1|animal), na.action="na.omit")
    #str(summary(model))
    mean.estimate <- summary(model)$coeff[1,1] + 
      summary(model)$coeff[2,1]/2 + 
      summary(model)$coeff[3,1]/2
    
    sd.estimate <- as.vector(sqrt(t(c(1, 1/2, 1/2))%*% summary(model)$vcov %*%c(1,1/2, 1/2)))  }
  
  
  return(c(mean.estimate, sd.estimate))
}

dim(depleted.dat)
dim(not.depleted.dat)
depleted.out <- laply(1:236, function(i)result(depleted.dat[i,], depleted = "TRUE"))

not.depleted.out <- laply(1:236, function(i)result(not.depleted.dat[i,], depleted = "FALSE"))

colnames(depleted.out) <- colnames(not.depleted.out) <- c("mean","sd")

summary(depleted.out) # depleted group
summary(not.depleted.out) # not.depleted group 
# mean_group <- as.data.frame(cbind(value =c(depleted.out[,"mean"], not.depleted.out[,"mean"]), group = rep(c("depleted", "whole"), each = dim(depleted.out)[1])))
# str(mean_group)
# qplot(mean_group)
# qplot(group, value, data = mean_group, geom= "boxplot")
#boxplot(depleted.out)


mean(depleted.out[,2] > not.depleted.out[,2])


wilcox.test(depleted.out[,2], not.depleted.out[, 2], alternative = "greater", paired = TRUE)
wilcox.test(depleted.out[,2], not.depleted.out[, 2], alternative = "two.sided", paired = TRUE)
?wilcox.test
wilcox.test(log(depleted.out[,2]/not.depleted.out[, 2]), alternative = "two.sided")
wilcox.test(log(depleted.out[,2]/not.depleted.out[, 2]), alternative = "greater")

# check for outlier
sum(abs(depleted.out[,1])>10)

# Delete outlier
depleted.out2 <- depleted.out[-which(abs(depleted.out[,1])>10),]

not.depleted.out2 <- not.depleted.out[-which(abs(depleted.out[,1])>10),]

summary(depleted.out2) # depleted group
summary(not.depleted.out2) # not.depleted group 


mean(depleted.out2[,2] > not.depleted.out2[,2])


wilcox.test(depleted.out2[,2], not.depleted.out2[, 2], alternative = "greater", paired = TRUE)
wilcox.test(depleted.out2[,2], not.depleted.out2[, 2], alternative = "two.sided", paired = TRUE)
wilcox.test(log(depleted.out2[,2]/not.depleted.out2[, 2]), alternative = "two.sided")
wilcox.test(log(depleted.out2[,2]/not.depleted.out2[, 2]), alternative = "greater")

depleted.out2 <- as.data.frame(depleted.out2)
not.depleted.out2 <- as.data.frame(not.depleted.out2)
depleted.out2$group <- "depleted"
not.depleted.out2$group <- "whole"
out2 <- rbind(depleted.out2, not.depleted.out2)

str(out2)
boxplot(out2$sd ~ out2$group)
boxplot(out$sd ~ out$group)
qplot(y = depleted.out2[,1], x = not.depleted.out2[,1], xlim = c(-8, 8), 
      ylim = c(-8, 8), main = "mean of each group(deleted outlier)",
      ylab = "depleted", 
      xlab = "whole")
qplot(y = depleted.out2[,2], x = not.depleted.out2[,2], xlim = c(-8, 8), 
      ylim = c(-8, 8), main = "sd of each group(deleted outlier)",
      ylab = "depleted", 
      xlab = "whole")


qplot(x = group, y = mean, 
        data = out2, geom = "boxplot", colour =group, 
        main = "mean for each group")

qplot(x = group, y = sd, 
      data = out2, geom = "boxplot", colour =group,
      main = "sd for each group")

  

depleted.out <- as.data.frame(depleted.out)
not.depleted.out <- as.data.frame(not.depleted.out)
depleted.out$group <- "depleted"
not.depleted.out$group <- "whole"
out <- rbind(depleted.out, not.depleted.out)


qplot(y = depleted.out[,1], x = not.depleted.out[,1], xlim = c(-8, 8), 
      ylim = c(-8, 8), main = "mean of each group(deleted outlier)",
      ylab = "depleted", 
      xlab = "whole")
qplot(y = depleted.out[,2], x = not.depleted.out[,2], xlim = c(-8, 8), 
      ylim = c(-8, 8), main = "sd of each group(deleted outlier)",
      ylab = "depleted", 
      xlab = "whole")


qplot(x = group, y = mean, 
      data = out, geom = "boxplot", colour =group, 
      main = "mean for each group")

qplot(x = group, y = sd, 
      data = out, geom = "boxplot", colour =group,
      main = "sd for each group")

