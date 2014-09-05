library("lme4")
library("plyr")
library("ggplot2")
library("MASS")
# read data

dat <- read.csv("239 tech paper_08_19_2014.csv")

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
#dat_used[3,]

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
mean_deplete <- apply(deplete, 1, mean, na.rm =TRUE)

# ?na.action
# ?apply
whole <- dat_final[,-1][,group=="Whole"]
dim(whole)

mean_whole <- apply(whole, 1, mean, na.rm = TRUE)
group[group == "Whole"]
head(whole[,c(3,4)])
head(deplete[,c(3,4)], 20)
dim(deplete)
length(which(deplete[,3] == deplete[,4]))/dim(deplete)[1]
length(which(whole[,3] == whole[,4]))/dim(whole)[1]

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

depleted_out <- laply(1:dim(dat_final)[1], function(i)out_model(deplete[i,], depleted = "TRUE"))

whole_out <- laply(1:dim(dat_final)[1], function(i)out_model(whole[i,], depleted = "FALSE"))
#sum(depleted_out[,1] < whole_out[,1])

colnames(depleted_out) <- colnames(whole_out) <- c("mean","sd")
sd_out <- data.frame(
  logsd = log(c(whole_out[, "sd"], depleted_out[,"sd"])), 
  sample = rep(c("whole", "depleted"), each = dim(whole_out)[1]))
head(sd_out)
write.table(sd_out, file = "log_sd.txt" ,
            sep = "\t", 
            row.names = FALSE)
p <- ggplot(sd_out, aes(sample, logsd))
p + geom_boxplot(aes(fill = sample)) + 
  ggtitle("Standard Error for Each Sample Type") + 
  scale_fill_discrete(name= "Sample Type") + 
  xlab("Sample Type") + 
  ylab("log(Standard Error)") + 
  theme(text = element_text(size=18))
#dev.off()
bp <- ggplot(data=sd_out, aes(x=sample, y=logsd, fill=sample)) + 
  geom_boxplot()
bp + guides(fill=FALSE)
head(sd_out)
wilcox.test(depleted_out[,2], whole_out[, 2], alternative = "greater", paired = TRUE)
wilcox.test(depleted_out[,2], whole_out[, 2], alternative = "two.sided", paired = TRUE)


wilcox.test(log(depleted_out[,2]/whole_out[, 2]), alternative = "two.sided")
wilcox.test(log(depleted_out[,2]/whole_out[, 2]), alternative = "greater")

## 
logsd = log(depleted_out[,"sd"]/whole_out[, "sd"])
boxplot(logsd)
abline(h = 0)

mixed_out <- data.frame(
  logsd = log(c(whole_out[, "sd"], depleted_out[,"sd"])), 
  lsmean = c(whole_out[, "mean"], depleted_out[,"mean"]),
  a_mean = c(mean_whole, mean_deplete),
  sample = rep(c("whole", "depleted"), each = dim(whole_out)[1]))
head(mixed_out)

p <- ggplot(mixed_out, aes(lsmean, logsd)) + geom_point()
# With one variable
p + facet_grid(. ~ sample)

p <- ggplot(mixed_out, aes(a_mean, logsd)) + geom_point()
# With one variable
p + facet_grid(. ~ sample)
cor(mixed_out$a_mean, mixed_out$lsmean)
#dim(mixed_out)[1]/2

