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
# 
# out_model <- function(x, depleted){ # x is the row of data (i.e., data of each protein spot)
#   if (depleted == "TRUE"){
#     cy <- as.factor(rep(c(5,5,3,3),4))
#   } else {
#     cy <- as.factor(rep(c(3,3,5,5),4))
#   }
#   animal <- as.factor(rep(1:8, each = 2))
#   line <- as.factor(rep(c(1,2), each = 8))
#   # check if all obsetvations for one Cy and all observation for line  are missing or not
#   
#   if ((sum(is.na(x[cy==3]))==8|sum(is.na(x[cy==5]))==8) & 
#         (sum(is.na(x[line==1]))==8|sum(is.na(x[line==2]))==8)){
#     model <- lmer(x~ (1|animal), na.action="na.omit")
#     mean_est <- summary(model)$coeff[,1]
#     sd_est <- as.vector(sqrt(summary(model)$vcov))
#     #str(summary(model))
#   } 
#   
#   # check if all observation for one cy missing and for line not missing
#   if ((sum(is.na(x[cy==3]))==8|sum(is.na(x[cy==5]))==8) & 
#         ((sum(is.na(x[line==1]))!=8)&(sum(is.na(x[line==2]))!=8))){
#     model <- lmer(x~ line + (1|animal), na.action="na.omit")
#     mean_est <- summary(model)$coeff[1,1] + summary(model)$coeff[2,1]/2
#     sd_est <- as.vector(sqrt(t(c(1,1/2)) %*%summary(model)$vcov%*%c(1,1/2)))
#   }
#   
#   # check if all observation for Cy not missing, Line missing
#   if ((sum(is.na(x[cy==3]))!=8&sum(is.na(x[cy==5]))!=8) & 
#         ((sum(is.na(x[line==1]))==8)|(sum(is.na(x[line==2]))==8))){
#     model <- lmer(x~ cy + (1|animal), na.action="na.omit")
#     mean_est <- summary(model)$coeff[1,1] + summary(model)$coeff[2,1]/2
#     sd_est <- as.vector(sqrt(t(c(1,1/2)) %*%summary(model)$vcov%*%c(1,1/2)))
#   }
#   
#   if ((sum(is.na(x[cy==3]))!=8& sum(is.na(x[cy==5]))!=8) & 
#         ((sum(is.na(x[line==1]))!=8)&(sum(is.na(x[line==2]))!=8))){ 
#     model <- lmer(x~ cy + line + (1|animal), na.action="na.omit")
#     #str(summary(model))
#     mean_est <- summary(model)$coeff[1,1] + 
#       summary(model)$coeff[2,1]/2 + 
#       summary(model)$coeff[3,1]/2
#     
#     sd_est <- as.vector(sqrt(t(c(1, 1/2, 1/2))%*% summary(model)$vcov %*%c(1,1/2, 1/2)))  }
#   
#   
#   return(c(mean_est, sd_est))
# }
# 
# #sum(whole_nocy_ind)
# depleted_out <- laply(1:dim(dat_final)[1], function(i)out_model(deplete[i,], depleted = "TRUE"))
# 
# whole_out <- laply(1:dim(dat_final)[1], function(i)out_model(whole[i,], depleted = "FALSE"))
# 
# colnames(depleted_out) <- colnames(whole_out) <- c("lsmean","logsd")
# sd_out <- data.frame(
#   logsd = log(c(whole_out[, "logsd"], depleted_out[,"logsd"])),
#   lsmean = c(whole_out[, "lsmean"], depleted_out[,"lsmean"]),
#   sample = rep(c("whole", "depleted"), each = dim(whole_out)[1]))
# 
# ## logsd plot ####
# 
# p <- ggplot(sd_out, aes(sample, logsd))
# p + geom_boxplot(aes(fill = sample)) + 
#   ggtitle("Standard Error for Each Sample Type") + 
#   scale_fill_discrete(name= "Sample Type") + 
#   xlab("Sample Type") + 
#   ylab("log(Standard Error)") + 
#   theme(text = element_text(size=18))
# 
# 
# wilcox.test(log(depleted_out[,2]), log(whole_out[, 2]), alternative = "greater", paired = TRUE)
# #wilcox.test(depleted_out[,2], whole_out[, 2], alternative = "two.sided", paired = TRUE)
# 
# 
# #wilcox.test(log(depleted_out[,2]/whole_out[, 2]), alternative = "two.sided")
# wilcox.test(log(depleted_out[,2]/whole_out[, 2]), alternative = "greater")
# 
# 
# 
# ## lsmean plot ####
# 
# p <- ggplot(sd_out, aes(sample, lsmean))
# 
# # p <- ggplot(subset(sd_out, , aes(sample, lsmean))
# 
# p + geom_boxplot(aes(fill = sample)) + 
#   ggtitle("LSmean for Each Sample Type") + 
#   scale_fill_discrete(name= "Sample Type") + 
#   xlab("Sample Type") + 
#   ylab("LSmean") + 
#   theme(text = element_text(size=18)) + 
#   ylim(-10, 10)
# 
# summary(depleted_out[,1])
# summary(whole_out[,1])
# wilcox.test((depleted_out[which(depleted_out[,1] >=-10),1]), (whole_out[which(depleted_out[,1] >=-10), 1]), 
#             alternative = "greater", paired = TRUE)
# 
# wilcox.test((depleted_out[which(depleted_out[,1] >=-10),1])- (whole_out[which(depleted_out[,1] >=-10), 1]), 
#             alternative = "greater")
# #wilcox.test(depleted_out[,2], whole_out[, 2], alternative = "two.sided", paired = TRUE)
# 
# 
# #wilcox.test(log(depleted_out[,2]/whole_out[, 2]), alternative = "two.sided")
# wilcox.test(log(depleted_out[,2]/whole_out[, 2]), alternative = "greater")
# 
# 
# ### Plot of logsd vs abundance ###
# qplot(lsmean, logsd, data=sd_out, colour=sample, 
#       xlim =c(-10, 10))
# 
# p1 <- ggplot(sd_out, aes(lsmean, logsd), colour = sample) + 
#   geom_point(aes(colour = sample)) + 
#   xlim(-10, 10)+ 
#   theme(text = element_text(size=18))+ 
#   ggtitle("Logsd vs. abundance of protein spots") +
#   geom_smooth()
# p1
# 
# 
# ##
# 
# deplete_nocy_ind <- laply(1:dim(dat_final)[1], function(i)(sum(is.na(deplete[i,][cy==3]))==8|sum(is.na(deplete[i,][cy==5]))==8) & 
#                             (sum(is.na(deplete[i,][line==1]))==8|sum(is.na(deplete[i,][line==2]))==8))
# #sum(deplete_nocy_ind)
# 
# whole_nocy_ind <- laply(1:dim(dat_final)[1], function(i)(sum(is.na(whole[i,][cy==3]))==8|sum(is.na(whole[i,][cy==5]))==8) & 
#                           (sum(is.na(whole[i,][line==1]))==8|sum(is.na(whole[i,][line==2]))==8))
# 
# 
# str(deplete)
# group[group=="Whole"] 
# head(deplete)
# 
# cy <-  rep(c(5,5,3,3),4)
# deplete_nocy_ind <- laply(1:dim(dat_final)[1], function(i)
#   (sum(is.na(deplete[i,][cy==3]))==8|sum(is.na(deplete[i,][cy==5]))==8))
# #sum(deplete_nocy_ind)
# which(deplete_nocy_ind)
# cy <-  rep(c(3,3,5,5),4)
# whole_nocy_ind <- laply(1:dim(dat_final)[1], function(i)
#   (sum(is.na(whole[i,][cy==3]))==8|sum(is.na(whole[i,][cy==5]))==8))
# #sum(whole_nocy_ind)
# which(whole_nocy_ind)
# whole[24,]
# deplete[24,]
# 
# 
# 
# cy <-  rep(c(5,5,3,3),4)
# deplete_noline_ind <- laply(1:dim(dat_final)[1], function(i)
#   (sum(is.na(deplete[i,][cy==3]))!=8& sum(is.na(deplete[i,][cy==5]))!=8) & 
#     ((sum(is.na(deplete[i,][line==1]))!=8)&(sum(is.na(deplete[i,][line==2]))!=8)))
# #sum(deplete_noline_ind)
# which(deplete_noline_ind)
# cy <-  rep(c(3,3,5,5),4)
# whole_nocy_ind <- laply(1:dim(dat_final)[1], function(i)
#   (sum(is.na(whole[i,][line ==1]))==8|sum(is.na(whole[i,][line ==2]))==8))
# #sum(whole_nocy_ind)
# which(whole_nocy_ind)
# whole[24,]
# deplete[24,]
# 
# ##Low RFI = 1209, 4810, 1906, 1807
# ##High RFI = 3106, 3908, 2107, 2712
# 
# ## linear mixed effect model 
# model <- lmer(deplete[1,]~ cy + line + (1|animal), na.action="na.omit")
# 
# str(summary(model))
####Model in Rnw. file####
out_model <- function(x, depleted){ # x is the row of data (i.e., data of each protein spot)
  if (depleted == "TRUE"){
    cy <- as.factor(rep(c(5,5,3,3),4))
  } else {
    cy <- as.factor(rep(c(3,3,5,5),4))
  }
  animal <- as.factor(rep(1:8, each = 2))
  line <- as.factor(rep(c(1,2), each = 8))
  # check if all obsetvations for one Cy are missing or not
  
  if ((sum(is.na(x[cy==3]))==8|sum(is.na(x[cy==5]))==8) & # if all cy is missing
        (sum(is.na(x[line==1]))==8|sum(is.na(x[line==2]))==8)){# if all Line is missing
    model <- lmer(x~ (1|animal), na.action="na.omit")
    s_model <- summary(model)
    mean_est <- s_model$coeff[,1]
    se_est <- as.vector(sqrt(s_model$vcov))
    #str(s_model)
  } 
  if ((sum(is.na(x[cy==3]))==8|sum(is.na(x[cy==5]))==8) & # if cy is missing
        ((sum(is.na(x[line==1]))!=8)&(sum(is.na(x[line==2]))!=8))){ # if line is not missing
    model <- lmer(x~ line + (1|animal), na.action="na.omit")
    s_model <- summary(model)
    mean_est <- s_model$coeff[1,1] + s_model$coeff[2,1]/2
    se_est <- as.vector(sqrt(t(c(1,1/2)) %*%s_model$vcov%*%c(1,1/2)))
  }
  
  if ((sum(is.na(x[cy==3]))!=8&sum(is.na(x[cy==5]))!=8) &  # if cy is not missing
        ((sum(is.na(x[line==1]))==8)|(sum(is.na(x[line==2]))==8))){ # if line is missing
    model <- lmer(x~ cy + (1|animal), na.action="na.omit")
    s_model <- summary(model)
    mean_est <- s_model$coeff[1,1] + s_model$coeff[2,1]/2
    se_est <- as.vector(sqrt(t(c(1,1/2)) %*%s_model$vcov%*%c(1,1/2)))
  }
  
  if ((sum(is.na(x[cy==3]))!=8& sum(is.na(x[cy==5]))!=8) & 
        ((sum(is.na(x[line==1]))!=8)&(sum(is.na(x[line==2]))!=8))){ 
    model <- lmer(x~ cy + line + (1|animal), na.action="na.omit")
    s_model <- summary(model)
    #str(s_model)
    mean_est <- s_model$coeff[1,1] + 
      s_model$coeff[2,1]/2 + 
      s_model$coeff[3,1]/2
    
    se_est <- as.vector(sqrt(t(c(1, 1/2, 1/2))%*% s_model$vcov %*%c(1,1/2, 1/2)))  }
  
  
  return(c(se_est, mean_est))
}
se_depleted <- laply(1:dim(dat_final)[1], function(i)out_model(deplete[i,], depleted = "TRUE")[1])

se_whole <- laply(1:dim(dat_final)[1], function(i)out_model(whole[i,], depleted = "FALSE")[1])

lsmean_depleted <- laply(1:dim(dat_final)[1], function(i)out_model(deplete[i,], depleted = "TRUE")[2])

lsmean_whole <- laply(1:dim(dat_final)[1], function(i)out_model(whole[i,], depleted = "FALSE")[2])


log_se <- data.frame(
  logse = log(c(se_whole, se_depleted)), 
  sample = rep(c("whole", "depleted"), each = length(se_whole)))
# write.table(log_se, file = "log_se.txt")

p <- ggplot(log_se, aes(sample, logse))
p + geom_boxplot(aes(fill = sample)) + 
  ggtitle("Standard Error for Each Sample Type") + 
  scale_fill_discrete(name= "Sample Type") + 
  xlab("Sample Type") + 
  ylab("log(Standard Error)") + 
  theme(text = element_text(size=11))



par(cex=.8)
plot(log(se_whole), log(se_depleted) ,
     xlim = c(-3, 4), ylim = c(-3, 4))
abline(a =0, b = 1)

difference_logse <- log(se_depleted) - log(se_whole)
average_logse <- (log(se_depleted) + log(se_whole))/2

par(cex=.8)
plot(average_logse, difference_logse,
     xlim = c(-3, 4), ylim = c(-3, 4))
abline(h = 0)
