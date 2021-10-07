# Randomized complete block design with subsampling

##=================
# model: 
# Y_ijl = MU + ALPHA_i + BETA_j + EPSHLON_ij + DELTA_ijl
# treatment: i = 1,....,k 
# block: j = 1,....,b
# subsample within each plot: l = 1,.....,s
# experiment error: epshlon~N(0,sigma^2_epshlon) 
# sampling error: delta~N(0,sigma^2_delta)
##=================


# Noticed that the data structure should be like this: 
# 'block' 'class1'  'subsample'  'value'
#    1        1         1         0.525
#    1        1         2         0.256
#    2        2         1         0.695
#    2        2         2         0.153

data=read.csv("C:/Users/User/Downloads/school/lecture/Statistics/RCBD_withsubsample_hw.csv", header=T, sep=",") ### 

attach(data)

alpha <- 0.05 ###
block <- data$block ###
class1 <- data$class1 ###
subsample<-data$subsample ###
value <- data$value ###


block<-factor(block)
class1<-factor(class1)
subsample<-factor(subsample)

nlevel_of_block<-nlevels(block)
nlevel_of_class1<-nlevels(class1)
nlevel_of_subsample<-nlevels(subsample)
total_n <- length(value)

level_block <- levels(block)
level_class1 <- levels(class1)
level_subsample <- levels(subsample)

df_block <- nlevel_of_block-1
df_class1 <- nlevel_of_class1-1
df_exp.error <- (nlevel_of_block-1)*(nlevel_of_class1-1)
df_subsample.error <- nlevel_of_block*nlevel_of_class1*(nlevel_of_subsample-1)

critical_value_block <- qf(1-alpha,df1 = df_block,df2 = df_exp.error)
critical_value_class1 <- qf(1-alpha,df1 = df_class1,df2 = df_exp.error)
critical_value_exp.error <- qf(1-alpha,df1 = df_exp.error,df2 =df_subsample.error )

### level_class <- levels(class)


#====================================
# ANOVA test: count by my own
#====================================
# 1. check for the data structure
if (is.data.frame(data)==TRUE){
  print("Data is loaded as data.frame. Go on.")
}else{
  print("Data is not loaded as data.frame. Please check the raw data and load it again.")
}

# 2. Count MSTotal 
mean_total <- mean(value)
SSTotal <- sum((data$value-mean_total)^2)
MSTotal <- SSTotal/(total_n-1)

# 3. Count MSB,MSTR,MSexp,MSse

level_name <- c()
mean_block <- c()
mean_treat <- c()
mean_cell <- matrix(nrow = nlevel_of_class1, ncol = nlevel_of_block)
sub_SSexp <- 0



   # 3.1 MSB
for (b in 1:nlevel_of_block){
  sub_data <- subset(data, block==level_block[b])
  mean_block[b] <- mean(sub_data$value)
}

   # 3.2 MSTR
for (t in 1:nlevel_of_class1){
  sub_data <- subset(data, class1==level_class1[t])
  mean_treat[t] <- mean(sub_data$value)
}

   # 3.3 MSExp.error
for (t in 1:nlevel_of_class1){
  for (b in 1:nlevel_of_block) {
    sub_data <- subset(data, class1==level_class1[t]&block==level_block[b])
    mean_cell[t,b]<-mean(sub_data$value)
    sub_SSexp <- (mean_cell[t,b] - mean_block[b] - mean_treat[t] + mean_total)^2+sub_SSexp
  }
}


   # 3.4 MSsubsample
value_in_cell <- c()
SSsubsample <- 0 

for (t in 1:nlevel_of_class1){
  for (b in 1:nlevel_of_block) {
      sub_data <- subset(data, class1==level_class1[t]&block==level_block[b])
      value_in_cell <- sub_data$value
      SSsubsample <- sum((value_in_cell-mean_cell[t,b])^2)+SSsubsample
  }
}


SSB <- sum((mean_block-mean_total)^2)*nlevel_of_class1*nlevel_of_subsample
MSB <- SSB/df_block

SSTR <- sum((mean_treat-mean_total)^2)*nlevel_of_block*nlevel_of_subsample
MSTR <- SSTR/df_class1

SSexp <- sub_SSexp*nlevel_of_subsample
MSexp <- SSexp/df_exp.error

MSsubsample <- SSsubsample/df_subsample.error

f_value_block <- MSB/MSexp
f_value_treatment <- MSTR/MSexp
f_value_exp.error <- MSexp/MSsubsample


# 4. Check: SSTotal = SSTR + SSB + SSexp + SSsubsample
cat(" SSTotal =",SSTotal,"\n","SSTR =",SSTR,"\n","SSB =",SSB,"\n","SSexp =",SSexp,"\n","SSsubsample =",SSsubsample,"\n","Check:SSTR + SSB + SSexp + SSsubsample =",(SSTR+SSB+SSexp+SSsubsample),"\n")


# 5. Show: MSTotal, MSTR, MSB, MSexp, MSsubsample
cat(" MSTotal =",MSTotal,"\n","MSTR =",MSTR,"\n","MSB =",MSB,"\n","MSexp =",MSexp,"\n","MSsubsample =",MSsubsample,"\n","F value (block)=",f_value_block,"\n","F value (Treatment)=",f_value_treatment,"\n","F value (exp.error)=",f_value_exp.error,"\n","critical value (block) =", critical_value_block,"  df1-> ",df_block,"  df2->",df_exp.error,"\n","critical value (Treatment) =", critical_value_class1,"  df1-> ",df_class1,"  df2->",df_exp.error,"\n","critical value (exp.error) =", critical_value_exp.error,"  df1-> ",df_exp.error,"  df2->",df_subsample.error,"\n")


#-----

#=============
# post hoc
#=============

# pairwise test with bonferroni
# pairwise.t.test(data$value, data$class1, p.adj = "bonferroni", pool.sd = FALSE)

# Bonferroni counted by my own
print(level_class1)
i=0
j=0

i=1 ###
j=4 ###
alpha_for_MC <- 0.05 ### if you want 95% CI
k <- nlevel_of_class1 ### check: which treatment you want to compare
n1 <- length(data$value[data$class==level_class1[i]]) ###
n2 <- length(data$value[data$class==level_class1[j]]) ###
cirital_t_bonferroni <- qt(alpha_for_MC/2*(choose(k,2)),df=df_exp.error)
Bonferroni <- cirital_t_bonferroni*sqrt(MSexp*(1/n1+1/n2))
cat(" Bonferroni = ",Bonferroni,"\n","mu1:",level_class1[i]," mu2:",level_class1[j],"\n","C.I.= ",(mean_treat[i]-mean_treat[j])-Bonferroni," , ",(mean_treat[i]-mean_treat[j])+Bonferroni,"\n")

#==========================
# Plot interaction 
#==========================

# factor1: the one you want to put in the x-axis
factor1 <- factor(class1) ### 
nlevel_of_factor1 <- length(levels(factor1))
# factor2: the one you want to distinct its levels as two or more saparated lines  
factor2 <- factor(data$block) ###
nlevel_of_factor2 <- length(levels(factor2))

factor1_mean <- c()
level_factor2 <- c()
level_factor1 <- c()
i=0
j=0
x=0
for (i in 1:nlevel_of_class1) {
  for (j in 1:nlevel_of_block) {
    x=x+1
    sub_plot <- subset(data, class1==level_class1[i]&block==level_block[j]) ###
    factor1_mean[x] <- mean(sub_plot$value) ###
    level_factor2[x] <- j
    level_factor1[x] <- i
  }
}
factor1_subset <- data.frame(factor1_mean,level_factor2,level_factor1) # col1: factor1 mean, col2:level of factor2

# Plot!!! 
color <- c("#636B46","#99CED4","#EEB6B7","#6E7376","#056571","#494E68","#98878F","#985E6D","#CDA34F","#E9E7DA") # 9 

plot(factor1_subset$level_factor1,factor1_subset$factor1_mean, xlab = "level_factor1", ylab = "mean value of factor1") ###
for (m in 1:nlevel_of_factor2) {
  lines(factor1_subset$factor1_mean[level_factor2==m],col=color[m]) ###
  points(factor1_subset$factor1_mean[level_factor2==m],col=color[m]) ###
}

