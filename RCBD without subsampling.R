# RCBD: NO subsampling 
# with multiple comparison

data=read.csv("FileName.csv", header=T, sep=",") ### 
# Noticed that the data structure should be like this: 
# 'class1' 'class1' 'value'
#    1        1      0.525
#    1        2      0.256
#    2        2      0.695

unbalanced <- 0 ###
class1 <- factor(data$class1) ###
class2 <- factor(data$class2) ###
value <- data$value ###
class1_n_level <- nlevels(class1)
class2_n_level <- nlevels(class2)
level_of_class1 <- levels(factor(class1))
level_of_class2 <- levels(factor(class2))
total_n <- length(value)

# ============
anova_result<-aov(value ~ factor(class1)+factor(class2), data=data)
summary(anova_result)

# Unbalanced group: 
# if (unblanced==1){
  #library(car)
  #Anova(anova_result, type = "III")
#}

#qf(0.95,1,16)

# ========================
# Multiple comparisons
# ========================

# Tukey's HSD
# Count by TukeyHSD()
# TukeyHSD(x, which, ordered = FALSE, conf.level = 0.95, ...) 
#          x: aov fit
#          which: which factor?

# Bonferroni method (pairwise comparison)
# pairwise.t.test(Y, Treatment, p.adj = "bonferroni")
# The result represents p-value of each comparison. If < 0.05, then there's significant difference.


if (class1_n_level>2){
  TukeyHSD(anova_result,'which', ordered = FALSE, conf.level = 0.95) ###
}else{
  pairwise.t.test(data$value, data$class1, p.adj = "bonferroni")
}
