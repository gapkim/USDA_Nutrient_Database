##############################
### Initial Setup

knitr::opts_chunk$set(comment=NA, echo=FALSE, warning=FALSE, message=FALSE,
                      fig.align="center")
options(digits=4)
rm(list=ls())

###
### Data cleaning
###

SR = read.table("ABBREV.txt", header=F, row.names=1, sep="^", quote="~")
SR = na.omit(SR) # remove rows with missing values
SR = SR[row.names(SR) != "13352",] # remove "duplicate" entry
row.names(SR) = SR[,1] # set more meaningful row names
SR = SR[,-1]

names(SR) = c("Water_(g)", "Energ_Kcal", "Protein_(g)", "Lipid_Tot_(g)", "Ash_(g)", "Carbohydrt_(g)", "Fiber_TD_(g)", "Sugar_Tot_(g)", "Calcium_(mg)", "Iron_(mg)", "Magnesium_(mg)", "Phosphorus_(mg)", "Potassium_(mg)", "Sodium_(mg)", "Zinc_(mg)", "Copper_(mg)", "Manganese_(mg)", "Selenium_(µg)", "Vit_C_(mg)", "Thiamin_(mg)", "Riboflavin_(mg)", "Niacin_(mg)", "Panto_Acid_(mg)", "Vit_B6_(mg)", "Folate_Tot_(µg)", "Folic_Acid_(µg)", "Food_Folate_(µg)", "Folate_DFE_(µg)", "Choline_Tot_(mg)", "Vit_B12_(µg)", "Vit_A_IU", "Vit_A_RAE", "Retinol_(µg)", "Alpha_Carot_(µg)", "Beta_Carot_(µg)", "Beta_Crypt_(µg)", "Lycopene_(µg)", "Lut+Zea_(µg)", "Vit_E_(mg)", "Vit_D_µg", "Vit_D_IU", "Vit_K_(µg)", "FA_Sat_(g)", "FA_Mono_(g)", "FA_Poly_(g)", "Cholestrl_(mg)", "GmWt_1", "GmWt_Desc1", "GmWt_2", "GmWt_Desc2", "Refuse_Pct")

SRp = SR[,c(1:46)] # restrict to just the nutrient variables

str(SRp)

#########################
### EDA
#########################

###
### Proximates
###

library(reshape2)
library(ggplot2)

proximates = subset(SRp, select=c("Water_(g)","Protein_(g)","Lipid_Tot_(g)",
                                "Carbohydrt_(g)","Ash_(g)"))

melt_proximates = melt(proximates)

# Figure 1: Boxplots of proximates in all foods

ggplot(data=melt_proximates) +
   geom_boxplot(aes(x=variable, y=value)) +
   facet_wrap(~variable, scale="free", ncol=5) +
   labs(x="Proximates", y="Values")

rownames(SRp)[order(-SRp$`Protein_(g)`)[1:5]]

###
### Minerals
###

minerals = subset(SRp, select=c("Calcium_(mg)", "Iron_(mg)", "Magnesium_(mg)", "Phosphorus_(mg)", "Potassium_(mg)", "Sodium_(mg)", "Zinc_(mg)", "Copper_(mg)", "Manganese_(mg)", "Selenium_(µg)"))

melt_minerals = melt(minerals)

# Figure 2: Bolxplots of minerals in all foods

ggplot(data=melt_minerals) +
   geom_boxplot(aes(x=variable, y=value)) +
   facet_wrap(~variable, scale="free", ncol=5) +
   labs(x="Minerals", y="Values")
   
rownames(SRp)[order(-SRp$`Calcium_(mg)`)[1:5]]

###
### Vitamins and other nutrients
###

vitamin = subset(SRp, select=c("Vit_A_RAE", "Vit_C_(mg)", "Vit_D_IU", "Vit_E_(mg)", "Folate_Tot_(µg)"))

melt_vitamin = melt(vitamin)

# Figure 3: Bolxplots of vitamins and other nutrients in all foods

ggplot(data=melt_vitamin) +
   geom_boxplot(aes(x=variable, y=value)) +
   facet_wrap(~variable, scale="free", ncol=5) +
   labs(x="Vitamins", y="Values")

#########################
### PCA
#########################

pr.out = prcomp(SRp, center=TRUE, scale=TRUE)  # With center=TRUE, scale=TRUE
pr.out$x = -pr.out$x
pr.out$rotation = -pr.out$rotation

pve = pr.out$sdev^2 / length(pr.out$sdev)  # Proportion of variance explained

library(ggplot2)
library(gridExtra)

# Figure 4: Proportion of variance explained and their cumulative sum from PCA

p1 = ggplot() +
   geom_line(aes(x=c(1:length(pr.out$sdev)), y=pve)) +
   geom_point(aes(x=c(1:length(pr.out$sdev)), y=pve), size=3) +
   geom_hline(yintercept=pve[7], color='red') +
   labs(x="Principal component", y="Proportion of variance explained")

p2 = ggplot() +
   geom_line(aes(x=c(1:length(pr.out$sdev)), y=cumsum(pve))) +
   geom_point(aes(x=c(1:length(pr.out$sdev)), y=cumsum(pve)), size=3) +
   geom_hline(yintercept=cumsum(pve)[7], color='red') +
   annotate("text", x=20, y=0.59, label=paste(round(cumsum(pve)[7],3))) +
   labs(x="Principal component", y="Cumulative proportion of variance explained")

grid.arrange(p1, p2, ncol=2)


###
### Heatmap
###

melt_pr.out = melt(pr.out$rotation[,c(1:20)])
colnames(melt_pr.out) = c("Nutrient", "PC", "Value")

# Figure 5: Heatmap of nutrients and first 20 principal components

ggplot(data=melt_pr.out) +
   geom_tile(aes(x=PC, y=Nutrient, fill=Value)) +
   scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0) +
   labs(x="Principal component")
   
library(reshape2)

melt_pr.out1 = melt(pr.out$rotation[,c(1:2)])
colnames(melt_pr.out1) = c("Nutrient", "PC", "Value")

# Figure 6: Full spectrum of the nutrient loadings in PC1 and PC2

ggplot(data=melt_pr.out1) +
   geom_col(aes(x=Nutrient, y=Value)) +
   facet_wrap(~PC, ncol=1) +
   theme(axis.text.x = element_text(angle = 90)) +
   labs(x="Nutrients", y="Loading")

#
# Sanity check
#
pve.mat = matrix(rep(pve, each = 46), nrow=length(pve))

nutrient.impact = apply(pr.out$rotation[,c(1:46)]^2 * pve.mat, 1, sum)
sum(nutrient.impact)  # Should equal 1
###
### Selection of nutrients
###

pve.mat = matrix(rep(pve, each = 7), nrow=length(pve))
nutrient.impact = apply(pr.out$rotation[,c(1:7)]^2 * pve.mat, 1, sum)
melt_NI = melt(nutrient.impact)

# Figure 7: Contribution of nutrients on the first 7 prinicpal components

ggplot(data=melt_NI) +
   geom_col(aes(x=reorder(rownames(melt_NI), -value), y=value)) +
   theme(axis.text.x = element_text(angle = 90)) +
   labs(x="Nutrients", y="Variable importance")

top6 = rownames(melt_NI)[order(-melt_NI$value)[1:6]]
top6

###
### Biplot
###

slopes = rep(NA, 6)

for(s in c(1:6)) {
   slopes[s] = pr.out$rotation[top6[s],2] / pr.out$rotation[top6[s],1]
}

# Figure 8: Biplot of 6 most significant nutrients impacting the first 6 principal components

ggplot() +
   geom_point(aes(x=pr.out$x[,1], y=pr.out$x[,2])) +
   geom_segment(aes(x=0, y=0, xend=-5, yend=slopes[1]*(-5)), 
                arrow = arrow(length = unit(0.03, "npc")), color="red") +
   annotate("text", x=-4, y=slopes[1]*(-5)+1, label=top6[1], color="red") +
   geom_segment(aes(x=0, y=0, xend=10, yend=slopes[2]*10), 
                arrow = arrow(length = unit(0.03, "npc")), color="red") +
   annotate("text", x=13, y=slopes[2]*10, label=top6[2], color="red") +
   geom_segment(aes(x=0, y=0, xend=12, yend=slopes[3]*12), 
                arrow = arrow(length = unit(0.03, "npc")), color="red") +
   annotate("text", x=12, y=slopes[3]*12-1, label=top6[3], color="red") +
   geom_segment(aes(x=0, y=0, xend=15, yend=slopes[4]*15), 
                arrow = arrow(length = unit(0.03, "npc")), color="red") +
   annotate("text", x=17, y=slopes[4]*15+1, label=top6[4], color="red") +
   geom_segment(aes(x=0, y=0, xend=3, yend=slopes[5]*3), 
                arrow = arrow(length = unit(0.03, "npc")), color="red") +
   annotate("text", x=3, y=slopes[5]*3-1, label=top6[5], color="red") +
   geom_segment(aes(x=0, y=0, xend=30, yend=slopes[6]*30), 
                arrow = arrow(length = unit(0.03, "npc")), color="red") +
   annotate("text", x=30, y=slopes[6]*30+2, label=top6[6], color="red") +
   labs(x="PC1", y="PC2")

###
### Interactive plots
###

# Figure 9: Interactive plot of water amount on the first three principal components

library(plotly)

pr.out_df = data.frame(cbind(pr.out$x[,c(1:7)], SRp))

plot_ly(pr.out_df, x =~PC1, y =~PC2, z =~PC3, color = ~Water_.g.,
        type = "scatter3d", mode = "markers", colors=c("blue", "red"),
        marker=list(size=3))

# Figure 10: Interactive plot of total lipid amount on PC2, PC6 and PC7

plot_ly(pr.out_df, x =~PC2, y =~PC6, z =~PC7, color = ~Fiber_TD_.g.,
        type = "scatter3d", mode = "markers", colors=c("blue", "red"),
        marker=list(size=3))

###
### Correlation
###

# Figure 11: Correlation plot of the dataset

library(corrplot)
SRp.cor = cor(SRp)
corrplot(SRp.cor, method="square", type="upper")

###
### Grouping of similar nutrients
###

sim.crit = 0.4
similar = list()
   
for(j in c(1:45)) {
   nutj = rownames(pr.out$rotation)[j]
   
   for(i in c(1:46)) {
      nuti = rownames(pr.out$rotation)[i]
      dist.scores = sum(abs(pr.out$rotation[j,c(1:7)] - pr.out$rotation[i,c(1:7)]))

      if ((dist.scores <= sim.crit) && (dist.scores != 0)) {
         similar[[nutj]] = append(similar[[nutj]], nuti)
      }
   }
}

for(n in c(1:length(similar))) {
   cat(names(similar[n]),":", similar[[n]],"\n")
}
