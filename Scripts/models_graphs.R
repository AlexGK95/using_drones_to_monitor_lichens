## This script runs generalized linear models for lichen abundance and relative coverage
## using species, relative elevation and their interaction as predictors. It also produces
## a pixel abundance-relative elevation and relative coverage-relative elevation graphs

### The first section focuses on pixel abundance-relative elevation

setwd("working_directory") # set working directory
getwd()
require(bigmemory) # use bigmemory to upload large files
require(ggplot2)

## import data

A1= read.big.matrix("data_path") # read All_p_elev.csv file
Pink= as.matrix(A1) # convert from big.matrix to matrix
A2= read.big.matrix("data_path") # read All_r_elev.csv file
Red= as.matrix(A2)
A3= read.big.matrix("data_path") # read All_b_elev.csv file
Black= as.matrix(A3)

## export pixel abuncance-relative elevation graph 

pdf("output_path", 11.7, 8.3) # export graph
plot(p2, col= rgb(1, 0, 0, 0.8), xlab="Relative elevation", ylab= "Pixel abundance", main="Relationship between lichen abundance and relative elevation", cex.lab=1.4, cex.main= 2)
legend("topleft", c("Red", "Pink", "Black"), fill=c("red", rgb(0, 0, 1, 1/2), "black"),  bty="n", title= expression(bold("Species")), y.intersp= 0.7, cex=1.4)
plot(p1, col= rgb(0, 0, 1, 1/2), add=T)
plot(p3, col= rgb(0, 0, 0, 1), add=T)
dev.off()

## isolate elevation values and their frequencies

Pink_elev= unique(Pink) # unique elevation values
Red_elev= unique(Red)
Black_elev= unique(Black)
All_elev= c(sort(Pink_elev, decreasing =TRUE), sort(Red_elev, decreasing  =TRUE), sort(Black_elev, decreasing  =TRUE)) # sort elevation values appropriately to fit with frequency values

Pink_freq=table(Pink) # frequency values
Pink_freq= rev(as.numeric(Pink_freq)) # isolate frequencies and reverse to fit with elevation values
Red_freq=table(Red) 
Red_freq= rev(as.numeric(Red_freq))
Black_freq=table(Black) 
Black_freq= rev(as.numeric(Black_freq))
All_freq= c(rev(Pink_freq), Red_freq, Black_freq)

## create dataframe

Species= c(rep("Pink", 163), rep("Red", 164), rep("Black", 164)) # create species column, repeating each species as many times it's represented in the data
All_data= data.frame(Species, All_elev, All_freq)

## run generalized linear model

model1= glm(All_freq~Species*scale(All_elev), data=All_data, family=poisson(link=log)) # response= pixel abundance, predictors= species*relative elevation, link funtion= log
Summary= summary(model1) 
Summary$coefficients
(model2$null.deviance-model2$deviance)/model2$null.deviance # proportion of deviance explained
par(mfrow=c(2,2))
plot(model1) # inspect residual plots
Coeff_back= exp(Summary$coefficients[,1]) # back-transform coefficients for interpretation
Results= cbind(Summary$coefficients, Coeff_back) # combine results and back-transformed parameters into dataframe
write.csv(Results,"output_path") # export results


### The second section focuses on relative coverage-relative elevation

rm(list=ls())

# import data

Data= read.csv("data_path", header= T) # read All_gradient_rpb.csv file

# produce relative coverage-relative elevation graph 

Data$Relative_altitude_perc= c("0-5%", "5-10%", "10-15%", "15-20%", "20-25%", "25-30%", "30-35%","35-40%", "40-45%","45-50%", "50-55%", "55-60%", "60-65%", "65-70%", "70-75%", "75-80%",  "80-85%", "85-90%", "90-95%", "95-100%", "0-5%", "5-10%", "10-15%", "15-20%", "20-25%", "25-30%", "30-35%","35-40%", "40-45%","45-50%", "50-55%", "55-60%", "60-65%", "65-70%", "70-75%", "75-80%",  "80-85%", "85-90%", "90-95%", "95-100%", "0-5%", "5-10%", "10-15%", "15-20%", "20-25%", "25-30%", "30-35%","35-40%", "40-45%","45-50%", "50-55%", "55-60%", "60-65%", "65-70%", "70-75%", "75-80%",  "80-85%", "85-90%", "90-95%", "95-100%") # add percentage ranges
Data$Relative_altitude_perc= as.factor(Data$Relative_altitude_perc) # convert to factor variable
Data$Relative_altitude_perc= factor(Data$Relative_altitude_perc, levels(Data$Relative_altitude_perc)[c(1,10,2:9, 11:20)]) # rearrange proportion ranges appropriately, so that each class is in the correct position numerically
pdf("output_path", 11.8, 7.3) # export graph
ggplot(Data, aes(fill=Species, y=Relative_abundance, x=Relative_altitude_perc))+geom_bar(position="dodge", stat="identity")+theme_bw()+scale_fill_manual(values = c("black", "purple", "red"))+xlab("Relative Altitude")+ylab("Relative coverage")+ggtitle("Relationship between relative lichen coverage and relative altitude")+scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))+theme(plot.title = element_text(hjust = 0.5))+theme(plot.title = element_text(size = 25, face = "bold"))+theme(axis.title= element_text(size=15))+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

# run generalized linear model

Data$Relative_altitude= as.numeric(Data$Relative_altitude) # convert relative altitude to numeric variable, to run model
Data$Relative_altitude= rep(seq(0.05,1,0.05), 3) # values must be adjusted to range from 0.05 to 1 for each species, as they change after conversion from factor to numeric
model2= glm(Data$Relative_abundance~Data$Species*scale(Data$Relative_altitude), family= quasibinomial(link=logit))
Summary2= summary(model2)
Summary2$coefficients
(model2$null.deviance-model2$deviance)/model2$null.deviance # proportion of deviance explained
par(mfrow=c(1,2))
plot(model2, which=c(1,2)) # inspect residual plots
Coeff_back2= (Summary2$coefficients[,1])/6 # back-transform coefficients for interpretation (divide by 6 as a rule of thumb)
Results2= cbind(Summary2$coefficients, Coeff_back2) # combine results and back-transformed parameters into dataframe
write.csv(Results2,"output_path") # export results
