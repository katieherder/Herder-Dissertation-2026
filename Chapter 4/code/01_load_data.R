################################################################################
# Chapter 4: Load Data
################################################################################

library(readxl)
setwd("C:/Users/katie/Desktop/Depression-NMA/Chapter 4")

#SSRI:
library(MBNMAdose)
data(ssri)

#GRISELDA
griselda_full <- read_excel("C:/Users/katie/Desktop/Depression-NMA/Chapter 4/data/Griselda_Fixed.xlsx",
                            sheet = "DATA",
                            skip = 2)
