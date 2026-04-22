library(readxl)

# Import pharmacologic networks

#SSRI:
library(MBNMAdose)
data(ssri)

#GRISELDA
griselda_full <- read_excel("C:/Users/katie/Desktop/Depression-NMA/data/Griselda_Fixed.xlsx",
                            sheet = "DATA",
                            skip = 2)

# Import nonpharmacologic network
Tian <- read_excel("C:/Users/katie/Desktop/Depression-NMA/data/Tian_Dataset.xlsx")


# Noetel dataset: Mixed interventions - NOT FOR SCOPE OF DISSERTATION ##########

#Available from https://osf.io/nzw6u/overview: data ->  cleaned_deprex
##### Noetel <- readRDS("C:/Users/katie/Desktop/Depression-NMA/data/cleaned_deprex.RDS")

#Tian dataset: Nonpharmacologic interventions


