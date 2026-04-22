
library(readxl)

# Read the Excel file
tian_data <- read_excel("C://Users//katie//Box//Herder_Dissertation_Dose-Response NMA//Tian_Dataset.xlsx")

tian_data3 <- read_excel("C://Users//katie//Box//Herder_Dissertation_Dose-Response NMA//Tian_Dataset2.xlsx")


# View the data
head(tian_data)
str(tian_data)

# Save data files
save(cleaned_deprex, file= "C://Users//katie//Box//Herder_Dissertation_Dose-Response NMA//Noetel.rdata")
save(tian_data, file= "C://Users//katie//Box//Herder_Dissertation_Dose-Response NMA//Tian.rdata")
save(tian_data, file= "C://Users//katie//Box//Herder_Dissertation_Dose-Response NMA//Tian.rdata")