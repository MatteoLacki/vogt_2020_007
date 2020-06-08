library(readxl)
library(stringr)
library(data.table)

Sample_List_2020 <- read_excel("data/Sample List 2020.xlsx", 
                               sheet = "Obelix (O)")
View(Sample_List_2020)
colnames(Sample_List_2020) = make.names(colnames(Sample_List_2020))
VOGT = as.data.table(Sample_List_2020[str_detect(Sample_List_2020$File.Text, 'Vogt'),])

View(VOGT)
