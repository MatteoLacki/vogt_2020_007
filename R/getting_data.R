library(readxl)
library(stringr)
library(data.table)
library(readr)
library(LFQBench2)

Sample_List_2020 = read_excel("data/Sample List 2020.xlsx", 
                              sheet = "Obelix (O)")
colnames(Sample_List_2020) = make.names(colnames(Sample_List_2020))
VOGT = as.data.table(Sample_List_2020[str_detect(Sample_List_2020$File.Text, 'Vogt'),])

VOGT[(nrow(VOGT)-4):nrow(VOGT),]
pre_isoquant = read_csv("data/pre_isoquant.csv")

R = read_wide_report('data/2020-007 samples combined ST_user designed 20200609-111429 reprocessed with samples annotated _quantification_report.xlsx',
                     skip=1, sheet="TOP3 quantification")
R_matteo = read_wide_report('data/2020-007 samples combined ST_user designed 20200609-111429 reprocessed with samples annotated _quantification_report_matteo.xlsx',
                            sheet="TOP3 quantification")

count_nonNAs = function(x) sum(!is.na(x)) 
nonNA_counts = sapply(R_matteo[,10:ncol(R_matteo)], count_nonNAs)
nonNA_counts_noPools = nonNA_counts[!str_detect(names(nonNA_counts), 'Pool')]
names(nonNA_counts_noPools) = str_sub( names(nonNA_counts_noPools) , end=11 )
nonNA_counts_noPools = data.table(File.Text=names(nonNA_counts_noPools), proteins=nonNA_counts_noPools)

VOGT_noPools = VOGT[!str_detect(File.Text, 'Pool')]
VOGT_noPools[,File.Text:=str_sub( File.Text, end=11 )] 
Filename2FileText = VOGT_noPools[,.(Filename, File.Text)]

Filename2FileText[nonNA_counts_noPools, on='File.Text']
nonNA_counts_noPools[Filename2FileText, on='File.Text']


