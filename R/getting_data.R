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
# R_matteo = read_wide_report('data/2020-007 samples combined ST_user designed 20200609-111429 reprocessed with samples annotated _quantification_report_matteo.xlsx',
#                             sheet="TOP3 quantification")

read_ss = function(path) as.data.table(read_excel(path, skip=25))
ss = rbind(read_ss("data/Sample Submission Sheet 2020-007_1-48.xlsx"),
           read_ss("data/Sample Submission Sheet 2020-007_49-60.xlsx")[1:12])
sdesc = read_wide_report("data/2020-007 samples combined ST_user designed 20200609-111429 reprocessed with samples annotated _quantification_report.xlsx", 
                         sheet = "project details", skip = 7, drop_na_columns = F)
sdesc = sdesc[,.(workflow_index,replicate_name,sample_description,acquired_name,identified_proteins)]
sdesc$idx = str_sub(sdesc$replicate_name, 10, 12)
ss$tag = str_replace(str_replace(str_sub(ss$`Sample Description`, 1, 4), '-', ''), ',','')
ss = cbind(ss, sdesc[1:60])

comparisons = list( c( "wt-cyto-no glu", "wt-cyto-glu"), 
                    c( "wt-mem-glu",     "wt-cyto-glu"),
                    c( "wt-mem-glu",     "wt-mem-no glu"),
                    c( "Δ-mem-no glu",   "Δ-mem-glu"),
                    c( "wt-mem-no glu",  "Δ-mem-no glu"),
                    c( "wt",  "ko") )
names(comparisons) = paste(sapply(comparisons,'[',1), sapply(comparisons,'[',2), sep=' VS ')
save(ss, comparisons, R, file='data/all_data.Rd')

# count_nonNAs = function(x) sum(!is.na(x)) 
# nonNA_counts = sapply(R[,10:ncol(R)], count_nonNAs)
# nonNA_counts_noPools = nonNA_counts[!str_detect(names(nonNA_counts), 'Pool')]
# names(nonNA_counts_noPools) = str_sub( names(nonNA_counts_noPools) , end=11 )
# nonNA_counts_noPools = data.table(File.Text=names(nonNA_counts_noPools), proteins=nonNA_counts_noPools)
# 
# VOGT_noPools = VOGT[!str_detect(File.Text, 'Pool')]
# VOGT_noPools[,File.Text:=str_sub( File.Text, end=11 )] 
# Filename2FileText = VOGT_noPools[,.(Filename, File.Text)]
# 
# Filename2FileText[nonNA_counts_noPools, on='File.Text']
# nonNA_counts_noPools[Filename2FileText, on='File.Text']

