library(stringr)
library(data.table)

load('data/all_data.Rd')
# ss$replicate name corresponds to columns in R
desc = colnames(R)[1:9]

x = c('accession', ss$replicate_name)
Rl = melt(R[,..x], id.vars = 'accession', variable.name = 'cond', value.name = "I")
Rl = merge(Rl, ss[,.(replicate_name, Group)], by.x='cond', by.y='replicate_name')
Rl = Rl[,-1]
Rl = Rl[,.(I=list(I)), by=.(accession, Group)]
Rw = dcast(Rl, accession~Group, value.var = 'I')

TenzerTest = function(A, B, cutoff=.5,...){
  no_A = sum(is.na(A))/length(A) >= cutoff
  no_B = sum(is.na(B))/length(B) >= cutoff
  if(!no_A & !no_B) return( t.test(A,B,...)$p.value )
  if(no_A & no_B) return(1)
  return(0)
}
TenzerTestLists = function(A,B,cutoff=.5,...) TenzerTest(A[[1]], B[[1]], cutoff, ...) 

ratio = function(A,B){
  no_A = sum(is.na(A))/length(A) >= .5
  no_B = sum(is.na(B))/length(B) >= .5
  med_A = median(A, na.rm=T)
  med_B = median(B, na.rm=T)
  med_A / med_B
}
ratioLists = function(A,B) ratio(A[[1]], B[[1]])

get_pvals = function(comp){
  comp = c('accession',comp)
  Rw_comp = Rw[,..comp]
  colnames(Rw_comp) = c('accession','A','B')
  res = Rw_comp[,.(pval=TenzerTestLists(A, B),
                   AoverB=ratioLists(A,B)), by=.(accession)]
  res[,`:=`(pval_bon=p.adjust(pval, method = 'bonferroni'), pval_BH=p.adjust(pval, method = 'BH'), pval_BY=p.adjust(pval, method = 'BY'))]
  res
}

r = rbindlist(lapply(comparisons, get_pvals), idcol='comp')
save(r, file='data/results.Rda')



