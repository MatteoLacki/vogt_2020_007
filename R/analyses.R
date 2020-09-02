# Overall task: compare the intensities between conditions.
library(stringr)
library(data.table)

load('data/all_data.Rd')
# ss$replicate name corresponds to columns in R
desc = colnames(R)[1:9]

x = c('accession', ss$replicate_name)
Rl = melt(R[,..x], id.vars = 'accession', variable.name = 'cond', value.name = "I")
Rl = merge(Rl, ss[,.(replicate_name, Group)], by.x='cond', by.y='replicate_name')
Rl$repl = as.integer(str_sub(Rl$cond, -1))
Rl = Rl[,-1]
Rw = dcast(Rl, accession+repl~Group, value.var = 'I')

TenzerTest = function(A, B, cutoff=.5, bogus_pval=1, highly_significant_pvalue=0, ...){
  no_A = sum(is.na(A))/length(A) >= cutoff
  no_B = sum(is.na(B))/length(B) >= cutoff
  if(!no_A & !no_B) return( t.test(A,B,...)$p.value )
  if(no_A & no_B) return(bogus_pval)
  return(highly_significant_pvalue)
}
ratio = function(A, B, fill_in_value=1){
  no_A = sum(is.na(A))/length(A) >= .5
  no_B = sum(is.na(B))/length(B) >= .5
  if(no_A){ A[is.na(A)] = fill_in_value }
  if(no_B){ B[is.na(B)] = fill_in_value }
  med_A = median(A, na.rm=T)
  med_B = median(B, na.rm=T)
  med_A / med_B
}

get_pvals = function(comp, Rw){
  comp = c('accession',comp)
  Rw_comp = Rw[,..comp]
  colnames(Rw_comp) = c('accession','A','B')
  res = Rw_comp[,.(pval = TenzerTest(A,B),
                   AoverB = ratio(A,B)),
                by=.(accession)]
  res[,`:=`(pval_bon= p.adjust(pval, method = 'bonferroni'),
            pval_BH = p.adjust(pval, method = 'BH'),
            pval_BY = p.adjust(pval, method = 'BY'))]
  res
}

r = rbindlist(lapply(comparisons, get_pvals, Rw=Rw), idcol='comp')
save(r, file='data/results.Rd')
rw = dcast(r, accession~comp, value.var = c('pval_bon','pval_BY', 'pval_BH', 'AoverB'), sep=": ")
rw = rw[order(accession)]
write.csv(rw, file='data/comparison.csv')

# Alternative procedure


fill_values = function(X, X_small_I, X_fill_in_value){
  X_NA_cnt = sum(is.na(X))
  if(X_NA_cnt == length(X)) return(replicate(length(X), X_fill_in_value))
  if(any(X[!is.na(X)] <= X_small_I)){
    X[is.na(X)] = X_fill_in_value
  }
  return(X)
}

fill_low_values = function(){
  RwFix = Rw[,.(accession, repl)]
  conditions = unique(unlist(comparisons, use.names = F))
  p = .001
  for(cond in conditions){
    cols = c('accession', cond)
    x = Rw[,..cols]
    colnames(x) = c('accession', 'X')
    x_small_I = quantile(x$X, p, na.rm=T)
    x_fill_in_value = min(x$X, na.rm=T)
    x[,fX:=fill_values(X, x_small_I, x_fill_in_value),by=accession]
    set(RwFix, j=cond, value=x$fX)  
  }
  return(RwFix)
}

RwFix = fill_low_values()

TenzerTest2 = function(A, B, ...){
  if(any(is.na(A))) return(-1) # some NAs persist: bogus
  if(any(is.na(B))) return(-2) # some NAs persist: bogus
  if(all(A == A[1]) & all(B == B[1])) return(1) # All outcomes are the same: t-test will fail.
  if(length(A) < 2 | length(B) < 2) return(-3) # impossible to estimate st-dev: t-test will fail
  
  return( t.test(A,B,...)$p.value )
}

ratio2 = function(A, B) median(A, na.rm=T) / median(B, na.rm=T)

get_pvals2 = function(comp, RwFix){
  comp = c('accession',comp)
  Rw_comp = RwFix[,..comp]
  colnames(Rw_comp) = c('accession','A','B')
  res = Rw_comp[,.(pval = TenzerTest2(A,B),
                   AoverB = ratio2(A,B)),
                by=.(accession)]
  res[,`:=`(pval_bon= p.adjust(pval, method = 'bonferroni'),
            pval_BH = p.adjust(pval, method = 'BH'),
            pval_BY = p.adjust(pval, method = 'BY'))]
  res
}

r2 = rbindlist(lapply(comparisons, get_pvals2, RwFix), idcol='comp')
save(r2, file='data/results2.Rd')
rw2 = dcast(r2, accession~comp, value.var = c('pval_bon','pval_BY', 'pval_BH', 'AoverB'), sep=": ")
rw2 = rw[order(accession)]
write.csv(rw2, file='data/comparison2.csv')





