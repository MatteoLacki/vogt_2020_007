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
  if(!no_A & !no_B) return( t.test(A,B,...)$p.value ) # the test neglects NAs: not good!
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

p_adjust = function(p, ...){
  r = p
  correction_needed = p > 0 & p < 1 & !is.na(p)
  r[correction_needed] = p.adjust(r[correction_needed], ...)
  return(r)
}

get_pvals = function(comp, Rw, test, ratio, adjustments=c('bonferroni','BH','BY')){
  comp = c('accession',comp)
  Rw_comp = Rw[,..comp]
  colnames(Rw_comp) = c('accession','A','B')
  res = Rw_comp[,.(pval = test(A,B), AoverB = ratio(A,B)), by=.(accession)]
  res[, (adjustments) := lapply(adjustments, function(m) p_adjust(pval, method=m)) ]
  return(res)
}
adjustments=c('bonferroni','BH','BY')
longify = function(r, adjustments=c('bonferroni','BH','BY')) 
  dcast(r, accession~comp, value.var = c(adjustments, 'AoverB'), sep=": ")[order(accession)]

r1 = rbindlist(lapply(comparisons, get_pvals, Rw, TenzerTest, ratio), idcol='comp')
save(r1, file='data/results1.Rd')
rw1 = longify(r1)
write.csv(rw1, file='data/comparison1.csv')


# another version
TenzerTestSafer = function(A, B, cutoff=.5, bogus_pval=1, highly_significant_pvalue=0, ...){
  A_NA_cnt = sum(is.na(A))
  B_NA_cnt = sum(is.na(B))
  no_A = A_NA_cnt / length(A) >= cutoff
  no_B = B_NA_cnt / length(B) >= cutoff

  if(length(A)-A_NA_cnt < 2 & length(B)-B_NA_cnt < 2) return(bogus_pval)
  if(A_NA_cnt == length(A) & B_NA_cnt == length(B)) return( t.test(A,B,...)$p.value )
  if(no_A & B_NA_cnt == length(B)) return(highly_significant_pvalue)
  if(no_B & A_NA_cnt == length(A)) return(highly_significant_pvalue)
  return(bogus_pval)
}

r2 = rbindlist(lapply(comparisons, get_pvals, Rw, TenzerTestSafer, ratio), idcol='comp')
save(r2, file='data/results2.Rd')
rw2 = longify(r2)
write.csv(rw2, file='data/comparison1.csv')  


# New procedure: 
#   infer small values if there is at least one small observed or all of them are small.
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
  if(all(A == A[1]) & all(B == B[1])) return(-3) # All outcomes are the same: t-test will fail.
  if(length(A) < 2) return(-4) # impossible to estimate st-dev in A: t-test will fail
  if(length(B) < 2) return(-5) # impossible to estimate st-dev in B: t-test will fail
  return( t.test(A,B,...)$p.value )
}
ratio2 = function(A, B) median(A, na.rm=T) / median(B, na.rm=T)

r3 = rbindlist(lapply(comparisons, get_pvals, RwFix, TenzerTest2, ratio2), idcol='comp')
save(r3, file='data/results3.Rd')
rw3 = longify(r3)
write.csv(rw3, file='data/comparison3.csv')


# further comparisons of the two approaches
r1$method = 'old'
r2$method = 'old_safer'
r3$method = 'new'
rr = rbind(r1, r2, r3)
rr = dcast(rr, comp+accession~method, value.var = c(adjustments,'pval', 'AoverB'))

# with(rr, plot(pval_old_safer, pval_new))
# with(rr[pval_new>=0], plot(pval_old, pval_new))
# with(rr, plot(AoverB_old, AoverB_new, asp=1))

NAs = Rl[,.(NAcnt=sum(is.na(I)), OBScnt=length(I)),by=.(Group,accession)]
plot(NAs[, .N, by=NAcnt][order(NAcnt)], type='b', ylim=c(0,5000)) # relatively few observations 

get_NA_stats = function(comp, Rw){
  comp = c('accession',comp)
  Rw_comp = Rw[,..comp]
  colnames(Rw_comp) = c('accession','A','B')
  res = Rw_comp[,.(A_NA_cnt=sum(is.na(A)), B_NA_cnt=sum(is.na(B))),by=.(accession)]
  return(res[,.N,by=.(A_NA_cnt, B_NA_cnt)])
}
# here be counts of 
NA_stats = rbindlist(lapply(comparisons, get_NA_stats, Rw), idcol='comp')
save(NAs, NA_stats, file='data/NA_stats.Rd')


# for corporation partners:
summary(r1)
any(is.na(r1$AoverB)) # no NAs
any(r1$AoverB == Inf) # no infnities
any(r1$AoverB == 0) # no infnities
sum(r1$AoverB == 1) # No of phoney replacements.
hist(log10(r1[AoverB!=1,AoverB]), breaks=100, main='A/B')

rw1_final = dcast(r1, accession~comp, value.var = c('BH', 'AoverB'), sep=": ")[order(accession)]
write.csv(rw1_final, file='data/comparison1_final.csv')

