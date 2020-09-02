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

TenzerTest = function(A, B, cutoff=.5,...){
  no_A = sum(is.na(A))/length(A) >= cutoff
  no_B = sum(is.na(B))/length(B) >= cutoff
  if(!no_A & !no_B) return( t.test(A,B,...)$p.value )
  if(no_A & no_B) return(1)
  return(0)
}

Rl[,.(p=c(0,.001), q=quantile(I, probs=c(0,.001), na.rm=T)), by=Group]

quantile(Rl$I, probs=c(0,.001), na.rm = T)

TenzerTest = function(A, B, cutoff=.5, fill_in_value=1,...){
  A_NA_cnt = sum(is.na(A))
  B_NA_cnt = sum(is.na(B))

  no_A = A_NA_cnt/length(A) >= cutoff
  no_B = B_NA_cnt/length(B) >= cutoff
  
  if(no_A){ A[is.na(A)] = fill_in_value }
  if(no_B){ B[is.na(B)] = fill_in_value }
  
  if(no_A & no_B) return(NA)
  if(A_NA_cnt > 0 & !no_A) return(NA) # some missing observations in A, but not enough to say all are 0.
  if(B_NA_cnt > 0 & !no_B) return(NA) # some missing observations in A, but not enough to say all are 0.
  
  
  return( t.test(A,B,...)$p.value )
}

X = c(NA, NA, 5, 100)
X_small_I = 10
X_fill_in_value = 1

fill_values = function(X, X_small_I, X_fill_in_value){
  X_NA_cnt = sum(is.na(X))
  if(X_NA_cnt == length(X)) return(replicate(length(X), X_fill_in_value))
  if(X[!is.na(X)] <= X_small_I){
    X[is.na(X)] = X_fill_in_value
  }
  return(X)
}

TenzerTest2 = function(A, B, A_small_I, B_small_I, A_fill_in_value, B_fill_in_value, ...){
  A = fill_values(A)
  B = fill_values(B)
  
  if(any(is.na(A))) return(NA) # some NAs persist: bogus
  if(any(is.na(B))) return(NA) # some NAs persist: bogus
  if(all(A == A[1]) & all(B == B[1])) return(1) # All outcomes are the same: t-test will fail.
  if(length(A) < 2) return(NA) # impossible to estimate st-dev for A: t-test will fail
  if(length(B) < 2) return(NA) # impossible to estimate st-dev for A: t-test will fail
  
  return( t.test(A,B,...)$p.value )
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

ratio(c(NA,NA,NA), 1:3)
ratio(1:3, c(NA,NA,NA))
ratio(1:3, c(NA,1,1))

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
save(r, file='data/results.Rda')
rw = dcast(r, accession~comp, value.var = c('pval_bon','pval_BY', 'pval_BH', 'AoverB'), sep=": ")
rw = rw[order(accession)]
write.csv(rw, file='data/comparison.csv')

