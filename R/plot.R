library(ggplot2)
library(data.table)
library(patchwork)
library(ggeasy)

load('data/results1.Rd')
load('data/results2.Rd')
load('data/results3.Rd')

erupt_volcano = function(r) ggplot(r) +
    geom_point(aes(x=AoverB, y=-log10(pval_BH)), size=.4) +
    scale_x_log10()+
    facet_grid(.~comp)+
    xlab('A / B') +
    ylab('-log10(pvalue[Benjamini & Hochberg])') 

thr_lines = function(thr=c(.001,.01,.05,.1)){
  lines = data.table(thr=thr)
  lines[,threshold := ordered(thr)]
  geom_hline(data=lines, aes(yintercept=-log10(thr), color=threshold))
}

r_good = r[!is.na(AoverB) & pval_BH > 0 & pval_BH < 1 & !is.na(pval_BH)]
r_good$method = 'old'
r2_good = r2[!is.na(AoverB) & pval_BH > 0 & pval_BH < 1 & !is.na(pval_BH)]
r2_good$method = 'new'

r1$method='old'
r2$method='old_safer'
r3$method='new'

erupt_volcano(rbind(r1, r3)) + facet_grid(method~comp) + thr_lines()


