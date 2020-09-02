library(ggplot2)
library(data.table)
library(patchwork)
library(ggeasy)

load('data/results.Rd')
load('data/results2.Rd')

erupt_volcano = function(r, thr=c(.001,.01,.05,.1)){
  volcano = ggplot(r[!is.na(AoverB) & pval_BH > 0]) +
    geom_point(aes(x=AoverB, y=-log10(pval_BH)), size=.4) +
    scale_x_log10()+
    facet_grid(.~comp)+
    xlab('A / B') +
    ylab('-log10(pvalue[Benjamini & Hochberg])') 
  
  lines = data.table(thr=thr)
  lines[,threshold := ordered(thr)]
  volcano + geom_hline(data=lines, aes(yintercept=-log10(thr), color=threshold))
}

(erupt_volcano(r) + ggtitle('Old Method')) / (erupt_volcano(r2) + ggtitle('New Method'))
