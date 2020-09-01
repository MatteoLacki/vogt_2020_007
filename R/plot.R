library(ggplot2)
library(data.table)

load('data/results.Rd')

volcano = ggplot(r[!is.na(AoverB) & pval_BH > 0]) +
  geom_point(aes(x=AoverB, y=-log10(pval_BH)), size=.4) +
  scale_x_log10()+
  facet_grid(.~comp)+
  xlab('A / B') +
  ylab('-log10(pvalue[Benjamini & Hochberg])')

lines = data.table(thr=c(.001,.01,.5,.1))
lines[,threshold := ordered(thr)]

volcano + geom_hline(data=lines, aes(yintercept=-log10(thr), color=threshold))
  


