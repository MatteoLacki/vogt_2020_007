library(ggplot2)
library(data.table)
library(patchwork)
library(ggeasy)

load('data/results1.Rd')
load('data/results2.Rd')
load('data/results3.Rd')
load('data/NA_stats.Rd')

erupt_volcano = function(r) ggplot(r) +
    geom_point(aes(x=AoverB, y=-log10(BH)), size=.4) +
    scale_x_log10() +
    facet_grid(.~comp)+
    xlab('A / B') +
    ylab('-log10(pvalue[Benjamini & Hochberg])') 

thr_lines = function(thr=c(.001,.01,.05,.1)){
  lines = data.table(thr=thr)
  lines[,threshold := ordered(thr)]
  geom_hline(data=lines, aes(yintercept=-log10(thr), color=threshold))
}

r1$method='old'
r2$method='old_safer'
r3$method='new'

erupt_volcano(rbind(r1, r3)) + facet_grid(method~comp) + thr_lines()

ggplot(NA_stats) + geom_tile(aes(x=A_NA_cnt, y=B_NA_cnt, fill=N)) + facet_wrap(~comp)
ggplot(NA_stats[A_NA_cnt>0 | B_NA_cnt>0]) + geom_tile(aes(x=A_NA_cnt, y=B_NA_cnt, fill=N)) + facet_wrap(~comp)

# for corporation partners:
erupt_volcano(r1) + thr_lines()
