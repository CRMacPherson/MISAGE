g1 = ggplot(subset(df2,gene == "IL13"),aes(x=age,y=expr,color=stim))  
g1 = g1 + geom_point(size=0.7,alpha=0.1,pch=16) #alpha=0.5
g1 = g1 + facet_grid(~stim)
g1 = g1 + stat_smooth(method = "lm",size=0.5)
g1 = g1 + scale_x_continuous(breaks = c(20,40,60))
g1 = g1
plot(g1)
