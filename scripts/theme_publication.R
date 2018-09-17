theme_Publication <- function(base_size=8, base_family="Helvetica", 
                              border.col = NA, border.size = 1.5,
                              grid.major.color = "white",
                              axis.tick.lwd = 0.5,
                              legend.pos = "bottom",
                              legend.dir = "horizontal",
                              legend.size = 0.3) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
  + theme(plot.title = element_text(#face = "bold",
                                    size = rel(1), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = border.col, size = border.size),
          axis.title = element_text(size = rel(1)),#face = "bold"
          axis.title.y = element_text(angle=90,vjust = 2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          #axis.line.x = element_line(colour="black"),
          #axis.line.y = element_line(colour="black"),
          axis.ticks = element_line(size=axis.tick.lwd),
          panel.grid.major = element_line(colour=grid.major.color,size=0.2),  #element_blank(), 
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = legend.pos,
          legend.direction = legend.dir,
          legend.key.size= unit(legend.size, "cm"),
          legend.margin = unit(0, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(5,5,5,5),"mm"),
          strip.background=element_rect(colour="grey20",fill="white"),
          #strip.background=element_rect(colour="white",fill="white"),
          #strip.text = element_text(face="bold",colour="grey20")
          strip.text = element_text(colour="grey20",size=rel(0.6))
  ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#999999","#377EB8","gold2","#4DAF4A","#984EA3","#E41A1C","#FF7F00")), ...)
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#999999","#377EB8","gold2","#4DAF4A","#984EA3","#E41A1C","#FF7F00")), ...)
  
}

scale_fill_Publication_SI <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#377EB8","gold2","#4DAF4A","#984EA3","#E41A1C","#FF7F00")), ...)
  
}

scale_colour_Publication_SI <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#377EB8","gold2","#4DAF4A","#984EA3","#E41A1C","#FF7F00")), ...)
  
}