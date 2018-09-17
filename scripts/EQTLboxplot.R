getEQTLboxplot <- function(stim.name, snp.name, gene.name, gene_response_data, maf, snp_pvalue, jitter = F){
  .libPaths(c(.libPaths(),'/pasteur/homes/bpiaseck/R/x86_64-unknown-linux-gnu-library/3.1/'))
    #library
    require("ggplot2")
    require("rhdf5")
    source("./scripts/theme_publication.R")
    #snp_data = read.table(sprintf("./data/snpData/%s.sel", snp.name), sep = ' ', header = TRUE)
    hdf = h5read("./data/snpData.h5", sprintf("/%s",snp.name))
    snp_data = data.frame(t(hdf$block0_values))
    colnames(snp_data) = hdf$block0_items
    if ("block1_items" %in% names(hdf)) {
      b1_df = data.frame(t(hdf$block1_values))
      if   (is.null(grep("rs", hdf$block1_items))) {snp_data = cbind(snp_data, b1_df)}
      else                                         {snp_data = cbind(b1_df, snp_data)}
      print ("---> block1 found")
    }
    colnames(snp_data) = c("IID", "snp")
    # get snp and expression data
    df.data = dplyr::left_join(gene_response_data, snp_data, by = c("id"="IID"))
    df.data = subset(df.data, !is.na(snp))
    #eqtl.res = data2save$eqtl.res
    
    #get allele data
    load("./data/imputations.minor.major.Rdata")
    all.minor = as.character(subset(IMP, snp == snp.name)$minor)
    all.major = as.character(subset(IMP, snp == snp.name)$major)
    
    ##function to change labels
    #names(pvals) = eqtl.res$stim
    #to_latin_name = as_labeller(c('null'=paste0("NS<br />p=",as.numeric(pvals["null"])),
    #                              'ecoli'=paste0("E.coli<br /> p=",as.numeric(pvals["ecoli"])),
    #                              'bcg'=paste0("BCG<br /> p=",as.numeric(pvals["bcg"])),
    #                              'saureus'=paste0("S.aureus<br /> p=",as.numeric(pvals["saureus"])),
    #                              'seb'=paste0("SEB<br /> p=",as.numeric(pvals["seb"])),
    #                              'candida'=paste0("C.albicans<br /> p=",as.numeric(pvals["candida"])),
    #                              'flu'=paste0("Influenza<br /> p=",as.numeric(pvals["flu"]))))
    to_latin_name = as_labeller(c('null'= "NS",
                                  'ecoli'= "E.coli",
                                  'bcg'= "BCG",
                                  'saureus'= "S.aureus",
                                  'seb'= "SEB",
                                  'candida'= "C.albicans",
                                  'flu'= "Influenza"))
                                  
    #order stimuli
    #df.data$stimulus = factor(df.data$stim, levels = c("null","ecoli","bcg","saureus","seb","candida","flu"))
    df.data$stimulus = factor(c("NS","E.coli","BCG","S.aureus","SEB","C.Albicans","Influenza")[factor(df.data$stim, levels = c("null","ecoli","bcg","saureus","seb","candida","flu"))], levels = c("NS","E.coli","BCG","S.aureus","SEB","C.Albicans","Influenza"))
    p = ggplot(data=df.data, aes(factor(snp), expr))
    #if (jitter){
    #  p = p + geom_boxplot(aes(fill=stimulus),size=0.1,outlier.shape=NA) 
    #  p = p + geom_jitter(position = position_jitter(width = .1), alpha = 0.5, size = 0.01, colour = "grey",pch=20) 
    #} else {
    #  
    #}
    alpha_facet = rep(0.25, 7)
    alpha_facet[which(c("NS", "E.coli", "BCG", "S.aureus", "SEB", "C.Albicans", "Influenza") == stim.name)] = 1
    p = p + geom_boxplot(aes(fill=stimulus, alpha=stimulus),size=0.1,outlier.shape=20,outlier.size=0.01) 
    p = p + facet_grid(. ~ stimulus) + scale_alpha_manual(values=alpha_facet)#, labeller = to_latin_name)
    p = p + scale_fill_Publication() + guides(fill=FALSE)
    #p = p + ggtitle(paste0(snp.name," (MAF=",as.character(maf),")"))
    p = p + xlab("")
    p = p + ylab("")
    p = p + scale_x_discrete(labels = c("1" = paste0(all.minor,all.major), "2"= paste0(all.minor,all.minor), "0" = paste0(all.major,all.major))) 
    p = p + theme_Publication(base_size=18,axis.tick.lwd = 0.1) + theme(strip.text = element_text(face="bold"), strip.background = element_blank())
    return(p)
}

