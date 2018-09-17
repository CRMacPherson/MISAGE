shinyServer(function(input, output, session) {
  ########### > REACTIVE VALUES ----
  rv = shiny::reactiveValues()
  rv$clickedGene = NULL
  rv$clickedGeneCount = 0
  rv$heatmapClicked = NULL
  ########### > NON-REACTIVE VALUES ----
  #- nv must always be used with shiny::isolate()
  nv = shiny::reactiveValues()
  nv$heatmapClicked = NULL
  ########### > FUNCTIONS ----
  getCov = shiny::reactive({
    load("./data/Results_expr.Rdata") #- loads data.frame `res` containing age and sex effects (beta, p-value, fdr)
    return (res)
  })
  getCovSex = shiny::reactive({
    res = getCov()
    thr.cs = 0.5 #threshold for gene scores
    thr.fdr = 0.01 #threshold for fdr
    cs = res[,grep("cs",colnames(res))]
    colnames(cs) = gsub("cs[.]","",colnames(cs))
    cs = ifelse(cs < thr.cs,0,1)
    cov = "sex"
    cov.name = toupper(cov)
    cov.data = res[,grep(paste0(cov.name,".*global[.]q"),colnames(res))]
    colnames(cov.data) = c("NS","E.coli","BCG","S.aureus","SEB","C.albicans","Influenza")
    cov.data = ifelse(cov.data > thr.fdr,0,1)  
    cov.data = cov.data * cs
    #select only genes that depend on age/sex in at least one condition
    #cov.data = cov.data[rowSums(cov.data) > 0,]
    #sort by specificity
    sort.index = order(128*rowSums(cov.data) + 64*cov.data[,1] + 32*cov.data[,2] + 16*cov.data[,3] + 8*cov.data[,4] + 4*cov.data[,5] + 2*cov.data[,6] + 1*cov.data[,7])
    cov.data = (cov.data*rowSums(cov.data))[sort.index, ]
    return (cov.data)
  })
  getCovAge = shiny::reactive({
    res = getCov()
    thr.cs = 0.5 #threshold for gene scores
    thr.fdr = 0.01 #threshold for fdr
    cs = res[,grep("cs",colnames(res))]
    colnames(cs) = gsub("cs[.]","",colnames(cs))
    cs = ifelse(cs < thr.cs,0,1)
    cov = "age"
    cov.name = toupper(cov)
    cov.data = res[,grep(paste0(cov.name,".*global[.]q"),colnames(res))]
    colnames(cov.data) = c("NS","E.coli","BCG","S.aureus","SEB","C.albicans","Influenza")
    cov.data = ifelse(cov.data > thr.fdr,0,1)  
    cov.data = cov.data * cs
    #select only genes that depend on age/sex in at least one condition
    #cov.data = cov.data[rowSums(cov.data) > 0,]
    #sort by specificity
    sort.index = order(128*rowSums(cov.data) + 64*cov.data[,1] + 32*cov.data[,2] + 16*cov.data[,3] + 8*cov.data[,4] + 4*cov.data[,5] + 2*cov.data[,6] + 1*cov.data[,7])
    cov.data = (cov.data*rowSums(cov.data))[sort.index, ]
    return (cov.data)
  })
  getCovNames = shiny::reactive({
    stimulusName = nv$heatmapClicked
    #- Get AGE data
    AGE = getCovAge()
    gene.group.index = apply(AGE > 0,1,sum) == as.integer(input$specifiicity_plot_selector)
    AGE = AGE[gene.group.index, stimulusName]
    AGE = names(AGE[AGE > 0])
    #- Get SEX data
    SEX = getCovSex()
    gene.group.index = apply(SEX > 0,1,sum) == as.integer(input$specifiicity_plot_selector)
    SEX = SEX[gene.group.index, stimulusName]
    SEX = names(SEX[SEX > 0])
    #- Make table
    COMBINED = sort(unique(c(SEX,AGE)))
    return (list(combined = COMBINED, sex = SEX, age = AGE, stimulusName = stimulusName))
  }) 
  getLabeler = shiny::reactive({
    labeler = as_labeller(c('null'="NS",
                            'ecoli'="E.coli",
                            'bcg'="BCG",
                            'saureus'="S.aureus",
                            'seb'="SEB",
                            'candida'="C.albicans",
                            'flu'="Influenza"))
    return (labeler)
  })
  getSigAssociations = shiny::reactive({
    load("./data/GeneSnpSignificantAssociations.Rdata")
    GENE.SNP.SIG.ASSOCIATIONS$MAF = round(GENE.SNP.SIG.ASSOCIATIONS$MAF,2)
    GENE.SNP.SIG.ASSOCIATIONS$Pvalue = round(GENE.SNP.SIG.ASSOCIATIONS$Pvalue, 15)
    return (GENE.SNP.SIG.ASSOCIATIONS)
  })
  getSigAssociationsByGene = shiny::reactive({
    data = subset(getSigAssociations(), Gene == rv$clickedGene)
    return (data)
  })
  getGeneResponseData = shiny::reactive({
    load("./data/gene_response_age.corrected.Rdata")
    return (df)
  })
  getVarianceExplained = shiny::reactive({
    load("./data/variance.explained.Rdata")
    return (variance.explained)
  })
  getRefSeqData = shiny::reactive({
    load("./data/refseq.Rdata")
    return (refseq)
  })
  getRefSeq = function(symbol) {
    data = getRefSeqData()
    return (as.character(data[data$Symbol == symbol, "RefSeq"]))
  }
  getGeneCard = shiny::reactive({
    refseq = getRefSeq(rv$clickedGene)
    data = mygene::getGene(refseq, fields = "all", scopes = "refseq")
    if (is.null(names(data))) {data = data[[1]]}
    return (data)
  })
  getGeneFeature = function(x) {
    data = getGeneCard()
    output = data[[x]]
    if (is.null(output)) {return ("NA")}
    else {return (output)}
  }
  triggerGeneModal = function() {
    rv$clickedGeneCount = rv$clickedGeneCount + 1
  }
  getPresenceQColor = function(value) {
    if (value > 0.9) {return ("green")}
    else if (value > 0.5) {return ("yellow")}
    else {return ("red")}
  }
  getChromLength = shiny::reactive({
    load("./data/chr_length.Rdata")
    return (chr.l)
  })
  ########### > PLOTS ----
  ########### >> PCA ----
  ###########
  output$pca = plotly::renderPlotly({
    load("./data/data_7stim.pca.age.sex.Rdata")
    mypc = dplyr::sample_frac(mypc, input$pca_subsample/100, replace = FALSE)
    mypc$group[mypc$group == "Null"] = "NS"
    mypc$group[mypc$group == "Candida"] = "C.albicans"
    mypc$group[mypc$group == "Flu"] = "Influenza"
    mypc$SEX = "Male"; mypc$SEX[mypc$sex == "F"] = "Female"
    axis = list(showgrid = TRUE, zeroline = FALSE, showline = FALSE, showticklabels = FALSE)
    scene = list(xaxis = c(axis,title="PC-2"), yaxis = c(axis,title="PC-3"), zaxis = c(axis,title="PC-1"), camera = list(eye = list(x = 0.5, y = -2.25, z = 0.5)))
    p = plotly::plot_ly(data = mypc, x = ~PC2, y = ~PC3, z = ~PC1, hoverinfo = 'text', source = "pca")
    if (input$pca_grouping == "Stimulus") {
      p = p %>% plotly::add_markers(marker = list(size=4), color = ~factor(group, levels = c("NS", "E.coli", "BCG", "S.aureus", "SEB", "C.albicans", "Influenza")), colors = c("NS" = "#999999", "E.coli" = "#377EB8", "BCG" = "gold2", "S.aureus" = "#4DAF4A", "SEB" = "#984EA3", "C.albicans" = "#E41A1C", "Influenza" = "#FF7F00"), text = ~sprintf("%s (%s aged %s-%s)", group, SEX, age, age+10))
    }
    else if (input$pca_grouping == "Sex") {
      p = p %>% plotly::add_markers(marker = list(size=4), color = ~SEX, colors = c("#1f78b4","#ff7f00"), text = ~sprintf("%s (%s aged %s-%s)", group, SEX, age, age+10))
    }
    else {
      p = p %>% plotly::add_markers(marker = list(size=4), color = ~as.factor(age), colors = c("#C51B7D", "#F1B6DA", "#D1D1D1", "#B8E186", "#4D9221"), text = ~sprintf("%s (%s aged %s-%s)", group, SEX, age, age+10)) #c("#a1d99b","#74c476","#41ab5d","#238b45","#006d2c")
    }
    p = p %>% plotly::layout(showlegend = TRUE, scene = scene, margin = list(b=0,l=0,r=0,t=0)) %>%
        plotly::config(showLink = FALSE, displayModeBar = FALSE)
  })
  ########### >> SPECIFICITY ----
  ###########
  get_specificity_pie = function (data, variable.name, selected, gene_selected, force_color = NULL) {
    if (gene_selected != "or, select gene") {
      if (variable.name == "AGE") {
        if (input$specifiicity_plot_geneFinder %in% rownames(getCovAge())) {selected = as.character(sum(getCovAge()[input$specifiicity_plot_geneFinder, ] > 0))}
        else {selected = "Select gene group"}
      }
      if (variable.name == "SEX") {
        if (input$specifiicity_plot_geneFinder %in% rownames(getCovAge())) {selected = as.character(sum(getCovSex()[input$specifiicity_plot_geneFinder, ] > 0))}
        else {selected = "Select gene group"}
      }
    }
    data = data.frame(table(apply(data > 0,2,sum)))
    colnames(data) = c("Group", "Freq")
    if (selected != "0") {
      data = subset(data, Group != 0)
      force_color = force_color[2:length(force_color)]
    }
    data$labels = sprintf("<b>Group %s</b><br />%s genes affected<br />by %s in %s stimuli", data$Group, data$Freq, variable.name, data$Group)
    if (is.null(force_color)) {
      data$colors = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf")
    } else {
      data$colors =  force_color
    }
    if (selected != "Select gene group") {
      selected = as.numeric(selected)
      pull = 0.1
      iSelected          = as.numeric(levels(data$Group)[data$Group]) == selected
      selected.and.above = as.numeric(levels(data$Group)[data$Group]) >= selected
      selected.angle = 360/sum(data$Freq) * data$Freq[iSelected]
      rotation = 180 - selected.angle/2
      data = dplyr::bind_rows(data[selected.and.above, ], data[!selected.and.above, ])
    }
    else {
      pull = 0
      rotation = 0
    }
    plotly::plot_ly(data, source = sprintf("source_specificity_pie_%s", variable.name)) %>%
      plotly::add_pie(hole=0.3, values = ~Freq, marker = list(colors = ~I(colors), line = list(color = "#ffffff", width = 1)), labels = ~labels, hoverinfo = "label", sort = FALSE, direction = "clockwise", text = ~paste(Group), rotation = rotation, pull = c(pull,rep(0,7))) %>%
      layout(title = NULL,  showlegend = F,
        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
      ) %>%
      plotly::config(showLink = FALSE, displayModeBar = FALSE)
  }
  output$specificity_plot_age = plotly::renderPlotly({
    get_specificity_pie(t(getCovAge()), "AGE", input$specifiicity_plot_selector, input$specifiicity_plot_geneFinder, force_color = c("#f0f0f0","#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5"))
  })
  output$specificity_plot_sex = plotly::renderPlotly({
    get_specificity_pie(t(getCovSex()), "SEX", input$specifiicity_plot_selector, input$specifiicity_plot_geneFinder, force_color = c("#f0f0f0","#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f"))
  })
  output$specificity_plot_age_heatmap = plotly::renderPlotly({
    data = t(getCovAge())
    if (input$specifiicity_plot_geneFinder != "or, select gene") {gene.group.index = colnames(data) == input$specifiicity_plot_geneFinder}
    else {gene.group.index = apply(data > 0,2,sum) == as.integer(input$specifiicity_plot_selector)}
    color = "#6baed6"#c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf")[as.integer(input$specifiicity_plot_selector) + 1]
    data = data[, gene.group.index]
    data[data > 0] = 1
    # count = 0
    # for (row in rownames(data)) {count = count + 1; data[row, ] = data[row, ] * count * 1000 + 1}
    # colorscale = c(c((0    )/7,"#f0f0f0"), c((0+0.4)/7,"#f0f0f0"),
    #                c((1-0.4)/7,"#999999"), c((1+0.4)/7,"#999999"),
    #                c((2-0.4)/7,"#377EB8" ), c((2+0.4)/7,"#377EB8" ),
    #                c((3-0.4)/7,"#eec900"  ), c((3+0.4)/7,"#eec900"  ),
    #                c((4-0.4)/7,"#4DAF4A"  ), c((4+0.4)/7,"#4DAF4A"  ),
    #                c((5-0.4)/7,"#984EA3" ), c((5+0.4)/7,"#984EA3" ),
    #                c((6-0.4)/7,"#E41A1C"  ), c((6+0.4)/7,"#E41A1C"  ),
    #                c((7-0.4)/7,"#FF7F00"  ), c((7    )/7,"#FF7F00"  ))
    # colorscale = c(c(0,"rgb(240,240,240)"), c(1, "rgb(255,127,0)"))
    cf = as.factor(data)
    if (!is.null(dim(data))) {
      plotly::plot_ly(x = colnames(data), y = rownames(data), z = data, source = "source_specificity_age", colors = c("#ffffff",color), opacity = 0.9, type = "heatmap", showscale=FALSE) %>%
        plotly::layout(xaxis = list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = TRUE, showgrid = TRUE, gridcolor = "black", categoryorder = "trace"), yaxis = list(categoryorder = "trace", showgrid = TRUE, gridcolor = "black"), showlegend = FALSE, margin = list(b=50,l=100,r=0,t=0)) %>%
        plotly::config(displayModeBar = FALSE, showLink = FALSE)
    }
    else {
      plotly::plot_ly(x = names(data), y = data, type = "bar", color = I("#6baed6"), height = 100) %>%
        plotly::layout(xaxis = list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = TRUE, showgrid = FALSE, categoryorder = "trace"), 
                       yaxis = list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)) %>%
        plotly::config(displayModeBar = FALSE, showLink = FALSE)
    }
  })
  output$specificity_plot_sex_heatmap = plotly::renderPlotly({
    data = t(getCovSex())
    if (input$specifiicity_plot_geneFinder != "or, select gene") {gene.group.index = colnames(data) == input$specifiicity_plot_geneFinder}
    else {gene.group.index = apply(data > 0,2,sum) == as.integer(input$specifiicity_plot_selector)}
    color = "#fc8d59"#c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf")[as.integer(input$specifiicity_plot_selector) + 1]
    data = data[, gene.group.index]
    data[data > 0] = 1
    if (!is.null(dim(data))) {
      plotly::plot_ly(x = colnames(data), y = rownames(data), z = data, source = "source_specificity_sex", colors = c("#ffffff", color), opacity = 0.9, type = "heatmap", showscale=FALSE) %>%
        plotly::layout(xaxis = list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = TRUE, showgrid = TRUE, gridcolor = "black", categoryorder = "trace"), yaxis = list(categoryorder = "trace", showgrid = TRUE, gridcolor = "black"), showlegend = FALSE, margin = list(b=50,l=100,r=0,t=0)) %>%
        plotly::config(displayModeBar = FALSE, showLink = FALSE)
    }
    else {
      plotly::plot_ly(x = names(data), y = data, type = "bar", color = I("#fc8d59"), height = 100) %>%
        plotly::layout(xaxis = list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = TRUE, showgrid = FALSE, categoryorder = "trace"), 
                       yaxis = list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)) %>%
        plotly::config(displayModeBar = FALSE, showLink = FALSE)
    }
  })
  output$specificity_plot_venn = plotly::renderPlotly({
    if (input$specifiicity_plot_selector == "Select gene group") {return (plotly::plot_ly() %>% plotly::add_text(x=1,y=1,text="Please select a gene group (above) to explore.") %>% plotly::config(displayModeBar = FALSE, showLink = FALSE))}
    #color = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf")[as.integer(input$specifiicity_plot_selector) + 1]
    color.sex = "#fc8d59"
    color.age = "#6baed6"
    SEX = t(getCovSex())
    gene.group.index = apply(SEX > 0,2,sum) == as.integer(input$specifiicity_plot_selector)
    SEX = colnames(SEX[, gene.group.index])
    AGE = t(getCovAge())
    gene.group.index = apply(AGE > 0,2,sum) == as.integer(input$specifiicity_plot_selector)
    AGE = colnames(AGE[, gene.group.index])
    SEX.and.AGE = length(intersect(SEX, AGE))
    SEX.not.AGE = length(setdiff(SEX, AGE))
    AGE.not.SEX = length(setdiff(AGE, SEX))
    plotly::plot_ly(x = c(1.5,2.5,3.5), y = c(2,2,2), mode = 'text', text = c(AGE.not.SEX, SEX.and.AGE, SEX.not.AGE), textposition = 'bottom middle', hoverinfo = 'text') %>%
      plotly::layout(title = "Age \u2229 Sex",
                     margin = list(b=200,l=0,r=0,t=50),
                     xaxis = list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
                     yaxis = list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
                     shapes = list(
                       list(type = 'circle', 
                            xref = 'x', x0 = 1, x1 = 3,
                            yref = 'y', y0 = 1, y2 = 3,
                            fillcolor = color.age, line = list(color = 'rgb(0,0,0)'),
                            opacity = 0.4),
                       list(type = 'circle', 
                            xref = 'x', x0 = 2, x1 = 4,
                            yref = 'y', y0 = 1, y2 = 3,
                            fillcolor = color.sex, line = list(color = 'rgb(0,0,0)'),
                            opacity = 0.4)
                      )) %>%
      plotly::config(displayModeBar = FALSE, showLink = FALSE)
  })
  output$specificity_plot_venn_list = DT::renderDataTable({
    #- Get data
    data = getCovNames()
    stimulusName = data$stimulusName
    AGE = data$age
    SEX = data$sex
    COMBINED = data$combined
    #- Make table
    data.frame(Gene = COMBINED, Age = COMBINED %in% AGE, Sex = COMBINED %in% SEX, "Age and Sex" = (COMBINED %in% AGE) * (COMBINED %in% SEX) == 1)}, selection = "single", filter="bottom", rownames = FALSE, options = list(dom = 'tp', rowCallback = DT::JS('function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {if (String(aData[1]) == "true") $("td:eq(1)", nRow).css("background-color", "#6baed6"); if (String(aData[2]) == "true") $("td:eq(2)", nRow).css("background-color", "#fc8d59");}')))
  getFullEffectTable = shiny::reactive({
    Sex = getCovSex()                ; Age = getCovAge()
    Sex = Sex[order(rownames(Sex)), ]; Age = Age[order(rownames(Age)), ]
    Sex[Sex > 0] = 1                 ; Age[Age > 0] = 2
    Sex.Age = Sex + Age
    Sex.Age[Sex.Age == 0] = "-"
    Sex.Age[Sex.Age == "1"] = "Sex"
    Sex.Age[Sex.Age == "2"] = "Age"
    Sex.Age[Sex.Age == "3"] = "Age & Sex"
    return (as.data.frame(Sex.Age))
  })
  output$specificity_plot_full_table = DT::renderDataTable({
    getFullEffectTable()
  }, selection = "single", filter = "top")
  output$specifiicity_plot_geneFinderUI = shiny::renderUI({
    load("data/gene.data.Rdata")
    all.genes = unique(levels(gene.data$probe))
    shiny::selectInput("specifiicity_plot_geneFinder", NULL, choices = c("or, select gene", all.genes))
  })
  observeEvent(input$specifiicity_plot_geneFinder, {
    if (input$specifiicity_plot_geneFinder == "or, select gene") {return ()}
    shiny::updateSelectInput(session, "specifiicity_plot_selector", selected = "Select gene group")
  })
  observeEvent(input$specifiicity_plot_selector, {
    if (input$specifiicity_plot_selector == "Select gene group") {return ()}
    shiny::updateSelectInput(session, "specifiicity_plot_geneFinder", selected = "or, select gene")
  })
  observeEvent(input$specificity_plot_venn_list_rows_selected, {
    #- Get data
    data = getCovNames()
    COMBINED = data$combined
    gene.selected = as.character(COMBINED[input$specificity_plot_venn_list_rows_selected])
    rv$clickedGene = gene.selected
    triggerGeneModal()
  })
  observeEvent(input$specifiicity_plot_fullTable, {
    shiny::showModal(shiny::modalDialog(title = "Age & Sex Effects Table", size = 'l', easyClose = TRUE, DT::dataTableOutput("specificity_plot_full_table")))
  })
  observeEvent(rv$heatmapClicked, {
    if (is.null(rv$heatmapClicked)) {return ()}
    nv$heatmapClicked = rv$heatmapClicked
    genes = getCovNames()$combined
    shiny::showModal(shiny::modalDialog(title = shiny::HTML(sprintf("%s Genes influenced by Age/Sex under <b>%s</b> stimulation.", length(genes), shiny::isolate(nv$heatmapClicked))), size = 'm', easyClose = TRUE,
      DT::dataTableOutput("specificity_plot_venn_list")
    ))
    rv$heatmapClicked = NULL
  }) #- Render the modal window for the specifity venn table
  ########### >> eQTL Distribution ----
  ###########
  output$eQTL_Distribution = plotly::renderPlotly({
    source("./scripts/theme_publication.R")
    load("./data/eQTL.GeneEffect.df.stim.Rdata")
    load("./data/eQTL.GeneEffect.df.null.Rdata")
    #- Subsample stim data
    #if (input$overview_geneSet == "All genes") {
    if (TRUE) {
      #df.stim = dplyr::sample_frac(df.stim, input$overview_subsample/100, replace = FALSE) REMOVE subsampling
      #df.null = dplyr::sample_frac(df.null, input$overview_subsample/100, replace = FALSE) REMOVE subsampling
      df.stim = dplyr::sample_frac(df.stim, 100/100, replace = FALSE)
      df.null = dplyr::sample_frac(df.null, 100/100, replace = FALSE)
    } else {
      df.stim = subset(df.stim, lab == input$overview_geneSet)
      df.null = subset(df.stim, lab == input$overview_geneSet)
    }
    if (dim(df.stim)[1] == 0) {
      return (plotly::plot_ly(x=1,y=1) %>% plotly::add_text(text="No eQTL data to be<br />shown for this gene") %>% 
                      plotly::config(displayModeBar = FALSE, showLink = FALSE)
             )
    }
    #- Get chromosome data
    chr.l = getChromLength()
    chr.labels = chr.l$chr
    chr.labels[c(19,21)] = ""
    total.length = sum(as.numeric(chr.l$length))
    autosom.length = sum(as.numeric(chr.l$length[1:22]))
    #- Create labeler for latin names of certain stimuli
    to_latin_name = as_labeller(c('null'="NS",
                                  'candida'="C.albicans",
                                  'ecoli'="E.coli",
                                  'saureus'="S.aureus",
                                  'flu'="Influenza",
                                  'bcg'="BCG",
                                  'seb'="SEB"))
    #- Create facet plot
    p2 = ggplot2::ggplot(df.stim,aes(x=x2,y=y2))
    p2 = p2 + ggplot2::scale_shape_identity() + scale_size_area(name = "", max_size = 5)
    #circles for null eQTLs
    p2 = p2 + ggplot2::geom_point(data=df.null, aes(fill=col, size=p, text=sprintf("NULL: %s (%s)", lab, cis.trans)),alpha=0.95,pch=21)
    #circles for stimulated eQTLs
    p2 = p2 + ggplot2::geom_point(aes(fill=col,size=p, text=sprintf("Gene name: %s<br />Effect type: %s-eQTL<br />-10*log(p-value): %s", lab, cis.trans, round(p,2))),alpha=0.6,pch=21)
    p2 = p2 + xlab("\n\neQTL position (chr)\n\n") + ylab("Gene position (chr)")
    p2 = p2 + facet_wrap(~stim, labeller=to_latin_name)
    p2 = p2 + scale_y_continuous(limits=c(0,total.length), breaks = chr.l$chr.start.pos, labels = chr.labels) 
    p2 = p2 + scale_x_continuous(limits=c(0,autosom.length), breaks = chr.l$chr.start.pos[1:22], labels = chr.labels[1:22])
    p2 = p2 + theme_Publication(grid.major.color = "grey88",axis.tick.lwd = 0.2) +
      scale_fill_Publication(guide=F)+scale_colour_Publication(guide=F) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=rel(0.9)),
            axis.text.y = element_text(size=rel(0.9)),
            strip.background = element_rect(fill = alpha("#377EB8",0.4)),
            strip.text = element_text(face="bold",size=rel(1.2)))
    p2 = plotly::ggplotly(p2, source = "source_eQTLDistribution", tooltip = c("text")) %>% 
         plotly::config(displayModeBar = FALSE, showLink = FALSE, hoverinfo = "text") %>%
         plotly::layout(showlegend = FALSE, margin = list(b=40,l=50,r=0,t=25), xaxis = list(title=""), yaxis = list(title=""))
    (p2)
  })
  ########### >> User events ----
  ###########
  observeEvent(input$overview_geneSet, {
    if (input$overview_geneSet == "All genes") {return ()}
    rv$clickedGene = input$overview_geneSet
    triggerGeneModal()
    shiny::updateSelectInput(session, "overview_geneSet", selected = "All genes")
  })
  observeEvent(plotly::event_data("plotly_click", source = "source_eQTLDistribution"), {
    event.data = plotly::event_data("plotly_click", source = "source_eQTLDistribution")
    load("./data/eQTL.GeneEffect.df.stim.Rdata")
    load("./data/eQTL.GeneEffect.df.null.Rdata")
    result = df.stim[df.stim$x2 == event.data$x & df.stim$y2 == event.data$y,]
    if (dim(result)[1] == 0) {return ()}
    rv$clickedGene = unique(result$lab)
    triggerGeneModal()
  })
  observeEvent(plotly::event_data("plotly_click", source = "source_specificity_age" ), {
    event.data = plotly::event_data("plotly_click", source = "source_specificity_age")
    rv$heatmapClicked = event.data$y
  })
  observeEvent(plotly::event_data("plotly_click", source = "source_specificity_sex" ), {
    event.data = plotly::event_data("plotly_click", source = "source_specificity_sex")
    rv$heatmapClicked = event.data$y
  })
  observeEvent(rv$clickedGeneCount, {
    if (is.null(rv$clickedGene)) {return ()}
    no.sig.associations = dim(dplyr::filter(getSigAssociations(), Gene == rv$clickedGene))[1]
    plurality = ifelse(no.sig.associations == 1, "", "s")
    shiny::showModal(shiny::modalDialog(title = sprintf("Gene symbol: %s", rv$clickedGene), size = 'l', easyClose = TRUE,
      shiny::tabsetPanel(
        shiny::tabPanel("Summary",
          shiny::fluidRow(
            shiny::column(width = 8,
              shiny::fluidRow(width = 12, shiny::htmlOutput("gene_card_presence")),
              shiny::br(),
              shiny::fluidRow(
                shiny::column(width=5, shiny::htmlOutput("gene_card_general"), shiny::a("View gene in GeneCards", href=sprintf("http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", rv$clickedGene), target="_blank")),
                shiny::column(width=6, shiny::htmlOutput("gene_card_summary")),
                shiny::column(width=1, shiny::p("", style = "margin-bottom: 290px;"))
              )
            ),
            shiny::column(width=4, 
                          shiny::fluidRow(shiny::column(width = 6, shiny::h3("Summary of Effects", style="text-align: center;")),
                                          shiny::column(width = 6, shiny::selectizeInput("gene_card_effects_select", "", choices = c("Default view", "Genetic", "Intrinsic", "Cellular", "Gen/Env/Cel")))),
                          plotly::plotlyOutput("gene_card_effects", height = "300px"))
          ),
          shiny::hr(),
          shiny::h3("Text mining data"), DT::dataTableOutput("gene_card_DT")
        ),
        shiny::tabPanel("Analysis of effects",
          shiny::helpText("Use the tabs below to select which aspect of the study you wish to explore further."),
          shiny::tabsetPanel(
            shiny::tabPanel("Effect of Age", plotly::plotlyOutput('examine_gene_age')),#, shiny::hr(), shinyBS::bsCollapse(shinyBS::bsCollapsePanel(title = "Options", shiny::checkboxInput("examine_gene_age_legend", "Show legend?", value = FALSE)))),
            shiny::tabPanel("Effect of Sex", plotly::plotlyOutput('examine_gene_sex')),#, shiny::hr(), shinyBS::bsCollapse(shinyBS::bsCollapsePanel(title = "Options", shiny::checkboxInput("examine_gene_sex_legend", "Show legend?", value = FALSE)))),
            shiny::tabPanel("Effect of Influenza across Age", plotly::plotlyOutput('examine_gene_inf'))
          )
        ),
        shiny::tabPanel("Significant eQTL associations",
          shiny::helpText(sprintf("%s eQTL%s found associated to this gene", no.sig.associations, plurality)),
          DT::dataTableOutput('examine_gene_eQTL_table'),
          if (no.sig.associations > 0) {
            shiny::div(
              plotly::plotlyOutput('examine_gene_eQTL_boxplot'),
              shiny::uiOutput('examine_gene_eQTL_boxplot_title')#,
              #shiny::hr(),
              #shinyBS::bsCollapse(
              #  shinyBS::bsCollapsePanel(title = "Options",
              #    shiny::checkboxInput("examine_gene_eQTL_legend", "Show legend?", value = TRUE)
              # )
              #)
            )}
          else {}
        )
      )
    ))
  }) #- Render the modal window for gene view
  ########### >>> In app help modals ----
  ###########
  observeEvent(input$help_specificity_age,{
    shiny::showModal(shiny::modalDialog(title = "Specificity of age effects on gene expression", size = 'm', easyClose = TRUE,
      shiny::tags$p("Specificity of age effect on gene expression across 7 conditions (including the NULL condition).")
    ))
  }) #- Render the modal window for help on specificity of age effects on gene expression
  observeEvent(input$help_specificity_sex,{
    shiny::showModal(shiny::modalDialog(title = "Specificity of sex effects on gene expression", size = 'm', easyClose = TRUE,
                                        shiny::tags$p("Specificity of sex effect on gene expression across 7 conditions.")
    ))
  }) #- Render the modal window for help on specificity of sex effects on gene expression
  observeEvent(input$help_specificity_detailed,{
    shiny::showModal(shiny::modalDialog(title = "Detailed map of gene-stimuli effects", size = 'm', easyClose = TRUE,
                                        shiny::tags$p("Specificity of age (left) and sex (right) effect on gene expression across conditions.")
    ))
  }) #- Render the modal window for help on detailed mappings of specificity of effects on gene expression
  observeEvent(input$help_eQTL_distribution,{
    shiny::showModal(shiny::modalDialog(title = "Genome wide distribution of eQTLs", size = 'm', easyClose = TRUE,
                                        shiny::tags$p("Distribution of significant eQTLs. On the x-axis is shown the genomic position of the SNP, on the y-axis is shown the position of the associated gene. Axis numbers indicate chromosome names. Circles corresponding to local eQTLs are on the diagonal, while off-axis circles correspond to trans-eQTLs. On each panel, eQTLs detected in the non-stimulated state are plotted first, and are overlaid with eQTLs detected in one of the stimulated condition.")
    ))
  }) #- Render the modal window for help on specificity of age effects on gene expression
  ########### >>> Gene card data ----
  ###########
  output$gene_card_general = shiny::renderText({
    sprintf("%s (%s)<br />%s", getGeneFeature("name"), getGeneFeature("symbol"), getGeneFeature("map_location"), getGeneFeature("type_of_gene")) 
  })
  output$gene_card_summary = shiny::renderText({
    data = getGeneCard()
    sprintf("%s", data$summary) 
  })
  output$gene_card_effects = plotly::renderPlotly({
    variance.explained = subset(getVarianceExplained(), GENE == rv$clickedGene)
    variance.explained$Genetic     = variance.explained$CIS + variance.explained$TRANS
    variance.explained$Environment = variance.explained$AGE + variance.explained$SEX
    variance.explained$Cellular    = variance.explained$TOTAL.CD45.
    if      (input$gene_card_effects_select == "Default view" ) {view = 1; overall = dplyr::select(variance.explained, STIM, TOTAL.CD45., AGE, SEX, CIS, TRANS); colnames(overall) = c("Stimulus", "Total.CD45", "Age", "Sex", "CIS", "TRANS")}
    else if (input$gene_card_effects_select == "Genetic"      ) {view = 2; overall = dplyr::select(variance.explained, STIM, CIS, TRANS); colnames(overall) = c("Stimulus", "CIS", "TRANS")}
    else if (input$gene_card_effects_select == "Intrinsic"    ) {view = 3; overall = dplyr::select(variance.explained, STIM, AGE, SEX); colnames(overall) = c("Stimulus", "Age", "Sex")}
    else if (input$gene_card_effects_select == "Cellular"     ) {view = 4; overall = dplyr::select(variance.explained, STIM, CD4., CD19., CD8B., NK, MONO, CD4.CD8., CD8B.CD4., NEUTROPHILS); colnames(overall) = c("Stimulus", "CD4", "CD19", "CD8B", "NK", "MONO", "CD4.CD8", "CD8B.CD4", "NEUTROPHILS")}
    else                                                        {view = 5; overall = dplyr::select(variance.explained, STIM, Genetic, Environment, Cellular); colnames(overall) = c("Stimulus", "Genetic", "Environment", "Cellular")}
    overall$Stimulus = factor(as.character(overall$Stimulus), levels = c("NS", "E.coli", "BCG", "S.aureus", "SEB", "C.albicans", "IVA"))
    overall[is.na(overall)] = 0
    p = plotly::plot_ly(overall, x = ~Stimulus, y = rep(0, 7), type = 'bar', text = "", hoverinfo = "text", showlegend = FALSE)
    if (view %in% c(1  )) {p = plotly::add_trace(p, y = ~Total.CD45 , name = 'CD45+'      , color = I("#636363"   ), text = ~round(Total.CD45 ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(1,3)) {p = plotly::add_trace(p, y = ~Age        , name = 'Age'        , color = I("skyblue1"  ), text = ~round(Age        ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(1,3)) {p = plotly::add_trace(p, y = ~Sex        , name = 'Sex'        , color = I("lightcoral"), text = ~round(Sex        ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(1,2)) {p = plotly::add_trace(p, y = ~CIS        , name = 'CIS'        , color = I("#e5c494"   ), text = ~round(CIS        ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(1,2)) {p = plotly::add_trace(p, y = ~TRANS      , name = 'TRANS'      , color = I("chocolate4"), text = ~round(TRANS      ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(4  )) {p = plotly::add_trace(p, y = ~CD4        , name = 'CD4+'       , color = I("gray1"     ), text = ~round(CD4        ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(4  )) {p = plotly::add_trace(p, y = ~CD19       , name = 'CD19+'      , color = I("gray15"    ), text = ~round(CD19       ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(4  )) {p = plotly::add_trace(p, y = ~CD8B       , name = 'CD8B+'      , color = I("gray24"    ), text = ~round(CD8B       ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(4  )) {p = plotly::add_trace(p, y = ~NK         , name = 'NK'         , color = I("gray36"    ), text = ~round(NK         ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(4  )) {p = plotly::add_trace(p, y = ~MONO       , name = 'MONO'       , color = I("gray48"    ), text = ~round(MONO       ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(4  )) {p = plotly::add_trace(p, y = ~CD4.CD8    , name = 'CD4-CD8-'   , color = I("gray65"    ), text = ~round(CD4.CD8    ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(4  )) {p = plotly::add_trace(p, y = ~CD8B.CD4   , name = 'CD8B+CD4+'  , color = I("gray72"    ), text = ~round(CD8B.CD4   ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(4  )) {p = plotly::add_trace(p, y = ~NEUTROPHILS, name = 'NEUTROPHILS', color = I("gray84"    ), text = ~round(NEUTROPHILS,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(5  )) {p = plotly::add_trace(p, y = ~Genetic    , name = 'Genetic'    , color = I("#e5c494"   ), text = ~round(Genetic    ,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(5  )) {p = plotly::add_trace(p, y = ~Environment, name = 'Environment', color = I("skyblue1"  ), text = ~round(Environment,2), hoverinfo = "text", showlegend = TRUE)}
    if (view %in% c(5  )) {p = plotly::add_trace(p, y = ~Cellular   , name = 'Cellular'   , color = I("#636363"   ), text = ~round(Cellular   ,2), hoverinfo = "text", showlegend = TRUE)}
    p = plotly::layout(p, xaxis = list(title="", showgrid = TRUE), yaxis = list(title = "expression<br /> variance explained", showgrid = TRUE), barmode = 'stack')
    plotly::config(p, displayModeBar = FALSE, showLink = FALSE)
  })
  output$gene_card_presence= shiny::renderUI({
    if (is.null(rv$clickedGene)) {return ()}
    load("data/gene.presence.Rdata")
    load("data/gene.data.Rdata")
    source("./scripts/custom.task.list.R")
    source("./scripts/custom.dashboardMenu.R")
    gene = rv$clickedGene
    data = subset(gene.presence, X == gene)
    min.data = min(dplyr::select(subset(gene.presence, X == gene), Null:IVA))
    if (min.data > 0.9) {badgeStatus = "success"}
    else if (min.data > 0.5) {badgeStatus = "warning"}
    else {badgeStatus = "danger"}
    probe = subset(gene.data, probe == gene)$nr.snps.in.probes
    rsIDs = stringr::str_extract_all(subset(gene.data, probe == gene)$snps.in.probes, "rs[0-9]+")[[1]]
    if (probe > 0) {
      notifications = lapply(rsIDs, function(rsID) {shinydashboard::notificationItem(rsID, icon = shiny::icon("exclamation-triangle"))})
      probe = customTaskList(type = "notifications", badgeStatus = "warning", custom_message = sprintf("%s probes contain %s SNPs with MAF > 0.05 in the 1000 Genome project; probe binding affinities may be affected.", gene, probe), .list = notifications)
    }
    else {
      probe = customTaskList(type = "notifications", badgeStatus = "success", custom_message = sprintf("%s probes are likely not affected by genetic variation.", gene))
    }
    return (
      custom.dashboardMenu(title = "Data and description",
        customTaskList(type = "tasks", badgeStatus = badgeStatus, custom_message = sprintf("Percentage of donors per stimulus where %s transcripts were detected:", gene),
          shinydashboard::taskItem(value = round(as.numeric(data$Null)*100), color = getPresenceQColor(data$Null), "Null"),
          shinydashboard::taskItem(value = round(as.numeric(data$E.coli)*100), color = getPresenceQColor(data$E.coli), "E.coli"),
          shinydashboard::taskItem(value = round(as.numeric(data$BCG)*100), color = getPresenceQColor(data$BCG), "BCG"),
          shinydashboard::taskItem(value = round(as.numeric(data$S.aureus)*100), color = getPresenceQColor(data$S.aureus), "S.aureus"),
          shinydashboard::taskItem(value = round(as.numeric(data$SEB)*100), color = getPresenceQColor(data$SEB), "SEB"),
          shinydashboard::taskItem(value = round(as.numeric(data$C.albicans)*100), color = getPresenceQColor(data$C.albicans), "C.albicans"),
          shinydashboard::taskItem(value = round(as.numeric(data$IVA)*100), color = getPresenceQColor(data$IVA), "Influenza")
        ),
        probe
    ))
  })
  output$gene_card_DT = DT::renderDataTable({
    data = getGeneCard()
    if ("generif" %in% names(data)) {
      output = as.data.frame(t(matrix(unlist(data$generif), nrow=length(unlist(data$generif[1])))))
      colnames(output) = c("Pubmed ID", "Context")
    }
    else {
      output = data.frame(Context = c("None"), "Pubmed ID" = c("None"))
    }
    output$`Pubmed ID` = sapply(output$`Pubmed ID`, function(pmid) {sprintf("<a href='https://www.ncbi.nlm.nih.gov/pubmed/%s/' target='_blank'>%s</a>",pmid,pmid)})
    return (output)
  }, options = list(dom = "tp"), filter = "top", escape = FALSE)
  ########### >>> Gene specific plots ----
  ###########
  output$examine_gene_inf = plotly::renderPlotly({
    source("./scripts/theme_publication.R")
    to_latin_name = getLabeler()
    df.subset = subset(getGeneResponseData(), gene == rv$clickedGene & stim == "flu")
    g1 = ggplot(df.subset, aes(x=factor(age_cat),y=expr))
    #g1 = g1 + geom_boxplot(fill="#FF7F00",outlier.shape=20,outlier.size=0.05,size=0.05)
    g1 = g1 + geom_boxplot(fill=alpha("#FF7F00", 0.8), outlier.shape=20, outlier.size=0.05, size=0.05)
    g1 = g1 + xlab("") + ylab("") + theme_Publication(base_size = 14) + theme(strip.background = element_blank()) + theme(axis.title.y = element_text(size = 8, angle = 90)) + theme(axis.title.x = element_text(size = 8, angle = 0))
    g1 = ggplotly(g1) %>%
      plotly::config(displayModeBar = FALSE, showLink = FALSE) %>%
      plotly::layout(showlegend = FALSE, margin = list(b=55,l=25,r=0,t=15), xaxis = list(title="Age group"), yaxis = list(title=""))
    return (g1)
  })
  output$examine_gene_age = plotly::renderPlotly({
    source("./scripts/theme_publication.R")
    to_latin_name = getLabeler()
    df.subset = subset(getGeneResponseData(), gene == rv$clickedGene & stim != "flu")
    #df.subset = dplyr::sample_frac(df.subset, input$examine_gene_age_subsample/100, replace = FALSE) REMOVE subsample
    #df.subset = dplyr::sample_frac(df.subset, 100/100, replace = FALSE)
    df.subset$stim = factor(c("NS","E.coli","BCG","S.aureus","SEB","C.albicans")[factor(df.subset$stim, levels = c("null","ecoli","bcg","saureus","seb","candida"))], levels = c("NS","E.coli","BCG","S.aureus","SEB","C.albicans"))
    g1 = ggplot(df.subset,aes(x=age,y=expr,color=stim))  
    g1 = g1 + geom_point(size=0.7,alpha=0.1,pch=16)
    g1 = g1 + facet_grid(~stim)#,labeller=to_latin_name)
    g1 = g1 + stat_smooth(method = "lm",size=0.5)
    g1 = g1 + scale_x_continuous(breaks = c(20,40,60))
    g1 = g1 + xlab("") + ylab("") + theme_Publication(base_size = 18) + scale_colour_Publication(guide=F) + theme(strip.background = element_blank()) + theme(strip.text.x = element_text(size = 12)) + theme(axis.title.y = element_text(size = 8, angle = 90)) + theme(axis.title.x = element_text(size = 8, angle = 0))
    g1 = ggplotly(g1) %>%
         plotly::config(displayModeBar = FALSE, showLink = FALSE) %>%
         plotly::layout(showlegend = FALSE, margin = list(b=30,l=30,r=0,t=20), xaxis = list(title=""), yaxis = list(title=""))
    return (g1)
  })
  output$examine_gene_sex = plotly::renderPlotly({
    source("./scripts/theme_publication.R")
    to_latin_name = getLabeler()
    df.subset = subset(getGeneResponseData(),gene == rv$clickedGene)
    df.subset$stim = factor(c("NS","E.coli","BCG","S.aureus","SEB","C.albicans","Influenza")[factor(df.subset$stim, levels = c("null","ecoli","bcg","saureus","seb","candida","flu"))], levels = c("NS","E.coli","BCG","S.aureus","SEB","C.albicans","Influenza"))
    g1 = ggplot(df.subset,aes(x=factor(sex),y=expr,fill=stim))
    g1 = g1 + geom_boxplot(outlier.shape=20,outlier.size=0.1, size = 0.1, color = "black")
    g1 = g1 + facet_grid(~stim)
    g1 = g1 + scale_x_discrete(labels = c("F","M"))
    g1 = g1 + xlab("") + ylab("") + theme_Publication(base_size = 18) + scale_fill_Publication(guide=F) + theme(strip.background = element_blank()) + theme(strip.text.x = element_text(size = 12)) + theme(axis.title.y = element_text(size = 8, angle = 90)) + theme(axis.title.x = element_text(size = 8, angle = 0))
    g1 = ggplotly(g1) %>%
      plotly::config(displayModeBar = FALSE, showLink = FALSE) %>%
      plotly::layout(showlegend = FALSE, margin = list(b=30,l=30,r=0,t=20), xaxis = list(title=""), yaxis = list(title=""))
    return (g1)
  })
  output$examine_gene_eQTL_boxplot = plotly::renderPlotly({
    if (is.null(input$examine_gene_eQTL_table_rows_selected)) {return (plotly::plot_ly() %>% plotly::add_text(x=1,y=1,text="Please select a<br />significant eQTL using<br />the table above.") %>% plotly::config(displayModeBar = FALSE, showLink = FALSE))}
    source("./scripts/EQTLboxplot.R")
    table = getSigAssociationsByGene()
    i = input$examine_gene_eQTL_table_rows_selected
    gene_response_data = getGeneResponseData(); gene_response_data = subset(gene_response_data, gene == rv$clickedGene)
    gene_snp_pvalue    = dplyr::filter(dplyr::select(table, SNP, Stimulus, Pvalue), SNP == table$SNP[i])
    p = getEQTLboxplot(table$Stimulus[i], table$SNP[i], table$Gene[i], gene_response_data, table$MAF[i], gene_snp_pvalue)
    p = plotly::plotly_build(p)
    #p$data <- lapply(p$data, FUN = function(x){
    #  print (x$marker)
    #  x$marker = list(size = 1)
    #  return(x)
    #})
    p %>% plotly::config(displayModeBar = FALSE, showLink = FALSE) %>%
      plotly::layout(showlegend = FALSE, margin = list(b=75,l=55,r=0,t=100),
                     yaxis = list(title = sprintf("%s Expression", table$Gene[i])))
  })
  output$examine_gene_eQTL_boxplot_title = shiny::renderUI({
    if (is.null(input$examine_gene_eQTL_table_rows_selected)) {return ("")}
    table = getSigAssociationsByGene()
    i = input$examine_gene_eQTL_table_rows_selected
    shiny::h3(sprintf("%s (MAF = %s)", table$SNP[i], table$MAF[i]), style="text-align: center;")
  })
  ########### >>> Gene specific tables ----
  ###########
  output$examine_gene_eQTL_table = DT::renderDataTable({
    subset(getSigAssociations(), Gene == rv$clickedGene)
  }, selection = "single", filter = "top", options = list(dom = "tp"))
  ########### >>> Search table ----
  ###########
  output$custom_ui_geneSelection = shiny::renderUI({
    load("data/gene.data.Rdata")
    all.genes = unique(levels(gene.data$probe))
    shiny::selectInput("overview_geneSet", NULL, choices = c("All genes", all.genes))
  })
  output$table_search = DT::renderDataTable({
    getSigAssociations()
  }, selection = "single", filter = "top",
     options = list(dom = "tp", columnDefs = list(list(width = '75px', targets = "_all"))))
  observeEvent(input$table_search_rows_selected, {
    data = getSigAssociations()
    gene.selected = as.character(data[input$table_search_rows_selected, ]$Gene)
    rv$clickedGene = gene.selected
    triggerGeneModal()
  })
  observeEvent(input$specificity_plot_full_table_rows_selected, {
    data = getFullEffectTable()
    gene.selected = rownames(data)[input$specificity_plot_full_table_rows_selected]
    rv$clickedGene = gene.selected
    triggerGeneModal()
  })
})