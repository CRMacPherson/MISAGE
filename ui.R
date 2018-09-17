########### > HEADER ----
###########
shinydashboard_header = dashboardHeader(title = shiny::p(tags$a(href='http://www.milieuinterieur.fr/en',
                                                       tags$img(src='./images/mi_logo.png', height="50", style="display: block; margin-left: auto; margin-right: auto;")), "Milieu Interieur"))
########### > SIDEBAR ----
###########
shinydashboard_sidebar = dashboardSidebar(sidebarMenu(
  shinydashboard::menuItem("Home"    , tabName = "tab_home"   , icon = shiny::icon("home"), selected = TRUE),
  shinydashboard::menuItem("Explore"    , tabName = "tab_overview"   , icon = shiny::icon("th"), startExpanded = TRUE,
    shinydashboard::menuSubItem("Age/Sex Effects"  , tabName = "tab_overview_gene", icon = shiny::icon("pie-chart")),
    shinydashboard::menuSubItem("eQTL Analysis"    , tabName = "tab_overview_eqtl", icon = shiny::icon("bar-chart")),
    shinydashboard::menuSubItem("Find Associations", tabName = "tab_overview_search", icon = shiny::icon("search"))
  ),
  #menuItem("Lists & Reports", tabName = "tab_lists", icon = icon("info-circle")),
  #menuItem("Download", tabName = "tab_download", icon = icon("download")),
  shinydashboard::menuItem("About"   , tabName = "tab_about"  , icon = icon("info-circle"    ))
))
########### > BODY ----
########### >> TAB: HOME ----
###########
tab.home     = tabItem(tabName = "tab_home",
  shiny::tags$head(shiny::tags$style(".content-wrapper, .right-side {background-color: white;}")),
  shiny::div(
    #shiny::div(shiny::img(src='./images/mi_logo.png', height = 100), style="display:inline-block; vertical-align:top"),
    shiny::div(shiny::br(), shiny::h1("The human immune response and the impact of age, sex, and genetics"), style="display:inline-block; vertical-align:center"),
  style = "text-align:center"),
  shiny::hr(),
  shiny::tags$p("Welcome! This work is brought to you by the Milieu Interieur Consortium and is published by", shiny::tags$a("Piasecka and Duffy et al., 2017.", href="http://www.pnas.org/content/115/3/E488", target="_blank"), align="center"),
  shiny::tags$p(shiny::tags$b("ABSTRACT"), "Immune responses are highly variable between individuals and populations, with this variance mediated by many factors including age, sex, and genetics. To dissect how each of these factors contributes to differential immune responses we recruited 1,000 healthy donors as part of the Milieu Intereiur cohort that were equally stratified by age (20-70 years old) and sex (50:50). Whole blood from each donor was stimulated in a standardized approach with Escherichia coli, BCG, Staphylococcus aureus, SEB, Candida albicans and Influenza virus and the transcriptomic response was analyzed by Nanostring gene expression arrays. We show that age affects gene expression in a stimuli specific manner, while the effect of sex is more common across conditions. Both age and sex contribute to the expression of many immune-related genes but at a low level, while in contrast, the genetic effects are much stronger but relevant to fewer genes. Genetic analysis identified multiple eQTLs, both local and trans, for all stimuli induced responses, including a novel, strong trans-eQTL located near the CR1 gene after C. albicans stimulation. Finally, we employed a structural equation model to describe how age and sex effects are mediated by different immune cell populations. These results lay the foundations for understanding how age, sex, and genetics may impact differentially induced immune responses, which may be considered for future clinical strategies proposing stratified or precision medicine approaches."),
  shiny::hr(),
  shiny::fluidRow( 
    shiny::column(width = 6, shiny::h1("Study outline", style="text-align: center;"), shiny::img(src="./images/process.png", style="display: block; margin-left: auto; margin-right: auto; max-height:550px; max-width:100%;")),
    shiny::column(width = 6, 
      shiny::h1("Transcriptomic response", style="text-align: center;"), 
      plotly::plotlyOutput("pca", width = "100%"),#, height="500"),
      shiny::helpText("Click-Drag to zoom, Shift-Click-Drag to rotate; Use 'Options' below to adjust subsampling."),
      shiny::hr(),
      shinyBS::bsCollapse(
        shinyBS::bsCollapsePanel(title = "Options",
          shiny::radioButtons("pca_grouping", "Color group", choices = c("Stimulus", "Sex", "Age"), inline = TRUE, width = "100%"),
          shiny::sliderInput("pca_subsample", "Subsample (%)", min = 1, max = 100, value = 10, step = 1)
        )
      )
    )
  )
)
########### >> TAB: DATA (SUPER-MENU) ----
########### >>> TAB: DATA-GENE ----
###########
tab.overview.gene = tabItem(tabName = "tab_overview_gene",
  shiny::h1("Specificity of age and sex effects on gene expression"),
  shiny::helpText("The effect of age is stimulus-specific with ~40% of genes presenting an age-dependent expression pattern in only one or two conditions. Conversely, sex effects on gene expression are more frequently shared across stimuli with ~45% of genes differing in their expression between sexes in at least six or seven conditions."),
  shiny::hr(),
  shiny::fluidRow(
    shiny::column(width = 5,
      shiny::h3(shiny::actionLink("help_specificity_age","",icon = shiny::icon("question-circle-o")), "Age effects on gene expression", style="text-align: center"),
      plotly::plotlyOutput("specificity_plot_age", width = "100%")
    ),
    shiny::column(width=2,
      shiny::div(style="padding-top: 50px",
        shiny::helpText("Selecting a gene group below will open interactive heatmaps showing which genes are affected by age/sex in which stimuli. Group-1 genes are affected in only 1 stimulus, group-2 in 2 stimuli, and group-7 in all stimuli."),
        shiny::selectInput("specifiicity_plot_selector", "", c("Select gene group","1","2","3","4","5","6","7")),
        shiny::uiOutput("specifiicity_plot_geneFinderUI"),
        shiny::actionButton("specifiicity_plot_fullTable", "View all effects", icon = shiny::icon("table"))
      )
    ),
    shiny::column(width=5,
      shiny::h3(shiny::actionLink("help_specificity_sex","",icon = shiny::icon("question-circle-o")), "Sex effects on gene expression", style="text-align: center"),
      plotly::plotlyOutput("specificity_plot_sex", width = "100%")
    )#,
    #shiny::h3(shiny::actionLink("help_specificity_heatmap","",icon = shiny::icon("question-circle-o")), "Some responsive feature to clicks on the left")
  ),
  shiny::conditionalPanel("input.specifiicity_plot_selector != 'Select gene group' || input.specifiicity_plot_geneFinder != 'or, select gene'", {
    shiny::fluidRow(
      shiny::column(width = 5, plotly::plotlyOutput("specificity_plot_age_heatmap", width = "100%")),
      shiny::column(width = 2, 
        shiny::conditionalPanel("input.specifiicity_plot_selector != 'Select gene group'", 
          shiny::h3(shiny::actionLink("help_specificity_detailed","",icon = shiny::icon("question-circle-o")), "Detailed map of gene-stimuli effects", style="text-align: center"),
          shiny::helpText("Clicking on any point in either heatmap will bring up a list of genes affected by age/sex in a particular stimulus."),
          plotly::plotlyOutput("specificity_plot_venn")
        )
      ),
      shiny::column(width = 5, plotly::plotlyOutput("specificity_plot_sex_heatmap", width = "100%"))
    )
  })
)
########### >>> TAB: DATA-eQTL ----
###########
tab.overview.eqtl = tabItem(tabName = "tab_overview_eqtl",
  shiny::h1("Overview of genotypic effects on gene expression"),
  shiny::helpText("135 Genes were found with local Cis eQTLs in the unstimulated state and around 140 genes presenting an eQTL after stimulation. On the other hand, trans-regulation was noticeably stronger after stimulation, particularly with E. coli, BCG, C. albicans and SEB."),
  shiny::hr(),
  shiny::fluidRow(
    shiny::column(width = 1),
    shiny::column(width = 10, shiny::h3(shiny::actionLink("help_eQTL_distribution","",icon = shiny::icon("question-circle-o")), "Genome wide distribution of eQTLs"), plotly::plotlyOutput("eQTL_Distribution", height = 700)),
    shiny::column(width = 1)
  )
)
########### >>> TAB: DATA-SEARCH ----
###########
tab.overview.search = tabItem(tabName = "tab_overview_search",
  shiny::helpText("Use the tools below to explore our dataset. The 'Quick gene view' will allow you to navigate to a particular gene of your choosing whereas the 'Table search' will give you more fine grained control in filtering through our genetic association results. Clicking on any gene in the 'Quick gene view' or on any row in the 'Table search' will bring up a detailed report for the respective gene. Genotypic data relevant to the gene will be found in this report."),
  shiny::h3("Quick gene view"),
  shiny::uiOutput("custom_ui_geneSelection"),
  shiny::hr(),
  shiny::h3("Table search"),
  DT::dataTableOutput("table_search")
)
########### >> TAB: DOWNLOAD ----
###########
tab.download = tabItem(tabName = "tab_download"    ,
  shiny::h1("Download study data"),
  shiny::helpText("Here, you will find access to primary and supplementary datasets. Much of this data may be explored in this app, but we welcome you to contact us for further information. When using data viewed or downloaded from this app please cite the article listed in the 'About' section of this app."),
  shiny::hr(),
  shiny::fluidRow(
    shiny::column(width = 1 , shiny::downloadButton('d1', label = "Download-1")),
    shiny::column(width = 11, "Main | Description-1...")
  ), shiny::hr(),
  shiny::fluidRow(
    shiny::column(width = 1 , shiny::downloadButton('d2', label = "Download-2")),
    shiny::column(width = 11, "Main | Description-2...")
  ), shiny::hr(),
  shiny::fluidRow(
    shiny::column(width = 1 , shiny::downloadButton('d3', label = "Download-3")),
    shiny::column(width = 11, "Main | Description-3...")
  ), shiny::hr(),
  shiny::fluidRow(
    shiny::column(width = 1 , shiny::downloadButton('d4', label = "Download-2")),
    shiny::column(width = 11, "Supp | Description-2...")
  )
)
########### >> TAB: ABOUT ----
###########
tab.about = tabItem(tabName = "tab_about",
  shiny::div(
    #shiny::div(shiny::img(src='./images/mi_logo.png', height = 100), style="display:inline-block; vertical-align:top"),
    shiny::div(shiny::br(), shiny::h1("The human immune response and the impact of age, sex, and genetics"), style="display:inline-block; vertical-align:center"),
    style = "text-align:center"),
  shiny::hr(),
  shiny::tags$p("This work is brought to you by the Milieu Interieur Consortium and is published by", shiny::tags$a("Piasecka and Duffy et al., 2017.", href="http://www.pnas.org/content/115/3/E488", target="_blank"), "with the following title, authorship, and affiliations:", align="center"),
  shiny::tags$p(shiny::tags$b("Impact of age, sex and genetics on transcriptional immune responses to bacterial, viral and fungal challenges."), "Barbara Piasecka",shiny::tags$sup("1,2*"),", Darragh Duffy",shiny::tags$sup("2,3,4*"),", Alejandra Urrutia",shiny::tags$sup("3,4,5"),", Helene Quach",shiny::tags$sup("1,6,7"),", Etienne Patin",shiny::tags$sup("1,6,7"),", Céline Posseme",shiny::tags$sup("2"),", Bruno Charbit",shiny::tags$sup("2"),", Vincent Rouilly",shiny::tags$sup("2"),", Cameron MacPherson",shiny::tags$sup(2),", Milena Hasan",shiny::tags$sup(2),", Benoit Albaud",shiny::tags$sup("8"),", David Gentien",shiny::tags$sup(8),", Jacques Fellay",shiny::tags$sup("9,10"),", Matthew L. Albert",shiny::tags$sup("2,3,4,5,#"),", Lluis Quintana-Murci",shiny::tags$sup("1,6,7,#"),", for The Milieu Intérieur Consortium",shiny::tags$sup("†"),"."),
  shiny::tags$ul(style="list-style: none;",
    shiny::tags$li("1  - Laboratory of Human Evolutionary Genetics, Department of Genomes & Genetics, Institut Pasteur, Paris, 75015, France"),
    shiny::tags$li("2  - Center for Translational Research, Institut Pasteur, Paris, 75015, France"),
    shiny::tags$li("3  - Laboratory of Dendritic Cell Immunobiology, Department of Immunology, Institut Pasteur, Paris, 75015, France"),
    shiny::tags$li("4  - INSERM U1223, Paris, 75015, France"),
    shiny::tags$li("5  - Department of Cancer Immunology, Genentech Inc, San Francisco, CA 94080, USA"),
    shiny::tags$li("6  - CNRS URA3012, Paris, 75015, France"),
    shiny::tags$li("7  - Center of Bioinformatics, Biostatistics and Integrative Biology, Institut Pasteur, Paris 75015, France"),
    shiny::tags$li("8  - Institut Curie, Centre de Recherche, Département de recherche translationnelle, Plateforme de Génomique, Paris, 75005, France"),
    shiny::tags$li("9  - School of Life Sciences, École Polytechnique Fédérale de Lausanne, Lausanne 1015, Switzerland"),
    shiny::tags$li("10 - Swiss Institute of Bioinformatics, Lausanne 1015, Switzerland"),
    shiny::tags$li("*  - Indicates co-first authors"),
    shiny::tags$li("#  - Indicates co-senior authors"),
    shiny::tags$li("†  - The Milieu Intérieur Consortium",shiny::tags$sup("¶")," is composed of the following team leaders: Laurent Abel (Hôpital Necker), Andres Alcover, Kalla Astrom (Lund University), Philippe Bousso, Pierre Bruhns, Ana Cumano, Caroline Demangel, Ludovic Deriano, James Di Santo, Françoise Dromer, Gérard Eberl, Jost Enninga, Jacques Fellay (EPFL, Lausanne) Antonio Freitas, Odile Gelpi, Ivo Gomperts-Boneca, Serge Hercberg (Université Paris 13), Olivier Lantz (Institut Curie), Claude Leclerc, Hugo Mouquet, Sandra Pellegrini, Stanislas Pol (Hôpital Côchin), Lars Rogge, Anavaj Sakuntabhai, Olivier Schwartz, Benno Schwikowski, Spencer Shorte, Vassili Soumelis (Institut Curie), Frédéric Tangy, Eric Tartour (Hôpital Européen George Pompidou), Antoine Toubert (Hôpital Saint-Louis), Marie-Noëlle Ungeheuer, Lluis Quintana-Murci",shiny::tags$sup("§"),", Matthew L. Albert",shiny::tags$sup("§")),
    shiny::tags$li("¶ - unless otherwise indicated, partners are located at Institut Pasteur, Paris"),
    shiny::tags$li("§ - co-coordinators of the Milieu Intérieur Consortium")
  ),
  shiny::hr(),
  shiny::tags$p("Data are made available from the electronic version of the publication (",shiny::a("here", href="http://www.pnas.org/content/115/3/E488", target="_blank"),")", align="center"),
  shiny::hr()
)
########### > AGGREGATE TABS ----
###########
shinydashboard_body = dashboardBody(
  tags$head(tags$style(HTML('.modal-lg {
                            width: 80%;}'))),
tabItems(
  tab.home,
  tab.overview.gene,
  tab.overview.eqtl,
  tab.overview.search,
  #tab.lists,
  #tab.download,
  tab.about
))
########### > GENERATE HTML ----
###########
dashboardPage(shinydashboard_header,
  shinydashboard_sidebar,
  shinydashboard_body,
  skin = "blue")
