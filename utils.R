### +-------------------------------------------------------------+
### |  Visualisation: Gene Sets                                   |
### +-------------------------------------------------------------+

# extract genes with abs(log2FC)
extract_regulated_sets <- function(list_of_degs, names_degs, fc_thresh=1.5) {
  # fc_thresh must be a positive number
  # list_of_degs is a list of dataframes
  # names_degs is a list of strings
  # they must have the same length
  reg_lvl <- paste0(' reg@log~2~FC≥', fc_thresh)
  breaks <- brk_manual(c(-fc_thresh, fc_thresh), left_vec = c(FALSE, TRUE))
  regulated_sets <- NULL
  for (l in 1:length(list_of_degs)) {
    degs_na <- dplyr::select(list_of_degs[[l]], c(log2FoldChange, padj, ensemblGeneID))
    # better not pass NAs to kiru
    degs <- na.omit(degs_na)
    # filter by log2FC threshold
    degs$reg <- 
      kiru(
        degs$log2FoldChange,
        breaks = breaks,
        extend = TRUE,
        labels=c("down", "non", "up")
      )
    # filter by p-val
    degs$reg <- as.character(degs$reg)
    x <- 1:nrow(degs)
    degs$reg <- ifelse(degs[x,'padj']<0.5, degs[x,'reg'], 'non')
    # recover NAs as non-regulataed
    added_nas <- degs_na[!(rownames(degs_na) %in% rownames(na.omit(degs_na))),]
    added_nas$reg <- 'non'
    degs <- rbind(degs, added_nas)
    # store reg status
    #rownames(degs) <- degs$ensemblGeneID
    degs[ , names_degs[[l]] ] <- degs$reg
    regulated_sets[[l]] <- dplyr::select(degs, c(ensemblGeneID, names_degs[[l]]) )
  }
  regulated_sets <- purrr::reduce(regulated_sets, full_join, by='ensemblGeneID')
  rownames(regulated_sets) <- regulated_sets$ensemblGeneID
  regulated_sets <- dplyr::select(regulated_sets, -ensemblGeneID)
  return(regulated_sets)
}

get_deg_logical <- function(regs, direction){
  # direction must be 'up' or 'down'
  # regs is the output of `extract_regulated_sets`
  deg_logical_set <- (regs == direction)[
    ( rowSums(regs==direction)!=0 ), 
  ]
  return(deg_logical_set)
}


### +-------------------------------------------------------------+
### |  Visualisation: MA plots                                    |
### +-------------------------------------------------------------+


# classify lists of genes as up/non/down-regulated
extract_regulated_sets2 <- function(list_of_degs, names_degs, fc_thresh=1.5, cols) {
  # fc_thresh must be a positive number
  # list_of_degs is a list of dataframes
  # names_degs is a list of strings
  # they must have the same length
  reg_lvl <- paste0(' reg@log~2~FC≥', fc_thresh)
  breaks <- brk_manual(c(-fc_thresh, fc_thresh), left_vec = c(FALSE, TRUE))
  regulated_sets <- NULL
  for (l in 1:length(list_of_degs)) {
    degs_na <- dplyr::select(list_of_degs[[l]], all_of(cols))
    # better not pass NAs to kiru
    degs <- na.omit(degs_na)
    # filter by log2FC threshold, create `reg`(ulated) col with up/down/non
    degs$reg <- 
      kiru(
        degs$log2FoldChange,
        breaks = breaks,
        extend = TRUE,
        labels=c("down", "non", "up")
      )
    # filter by p-val
    degs$reg <- as.character(degs$reg) # remove factor
    x <- 1:nrow(degs)
    degs$reg <- ifelse(degs[x,'padj']<0.05, degs[x,'reg'], 'non')
    # recover NAs as non-regulataed
    added_nas <- degs_na[!(rownames(degs_na) %in% rownames(na.omit(degs_na))),]
    added_nas$reg <- 'non'
    degs <- rbind(degs, added_nas)
    # store reg status
    degs[ , names_degs[[l]] ] <- degs$reg
    regulated_sets[[l]] <- dplyr::select(degs, c(ensemblGeneID, names_degs[[l]]) )
  }
  regulated_sets <- purrr::reduce(regulated_sets, full_join, by='ensemblGeneID')
  rownames(regulated_sets) <- regulated_sets$ensemblGeneID
  regulated_sets <- dplyr::select(regulated_sets, -ensemblGeneID)
  return(regulated_sets)
}

# custom MA plot
ggmaplot2 <- function(deg, markers, fc_thresh = 1.5,
                      md_label='*genotype^—^*', repulsion) {
  # deg is the output of DESeq2 with an added column of gene symbols
  # gene.symbols is that column
  
  make_goilist_ggmaplot2 <- function(deg, fc_thresh, markerlist) {
    goi_list <- deg %>% dplyr::filter(gene_symbol %in% markerlist &
                                        abs(log2FoldChange) > fc_thresh &
                                        padj<0.05)
    return(goi_list)
  }
  
  # standard MA plot with `ggpubr`, using the CBD1 palette from `cetcolor` 
  ma <- ggmaplot(
    deg,
    fdr = 0.05,
    fc = fc_thresh,
    genenames = deg$gene_symbols,
    size = 2,
    alpha = 0.5,
    seed = NA,
    font.label = c(16, "bold", "black"),
    label.rectangle = FALSE,
    palette = c(cet_pal(n = 3, name = "cbd1", alpha = 1)[3], # "#A89008"
                cet_pal(n = 3, name = "cbd1", alpha = 1)[1], # "#3A90FE"
                "#AAAAAA"),
    top = 0,
    main = NULL,
    xlab = "log~2~(mean expression)",
    ylab = "log~2~(fold change)",
    ggtheme = theme_linedraw(),
    legend = 'top'
  )
  
  annotations <- data.frame(
    xpos = floor(min(ggplot_build(ma)$layout$panel_params[[1]]$x.range)),
    ypos = ceiling(max(abs(ggplot_build(ma)$layout$panel_params[[1]]$y.range))),
    annotateText = md_label,
    hjustvar = 0,
    vjustvar = 1)
  
  goi_list <- make_goilist_ggmaplot2(deg, fc_thresh, unlist(markers$xmarkers))
  
  ma <- ma +
     # mark special genes if upreg
     geom_text_repel(data = deg %>% filter(gene_symbol %in% goi_list$gene_symbol &
                                             log2FoldChange > 0),
                     aes(x = log2(baseMean),
                         y = log2FoldChange,
                         label = gene_symbol,
                         segment.square  = FALSE,
                         segment.inflect = TRUE),
                     color         = '#6E5C03', # ==hue, <ligthness than than "#A89008"
                     fontface = 'bold',
                     segment.alpha = 0.8,
                     segment.linetype = 3,
                     hjust = 0,
                     # positioning
                     box.padding   = repulsion$box.padding,
                     point.padding = repulsion$point.padding,
                     nudge_x       = repulsion$nudge_x.up,
                     nudge_y       = repulsion$nudge_y.up,
                     force         = repulsion$force,
                     force_pull    = repulsion$force_pull,
                     max.overlaps  = Inf,
                     xlim          = repulsion$xlims.up,    # NA repels from edges
                     ylim          = repulsion$ylims.up,
                     seed          = repulsion$seed.up) +
     # mark special genes if downreg
     geom_text_repel(data = deg %>% filter(gene_symbol %in% goi_list$gene_symbol &
                                             log2FoldChange < 0),
                     aes(x = log2(baseMean),
                         y = log2FoldChange,
                         label = gene_symbol,
                         segment.square  = FALSE,
                         segment.inflect = TRUE),
                     color         = '#0565AF', # ==hue, <ligthness than "#3A90FE"
                     fontface = 'bold',
                     segment.alpha = 0.8,
                     segment.linetype = 3,
                     hjust = 0,
                     # positioning
                     box.padding   = repulsion$box.padding,
                     point.padding = repulsion$point.padding,
                     nudge_x       = repulsion$nudge_x.dn,
                     nudge_y       = repulsion$nudge_y.dn,
                     force         = repulsion$force,
                     force_pull    = repulsion$force_pull,
                     max.overlaps  = Inf,
                     xlim          = repulsion$xlims.dn,             # NA repels from edges
                     ylim          = repulsion$ylims.dn,
                     seed          = repulsion$seed.dn) +
     # genotype label
     geom_richtext(
       data = annotations,
       aes(x = xpos,y = ypos,
           hjust = hjustvar, vjust = vjustvar,
           label = annotateText,
           fontface = 'bold',
           size = 9),
       fill = NA,
       label.color = NA,
       label.padding = grid::unit(rep(0, 4), "pt"),
       label.margin = grid::unit(rep(0, 4), "pt")
     ) +
     # legend
     guides(size = 'none') +
     # apply markdown formatting, etc
     theme(# markdown
       axis.title.x = element_markdown(size = 14, face = 'bold'),
       axis.title.y = element_markdown(size = 14, face = 'bold'),
       axis.text = element_text(size = 10),
       # grid and panel
       panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
       panel.grid = element_blank(),
       # legend
       legend.margin = margin(b = -10),
       legend.key.size = unit(1.2, 'cm'),
       legend.text = element_text(size = 12, face = 'bold')
       )
  return(ma)
}

# custom MA plot with labels coloured by cell type
ggmaplot3 <- function(deg, markers, fc_thresh = 1.5,
                      md_label='*genotype^—^*', repulsion) {
  # deg is the output of DESeq2 with an added column of gene symbols
  # gene.symbols is that column
  
  goi_list  <- tidyr::unnest_longer(markers, xmarkers) %>%
    dplyr::rename(gene_symbol = xmarkers)
  
  # standard MA plot with `ggpubr`, using the CBD1 palette from `cetcolor` 
  ma <- ggmaplot(
    deg,
    fdr = 0.05,
    fc = fc_thresh,
    genenames = deg$gene_symbols,
    size = 2,
    alpha = 0.5,
    seed = NA,
    font.label = c(16, "bold", "black"),
    label.rectangle = FALSE,
    palette = c(cet_pal(n = 3, name = "cbd1", alpha = 1)[3], # "#A89008"
                cet_pal(n = 3, name = "cbd1", alpha = 1)[1], # "#3A90FE"
                "#AAAAAA"),
    top = 0,
    main = NULL,
    xlab = "log~2~(mean expression)",
    ylab = "log~2~(fold change)",
    ggtheme = theme_linedraw(),
    legend = 'top'
  )
  
  annotations <- data.frame(
    xpos = floor(min(ggplot_build(ma)$layout$panel_params[[1]]$x.range)),
    ypos = ceiling(max(abs(ggplot_build(ma)$layout$panel_params[[1]]$y.range))),
    annotateText = md_label,
    hjustvar = 0,
    vjustvar = 1)
  
  genlabsup <- deg %>%
    filter(gene_symbol %in% goi_list$gene_symbol &
             log2FoldChange > fc_thresh &
             padj < 0.05) %>%
    dplyr::select(c(gene_symbol, baseMean, log2FoldChange)) %>%
    dplyr::left_join(dplyr::select(goi_list, !'celltype'), by = 'gene_symbol')
  
  genlabsdn <- deg %>%
    filter(gene_symbol %in% goi_list$gene_symbol &
             log2FoldChange < -fc_thresh &
             padj < 0.05) %>%
    dplyr::select(c(gene_symbol, baseMean, log2FoldChange)) %>%
    dplyr::left_join(dplyr::select(goi_list, !'celltype'), by = 'gene_symbol')
  
  ma <- ma +
    # mark special genes if upreg
    geom_text_repel(data = genlabsup,
                    aes(x = log2(baseMean),
                        y = log2FoldChange,
                        label = gene_symbol,
                        segment.square  = FALSE,
                        segment.inflect = TRUE),
                    colour = genlabsup$cellcolour,
                    segment.color = '#6E5C03', # ==hue, <ligthness than than "#A89008"
                    fontface = 'bold',
                    segment.alpha = 0.8,
                    segment.linetype = 3,
                    hjust = 0,
                    # positioning
                    box.padding   = repulsion$box.padding,
                    point.padding = repulsion$point.padding,
                    nudge_x       = repulsion$nudge_x.up,
                    nudge_y       = repulsion$nudge_y.up,
                    force         = repulsion$force,
                    force_pull    = repulsion$force_pull,
                    max.overlaps  = Inf,
                    xlim          = repulsion$xlims.up,    # NA repels from edges
                    ylim          = repulsion$ylims.up,
                    seed = repulsion$seed.up) +
    # mark special genes if downreg
    geom_text_repel(data = genlabsdn,
                    aes(x = log2(baseMean),
                        y = log2FoldChange,
                        label = gene_symbol,
                        segment.square  = FALSE,
                        segment.inflect = TRUE),
                    colour = genlabsdn$cellcolour,
                    segment.color = '#0565AF', # ==hue, <ligthness than "#3A90FE"
                    fontface = 'bold',
                    segment.alpha = 0.8,
                    segment.linetype = 3,
                    hjust = 0,
                    # positioning
                    box.padding   = repulsion$box.padding,
                    point.padding = repulsion$point.padding,
                    nudge_x       = repulsion$nudge_x.dn,
                    nudge_y       = repulsion$nudge_y.dn,
                    force         = repulsion$force,
                    force_pull    = repulsion$force_pull,
                    max.overlaps  = Inf,
                    xlim          = repulsion$xlims.dn,             # NA repels from edges
                    ylim          = repulsion$ylims.dn,
                    seed = repulsion$seed.dn) +
    # genotype label
    geom_richtext(
      data = annotations,
      aes(x = xpos,y = ypos,
          hjust = hjustvar, vjust = vjustvar,
          label = annotateText,
          fontface = 'bold',
          size = 9),
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt"),
      label.margin = grid::unit(rep(0, 4), "pt")
    ) +
    # legend
    guides(size = 'none') +
    # apply markdown formatting, etc
    theme(# markdown
      axis.title.x = element_markdown(size = 14, face = 'bold'),
      axis.title.y = element_markdown(size = 14, face = 'bold'),
      axis.text = element_text(size = 10),
      # grid and panel
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      # legend
      legend.margin = margin(b = -10),
      legend.key.size = unit(1.2, 'cm'),
      legend.text = element_text(size = 12, face = 'bold')
    )
  return(ma)
}

# gene-based layered heatmaps

# horizontal
layer.heatmaph <- function(genehm.df, cluster=FALSE, arr=NULL) {
  logFC.df <- genehm.df %>% dplyr::select(-p.adjust)
  padj.df <- genehm.df %>% dplyr::select(-log2FoldChange)
  if (cluster) {
    vectors <- genehm.df %>%
      dplyr::select(gene_symbol, log2FoldChange, condition, cellcolour) %>%
      pivot_wider(names_from = condition, values_from = log2FoldChange)
    clust <- hclust(dist(as.matrix(vectors[3:6])))
    logFC.df <- logFC.df %>%
      mutate(gene_symbol = factor(gene_symbol, levels=vectors$gene_symbol[clust$order]))
  } else if (!cluster & !is.null(arr)) {
    logFC.df <- logFC.df %>% arrange(!!as.name(arr), gene_symbol, condition) 
    padj.df <- padj.df %>% arrange(!!as.name(arr), gene_symbol, condition)
    clust <- data.frame(order = 1:length(unique(logFC.df$gene_symbol)))
  }
  xlab.colours <- logFC.df$cellcolour[seq(1, nrow(logFC.df), 4)][clust$order]
  p <- ggplot(logFC.df, aes(x=gene_symbol, y=condition)) +
    # plot statistic (log2fc)
    geom_tile(aes(fill=log2FoldChange), width=1) +
    scale_fill_gradient2(low = cet_pal(3, name='cbd1')[1],
                         mid = cet_pal(3, name='cbd1')[2],
                         high = cet_pal(3, name='cbd1')[3],
                         midpoint = 0) +
    # plot p-value
    # statistically significant with colour and shape `*`
    geom_point(data=subset(padj.df, p.adjust<0.05),
               aes(x=gene_symbol, y=condition, colour=-log10(p.adjust)),
               size=3, shape=8, stroke=1.5, alpha=1) +
    scale_colour_gradient(low = cet_pal(3, name='d2')[2],
                          high = cet_pal(3, name='d2')[3]) +
    coord_equal() +
    theme_bw()
  if (cluster) {
    p <- p +
      scale_x_discrete(labels = levels(logFC.df$gene_symbol),
                       position = "top",
                       expand = expansion(mult = 0, add = 0))
  } else {
    p <- p +
      scale_x_discrete(position = "top",
                       expand = expansion(mult = 0, add = 0))
  }
  p <- p +
    scale_y_discrete(expand = expansion(mult = 0, add = 0)) +
    labs(fill = 'log~2~FC',
         colour = '-log~10~(_p.adj_)<br><span style = "font-size:8pt;">_p.adj_<0.05</span>') +
    theme(axis.text.x = element_markdown(angle=25, hjust=0,
                                         face='bold', size=10.5,
                                         colour = xlab.colours),
          axis.text.y = element_markdown(hjust=1, face='bold', size=12),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_markdown(hjust=0.5, vjust=0.75),
          legend.direction = 'horizontal',
          legend.position = 'bottom',
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 1))
  return(p)
}

# vertical
layer.heatmapv <- function(genehm.df, cluster=FALSE, arr=NULL) {
  logFC.df <- genehm.df %>% dplyr::select(-p.adjust)
  padj.df <- genehm.df %>% dplyr::select(-log2FoldChange)
  if (cluster) {
    vectors <- genehm.df %>%
      dplyr::select(gene_symbol, log2FoldChange, condition, cellcolour) %>%
      pivot_wider(names_from = condition, values_from = log2FoldChange)
    clust <- hclust(dist(as.matrix(vectors[3:6])))
    logFC.df <- logFC.df %>%
      mutate(gene_symbol = factor(gene_symbol, levels=vectors$gene_symbol[clust$order]))
  } else if (!cluster & !is.null(arr)) {
    logFC.df <- logFC.df %>% arrange(!!as.name(arr), gene_symbol, condition) 
    padj.df <- padj.df %>% arrange(!!as.name(arr), gene_symbol, condition)
    clust <- data.frame(order = 1:length(unique(logFC.df$gene_symbol)))
  }
  ylab.colours <- logFC.df$cellcolour[seq(1, nrow(logFC.df), 4)][clust$order]
  p <- ggplot(logFC.df, aes(x=condition, y=gene_symbol)) +
    # plot statistic (log2fc)
    geom_tile(aes(fill=log2FoldChange), width=1) +
    scale_fill_gradient2(low = cet_pal(3, name='cbd1')[1],
                         mid = cet_pal(3, name='cbd1')[2],
                         high = cet_pal(3, name='cbd1')[3],
                         midpoint = 0) +
    # plot p-value
    # statistically significant with colour and shape `*`
    geom_point(data=subset(padj.df, p.adjust<0.05),
               aes(x=condition, y=gene_symbol, colour=-log10(p.adjust)),
               size=3, shape=8, stroke=1.5, alpha=1) +
    scale_colour_gradient(low = cet_pal(3, name='d2')[2],
                          high = cet_pal(3, name='d2')[3]) +
    coord_equal() +
    theme_bw()
  if (cluster) {
    p <- p +
      scale_y_discrete(labels = levels(logFC.df$gene_symbol),
                       expand = expansion(mult = 0, add = 0),
                       position = 'right')
  } else {
    p <- p +
      scale_y_discrete(expand = expansion(mult = 0, add = 0),
                       position = 'right')
  }
  p <- p +
    scale_x_discrete(expand = expansion(mult = 0, add = 0),
                     position = 'bottom') +
    labs(fill = 'log~2~FC',
         colour = '-log~10~(_p.adj_)<br><span style = "font-size:8pt;">_p.adj_<0.05</span>') +
    theme(axis.text.y = element_markdown(angle=0, hjust=0,
                                         face='bold', size=12,
                                         colour = ylab.colours),
          axis.text.x = element_markdown(hjust=0.5, face='bold', size=20,
                                         angle=25),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_markdown(hjust=0.5, vjust=0.75),
          legend.direction = 'vertical',
          legend.position = 'left',
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 1))
  return(p)
}


### +-------------------------------------------------------------+
### |  Gene Set Enrichment Analysis                               |
### +-------------------------------------------------------------+


# To get a differentially expressed gene set for Over-Representation Analysis
make_degset <- function(deg, up, fc_thresh) {
  converter <- ifelse(up,1,-1)
  degset <- deg %>%
    filter(padj<0.05) %>%
    filter(abs(log2FoldChange) > fc_thresh) %>%
    filter(log2FoldChange*converter > 0)
  return(degset)
}

# To get a ranked gene list for Gene Set Enrichment Analysis
make_degrank <- function(deg, mode='log2fc', key='gene_symbol') {
  if (!(key %in% c('gene_symbol', 'ensemblGeneID'))) {
    stop("The `key` parameter needs to be one of 'gene_symbol', 'ensemblGeneID'")
  }
  # mode: whether to use padj or log2FC
  if (mode=='log2fc') {
    deg_sorted <- deg %>% arrange(desc(log2FoldChange))
    degrank <- deg_sorted$log2FoldChange
    names(degrank) <- deg_sorted[,key]
  } else if (mode=='padj') {
    deg$p.rank <- sign(deg$log2FoldChange) * -log10(deg$padj)
    deg_sorted <- deg %>% arrange(desc(p.rank))
    degrank <- deg_sorted$p.rank
    names(degrank) <- deg_sorted[,key]
  } else {
    stop('the `mode` parameter can only be either `log2fc` or `padj`')
  }
  return(degrank)
}

# To import the gene sets in the gmx files
import_from_gmx <- function(gmxfile) {
  # read a GMX file and turn it into a df input appropriate for `clusterProfiler`
  df <- read.csv(gmxfile, header=TRUE, sep='\t')
  # remove 'na's in 1st row
  df <- df[ !df[1]=='na', ]
  # pivot
  df <- df %>% pivot_longer(cols = everything(),
                            names_to = 'term',
                            values_to = 'gene')
  # remove rows with no genes
  df <- df %>% subset(gene != '')
  df <- df %>% mutate(gene = str_replace(gene, "FBGN", "FBgn"),
                      term = str_replace(term, ".only", "-only"))
  df <- df %>% mutate(term = str_replace(term, "\\.", " "))
  df <- df %>% mutate(term = str_replace(term, " genes", ""))
  return(df)
}

# To prepare GSEA results for plotting as layered heatmaps.
gseCP_summarise <- function(gmx, gse_list, conditions, sets.as.factors, cluster=FALSE, nsig.out=FALSE) {
  # purr::map to convert the gse_list from S4 objects to their @result slots 
  gseCP_list <- map( gse_list, \(x) dplyr::select(x@result, NES, p.adjust, ID) )
  # name them to associate conditions with the data
  # add condition as an extra column
  gseCP_list <- lapply( 1:length(gseCP_list), \(x) cbind(gseCP_list[[x]],
                                                         condition = conditions[[x]]) )
  df <- bind_rows(gseCP_list)
  df$condition <- factor(df$condition, levels = conditions)
  df$ID <- factor(df$ID, levels = sets.as.factors)
  # apply filtering
  if (nsig.out) {
      # get the IDs for which there is at least one condition with significant enrichment
      sign_lgl <- lapply( lapply(unique(df$ID),
                                 \(x) filter(df, ID==x)),
                          \(x) any(x$p.adjust<0.05) )
      df_bycondition <- lapply(unique(df$ID),
                               \(x) filter(df, ID==x) )[unlist(sign_lgl)]
      df <- bind_rows(df_bycondition)
      # to avoid passing on non-filtered terms
      df$ID <- as.character(df$ID)
      levels(df$ID) <- factor(unique(df$ID))
  }
  # apply clustering
  if (cluster & length(sets.as.factors)>2) { # `hclust` must have n >= 2 objects to cluster
    vectors <- df %>%
      dplyr::select(ID, NES, condition) %>%
      pivot_wider(names_from = condition, values_from = NES) %>%
      mutate_all(replace_na, 100)
    clust <- hclust(dist(as.matrix(vectors[2:length(gseCP_list)])))
    df$ID <- factor(df$ID, levels=vectors$ID[clust$order])
  } else if (cluster & length(sets.as.factors < 3)) {
    warning("`hclust` must have n >= 2 objects to cluster. NES columns will not be clustered.")
  }
  return(df)
}

# create a heatmap with NES in colour and coloured point as padj.
layer.heatmap <- function(layerhm.df, subt) {
  NES.df <- layerhm.df %>% dplyr::select(-p.adjust)
  padj.df <- layerhm.df %>% dplyr::select(-NES)
  p <- ggplot(NES.df, aes(x=ID, y=condition)) +
    # plot statistic (NES)
    geom_tile(aes(fill=NES), width=1) +
    scale_fill_gradient2(low = cet_pal(3, name='cbd1')[1],
                         mid = cet_pal(3, name='cbd1')[2],
                         high = cet_pal(3, name='cbd1')[3],
                         midpoint = 0) +
    # plot p-value
    # statistically significant with colour and shape `*`
    geom_point(data=subset(padj.df, p.adjust<0.05),
               aes(x=ID, y=condition, colour=-log10(p.adjust)),
               size=3, shape=8, stroke=1.5, alpha=1) +
    #size=5, shape=23, alpha=1) +
    # non-significant in dark gray and shape `x`
    # geom_point(data=subset(padj.df, p.adjust>0.05),
    #            aes(x=ID, y=condition),
    #            size=3, shape=4, stroke=2, alpha=1, colour='gray50') +
    scale_colour_gradient(low = cet_pal(3, name='d2')[2],
                          high = cet_pal(3, name='d2')[3]) +
    coord_equal() +
    theme_bw() +
    ggtitle("Normalised Enrichment Scores",
            subtitle = subt) +
    scale_x_discrete(labels=toupper(levels(NES.df$ID)),
                     position = "top",
                     expand = expansion(mult = 0, add = 0)) +
    scale_y_discrete(expand = expansion(mult = 0, add = 0)) +
    labs(fill = 'NES',
         colour = '-log~10~(_p.adj_)<br><span style = "font-size:8pt;">_p.adj_<0.05</span>') +
    theme(plot.subtitle = element_markdown(),
          axis.text.x = element_markdown(angle=25, hjust=0,
                                         face='bold', size=10.5),
          axis.text.y = element_markdown(hjust=1, face='bold', size=12),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_markdown(hjust=0.5, vjust=0.75),
          legend.direction = 'horizontal',
          legend.position = 'bottom',
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 1))
  return(p)
}

# For GLAD sub.(sub.)group terms
refine_glad_by <- function(dataset, col.name) {
  # refined by col.name
  subgrouping <- dataset %>%
    # select terms with subgroups
    filter(!is.na(.data[[col.name]])) %>%
    # select the term and subgroup cols
    dplyr::select(term, {{col.name}}) %>%
    # get unique pairs of terms/subgroup
    distinct()
  # re-make term/gene cols as per subgroups now
  gmx <- lapply(
    unique(subgrouping$term),
    \(x) dataset %>%
      filter( !is.na(.data[[col.name]]) & term=={{x}} ) %>%
      dplyr::select({{col.name}}, FBgn) %>%
      rename(gene = FBgn, term = {{col.name}})
  )
  names(gmx) <- unique(subgrouping$term)
  # remove terms that only have one subgroup
  gmx <- gmx[ sapply(gmx, \(x) length(unique(x$term)))>1 ]
  return(gmx)
}

subglad_gsea <- function(deg, gmx, perms=1000) {
  g <- lapply(
    1:length(gmx), \(x) do.call(
      GSEA,
      c(list(geneList=make_degrank(deg, mode='log2fc', key='ensemblGeneID'),
             nPermSimple = perms,
             TERM2GENE = gmx[[x]]),
        GSEAparams)
    )
  )
  names(g) <- names(gmx)
  return(g)
}


### +-------------------------------------------------------------+
### |  DNA motif enrichment                                       |
### +-------------------------------------------------------------+


class.swarm <- function(df, # dataframe
                        categories, examples, values, labels, # data variables
                        dotparams, palette, texts, emp.ratio=38/64) { # plot variables
  
  # prepare for ggplot evaluation
  categories_s <- ensym(categories)
  examples_s   <- ensym(examples)
  values_s     <- ensym(values)
  labels_s     <- ensym(labels)

  # unravel parameter bundles
  binwidth <- dotparams$binwidth
  dotsize  <- dotparams$dotsize
  alpha    <- dotparams$alpha
  stroke   <- dotparams$stroke
  x_label      <- texts$x_label
  target_class <- texts$target_class
  obj.regex    <- texts$tfs.regex
  obj.var      <- texts$tfs.var
  obj.sel      <- texts$tfs.select
  
  # get the beeswarm
  p <- ggplot(df, aes(x=!!categories_s, y=!!values_s, fill=!!labels_s, colour=!!labels_s)) + 
    geom_dotplot(binaxis     ='y',
                 stackdir    ='center',
                 stackgroups = TRUE,
                 method      = 'histodot',
                 binwidth    = binwidth,
                 dotsize     = dotsize,
                 alpha       = alpha,
                 stroke      = stroke) +
    scale_fill_manual(name   = '',
                      breaks = levels( df[, {{labels}}] ),
                      values = palette,
                      labels = levels( df[, {{labels}}] ),
                      guide  = guide_legend(nrow = 1)) +
    scale_colour_manual(name   = '',
                        breaks = levels( df[, {{labels}}] ),
                        values = palette,
                        labels = levels( df[, {{labels}}] ),
                        guide  = guide_legend(nrow = 1)) +
    xlab(x_label) +
    theme(axis.text.x     = element_markdown(face = 'bold', size = 10),
          aspect.ratio    = 0.5,
          legend.position = 'bottom')

  # prepare labels for specific genes
  # using https://stackoverflow.com/questions/44991607
  for (j in 1:length(obj.regex)) {
    df <- df %>%
      mutate('{obj.var[[j]]}' := if_else(
        (str_detect( as.character( df[, examples] ), obj.regex[[j]] ) & !!labels_s==target_class),
        obj.var[[j]], NA) )
  }
  df <- df %>%
    unite('obj.labs', all_of(obj.var), sep = ', ', na.rm = TRUE, remove=FALSE)
  
  # get the lines of ggrepel text where we want:
  built <- ggplot_build(p)
  point.pos <- built$data[[1]]
  size <- dev.size(units = 'px')
  extent <- with(built$layout$panel_params[[1]], abs(c(diff(x.range), diff(y.range))))
  bw <- point.pos$binwidth[[1]]
  # por la cuenta de la vieja, al final
  xtext <- point.pos$x + point.pos$stackpos * bw * (size[2] / size[1]) * (extent[1] / extent[2]) * (emp.ratio)
  ytext <- point.pos$y
  
  # add labels
  p <- p +
    # it would be helpful to allow user control of more parameters here
    geom_text_repel(
      # reindex - point.pos is organised by X values, then by FILL colour
      # which equates to ordering obj.lab as X, then labels, then values
      aes(label = arrange(df, {{categories}}, {{labels}}, {{values}})[, obj.sel[[1]] ] ,
          x = xtext, y = ytext,
          colour = arrange(df, {{categories}}, {{labels}}, {{values}})[, {{labels}}] ),
      size = 3,
      # general
      max.overlaps = Inf,
      # position
      direction = 'both',
      nudge_x = -0.25,
      nudge_y = 1.5,
      # segment
      min.segment.length = 0,
      segment.color = 'gray30',
      segment.alpha = 0.75,
      segment.size = 0.2) +
    
    # repeat, bc I don't know how to loop in ggplot2
    geom_text_repel(
      aes(label = arrange(df, {{categories}}, {{labels}}, {{values}})[, obj.sel[[2]] ] ,
          x = xtext, y = ytext,
          colour = arrange(df, {{categories}}, {{labels}}, {{values}})[, {{labels}}] ),
      size = 3,
      max.overlaps = Inf,
      direction = 'both',
      nudge_x = -0.25,
      nudge_y = 1.5,
      min.segment.length = 0,
      segment.color = 'gray30',
      segment.alpha = 0.75,
      segment.size = 0.2) +
  guides(fill = guide_legend(nrow=1), colour = guide_legend(nrow=1))
  return(p)
}


### +-------------------------------------------------------------+
### |  Cleaner downloads                                          |
### +-------------------------------------------------------------+


is.valid.path <- function(path) {
  require(checkmate)
  valid <- c(checkPathForOutput(path)==TRUE,
             file.exists(path),
             dir.exists(path))
  if(any(valid)) return(TRUE)
  else return(FALSE)
}

databringr <- function(from, to) {
  
  # check 'from'
  pattern <- "(https?|http?|ftp)://[^ /$.?#].[^\\s]*"
  if (!stringr::str_detect(from, pattern)) { stop(
    '`from` must be character string containing a valid URL')}
  # check 'to'
  fsep <- .Platform$file.sep
  pathels <- str_split(to, fsep)[[1]]
  topath <- paste0(pathels[1:length(pathels)-1], collapse = fsep)
  if(!dir.exists(topath) & is.valid.path(topath)) {
    dir.create(topath)
  }
  if(is.valid.path(to)) {
  m <- Sys.which(c('wget', 'curl', 'libcurl', 'internal', 'wininet'))
  for(method in names(m[nchar(m)>1])) {
    if(!file.exists(to)) {
      download.file(url = from, destfile = to,
                    method = method, quiet = TRUE)
      }
    }
  } else { stop(
    '`to` must be character string containing a valid path in your system')}
}
