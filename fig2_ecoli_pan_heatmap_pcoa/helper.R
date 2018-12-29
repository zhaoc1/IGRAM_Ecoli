library(qiimer)
library(ape)
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(kylemisc)
library(tidyr)
library(broom)
library(pheatmap)

library(lme4)
library(lmerTest)

library(subfunc)
library(forcats)
library(lubridate)
library(magrittr)
library(pander)

ggplot_contam_taxon <- function(taxon_toPlot, color_by = Proportion, nrow = 2){
  props_df %>%
    mutate(TaxonLabel = sub(" ", "\n", Taxon)) %>%
    filter(Taxon %in% taxon_toPlot) %>%
    ggplot(aes(x = library_concentration_ng_ul, y = Proportion, color = get(color_by))) +
    geom_point() + 
    geom_smooth(method='lm') +
    scale_y_log10() +
    facet_wrap( ~ TaxonLabel, nrow = nrow) +
    scale_fill_brewer(palette = "Set1", name=color_by) +
    ggtitle("Potential Contaminated OTUs")
}

logit <- function (p) log(p / (1 - p))

heatmap_clustered <- function(genus_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.1,option=1){
  color = saturated_rainbow(101)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  heatmap_props <- genus_props[,heatmap_s$SampleID]
  
  # filter: samples are filtered to retain OTUs with a relative abundance > 1% in at least one sample
  if (option == 1)
    rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) #
  else
    rows_to_keep <- apply(heatmap_props,1,max) >= thre #
  heatmap_props <- heatmap_props[rows_to_keep,]
  
  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps
  
  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, filename = fname, fontsize_col = 8, fontsize_row = 8,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, fontsize_col = 8, fontsize_row = 8, cellheight = 8, cellwidth = 8)
}

heatmap_grouped <- function(genus_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.1,option=1){
  
  color = saturated_rainbow(101)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  heatmap_props <- genus_props[,heatmap_s$SampleID]
  
  if (option == 1)
    rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) #
  else
    rows_to_keep <- apply(heatmap_props,1,max) >= thre #
  heatmap_props <- heatmap_props[rows_to_keep,]
  
  ## group the SampleIDs
  heatmap_s %<>% arrange_(.dots=grps)
  heatmap_props <- heatmap_props[, heatmap_s$SampleID]
  
  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps
  
  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, filename = fname, fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
}
heatmap_clustered_annot <- function(genus_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.1,option=1, annot=mycolor){
  color = saturated_rainbow(101)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  heatmap_props <- genus_props[,heatmap_s$SampleID]
  
  # filter: samples are filtered to retain OTUs with a relative abundance > 1% in at least one sample
  if (option == 1)
    rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) #
  else
    rows_to_keep <- apply(heatmap_props,1,max) >= thre #
  heatmap_props <- heatmap_props[rows_to_keep,]
  
  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps
  
  ## change the color of annotation
  mycolors <- annot
  names(mycolors) <- levels(heatmap_s$StudyGroup)
  mycolors <- list(StudyGroup = mycolors)
  
  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, annotation_colors=mycolors, color = color, breaks = breaks, filename = fname, fontsize_col = 8, fontsize_row = 8,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, annotation_colors=mycolors, color = color, breaks = breaks, fontsize_col = 8, fontsize_row = 8, cellheight = 8, cellwidth = 8)
}

taxonomic_colors <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", 
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
taxonomic_colors2 <- taxonomic_colors[1:6]

taxonomic_barchart <- function(cts, props,s_toTest,fill=taxonomic_colors, grp=c("RunNumber", "study_group"), option=1, nrow=2){
  ## Helps us to know whether authetic difference or not!
  
  ## subset the matrix first
  props <- props[,s_toTest$SampleID]
  cts <- cts[,s_toTest$SampleID]
  
  if (option == 1)
    row_vals <- apply(props, 1, max) # top dominant
  else 
    row_vals <- apply(props, 1, mean) # top abundant
  
  num_taxa <- length(fill)
  top_rows <- order(row_vals, decreasing = TRUE)[seq(num_taxa)]
  top_taxa <- rownames(props)[top_rows]
  
  top_taxa_labels <- ifelse(rownames(props) %in% top_taxa, rownames(props), "Other")
  top_props <- rowsum(props, top_taxa_labels)
  
  ## very nice: wide to long
  ## sort the SampleID by NumReads
  props_df <- s_toTest %>%
    merge(t(top_props), by.x = "SampleID", by.y = "row.names") %>%
    arrange(`Read.Counts`) %>% #arrange_(.dots = top_taxa[1])
    mutate(SampleID = factor(SampleID, levels = .$SampleID)) %>%
    melt(colnames(s_toTest), rownames(top_props), value.name = "Proportion", variable.name = "Taxon") %>%
    mutate(Taxon = as.factor(Taxon))
  
  ## reorder: generic function, default FUN=mean...
  props_df$Taxon <- reorder(props_df$Taxon, props_df$Proportion)
  props_df %<>% mutate(Taxon = fct_relevel(Taxon, "Other"))
  props_df %<>% mutate(Taxon = fct_rev(props_df$Taxon))
  
  props_df %<>%
    group_by(SampleID) %>% 
    arrange(Taxon) %>%
    ungroup()
  
  grp1 <- grp[1]
  grp2 <- grp[2]
  props_df %>% 
    filter(SampleID %in% s_toTest$SampleID) %>%
    ggplot(aes(x=SampleID, y=Proportion, fill=Taxon, width=0.9)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c(fill, "#CCCCCC")) + #, breaks = NULL, drop=TRUE
    scale_y_continuous(limits = c(0,1.01), expand=c(0,0.1)) + 
    facet_wrap(factor(get(grp1)) ~ factor(get(grp2)), scales = "free_x", nrow=nrow) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
  #theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
}

make_pcoa_plot_discrete <- function(uu, s, shape_by, color_by, gtitle="") {
  uu_pcoa <- pcoa(uu)
  uu_df <- merge(s, uu_pcoa$vectors[, 1:5], by.x="SampleID", by.y="row.names")
  uu_pct <- round(uu_pcoa$values$Relative_eig * 100)
  
  nshape <- nlevels(s[,shape_by])
  g_uu = uu_df %>% 
    ggplot(aes(x=Axis.1, y=Axis.2, color=get(color_by), shape= get(shape_by))) +
    coord_equal() +
    theme_bw() +
    xlab(paste0("PCoA axis 1 (", uu_pct[1], "%)")) +
    ylab(paste0("PCoA axis 2 (", uu_pct[2], "%)")) +
    #scale_shape_discrete(name=sub("_", " ", shape_by)) +
    scale_shape_manual(values = 1:nshape, name=shape_by) +
    scale_color_discrete(name=color_by) +
    #scale_color_brewer(name=color_by) +
    #scale_color_gradient(name = color_by,low = 'red', high = 'blue') +
    ggtitle(gtitle) +
    geom_point()
  
  return(g_uu)
}

make_pcoa_plot_discrete_barcodeID <- function(uu, s, shape_by, color_by, gtitle="") {
  uu_pcoa <- pcoa(uu)
  uu_df <- merge(s, uu_pcoa$vectors[, 1:5], by.x="BarcodeID", by.y="row.names")
  uu_pct <- round(uu_pcoa$values$Relative_eig * 100)
  
  nshape <- nlevels(s[,shape_by])
  g_uu = uu_df %>% 
    ggplot(aes(x=Axis.1, y=Axis.2, color=get(color_by), shape= get(shape_by))) +
    coord_equal() +
    theme_bw() +
    xlab(paste0("PCoA axis 1 (", uu_pct[1], "%)")) +
    ylab(paste0("PCoA axis 2 (", uu_pct[2], "%)")) +
    #scale_shape_discrete(name=sub("_", " ", shape_by)) +
    scale_shape_manual(values = 1:nshape, name=shape_by) +
    #scale_color_discrete(name=color_by) +
    scale_color_brewer(palette = "Set1", name=color_by) +
    #scale_color_gradient(name = color_by,low = 'red', high = 'blue') +
    ggtitle(gtitle) +
    geom_point()
  
  return(g_uu)
}
make_pcoa_plot_discrete_barcodeID_matchid <- function(uu, s, shape_by, color_by, gtitle="") {
  uu_pcoa <- pcoa(uu)
  uu_df <- merge(s, uu_pcoa$vectors[, 1:5], by.x="BarcodeID", by.y="row.names")
  uu_pct <- round(uu_pcoa$values$Relative_eig * 100)
  
  nshape <- nlevels(s[,shape_by])
  g_uu = 
    uu_df %>% 
    ggplot(aes(x=Axis.1, y=Axis.2, color=get(color_by), shape=get(shape_by), group = match_id)) +
    coord_equal() +
    theme_bw() +
    xlab(paste0("PCoA axis 1 (", uu_pct[1], "%)")) +
    ylab(paste0("PCoA axis 2 (", uu_pct[2], "%)")) +
    scale_shape_manual(values = 1:nshape, name=shape_by) +
    scale_color_discrete(name=color_by) +
    ggtitle(gtitle) +
    geom_point() + 
    geom_line()
  
  return(g_uu)
}


make_pcoa_plot_continuous <- function(uu, s, shape_by, color_by, gtitle="") {
  uu_pcoa <- pcoa(uu)
  uu_df <- merge(s, uu_pcoa$vectors[, 1:5], by.x="SampleID", by.y="row.names")
  uu_pct <- round(uu_pcoa$values$Relative_eig * 100)
  
  g_uu = uu_df %>% 
    ggplot(aes(x=Axis.1, y=Axis.2, color=get(color_by), shape= get(shape_by))) +
    coord_equal() +
    theme_bw() +
    xlab(paste0("PCoA axis 1 (", uu_pct[1], "%)")) +
    ylab(paste0("PCoA axis 2 (", uu_pct[2], "%)")) +
    scale_shape_discrete(name=sub("_", " ", shape_by)) +
    scale_color_gradient(name = color_by,low = 'red', high = 'blue') +
    ggtitle(gtitle) +
    geom_point()
  
  return(g_uu)
}

make_pcoa_plot_continuous_barcodeID <- function(uu, s, shape_by, color_by, gtitle="") {
  uu_pcoa <- pcoa(uu)
  uu_df <- merge(s, uu_pcoa$vectors[, 1:5], by.x="BarcodeID", by.y="row.names")
  uu_pct <- round(uu_pcoa$values$Relative_eig * 100)
  
  g_uu = uu_df %>% 
    ggplot(aes(x=Axis.1, y=Axis.2, color=get(color_by), shape= get(shape_by))) +
    coord_equal() +
    theme_bw() +
    xlab(paste0("PCoA axis 1 (", uu_pct[1], "%)")) +
    ylab(paste0("PCoA axis 2 (", uu_pct[2], "%)")) +
    scale_shape_discrete(name=sub("_", " ", shape_by)) +
    scale_color_gradient(name = color_by,low = 'red', high = 'blue') +
    ggtitle(gtitle) +
    geom_point()
  
  return(g_uu)
}

tidy_rank_test_anca <- function(stemp, comp="Neg vs MPO"){
  ret <- stemp %>%
    droplevels() %>%
    #mutate(anca_elisa = factor(anca_elisa, levels=c("Neg","MPO"))) %>%   
    group_by(Taxon) %>%
    do(tidy(wilcox.test(Proportion ~ anca_elisa, data=.))) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    mutate(comparision = comp) %>%
    filter(p.value <= 0.05) %>%
    dplyr::select(-method) %>%
    dplyr::select(-alternative)
  
  if (dim(ret)[1] > 0){
    taxa <- ret$Taxon %>% droplevels() %>% as.character()
    fig <- stemp %>%
      filter(Taxon %in% taxa) %>%
      droplevels() %>%
      ggplot(aes(x=anca_elisa, y=Proportion, color=anca_elisa)) + 
      geom_boxplot(coef = 10) +
      geom_quasirandom(group = 1) +
      scale_color_brewer(palette = "Set2") +
      scale_y_log10() + 
      facet_wrap(~ TaxonLabel, scales = "free_y", nrow = 1) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw()
    print(fig)
  }
  
  ret
  
}
tidy_rank_test_gpa <- function(stemp, comp="Neg vs MPO"){
  ret <- stemp %>%
    droplevels() %>%
    group_by(Taxon) %>%
    do(tidy(wilcox.test(Proportion ~ GPA_activeness, data=.))) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    mutate(comparision = comp) %>%
    filter(p.value <= 0.05) %>%
    dplyr::select(-method) %>%
    dplyr::select(-alternative)
  
  if (dim(ret)[1] > 0){
    taxa <- ret$Taxon %>% droplevels() %>% as.character()
    fig <- stemp %>%
      filter(Taxon %in% taxa) %>%
      droplevels() %>%
      ggplot(aes(x=GPA_activeness, y=Proportion, color=GPA_activeness)) + 
      geom_boxplot(coef = 10) +
      geom_quasirandom(group = 1) +
      scale_color_brewer(palette = "Set2") +
      scale_y_log10() + 
      facet_wrap(~ TaxonLabel, scales = "free_y", nrow = 1) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw()
    print(fig)
  }
  
  ret
  
}

quantity_duration <- function(yrs, step=0.5){
  # step size 5 
  yrs <- as.numeric(yrs)
  
  low <- floor(range(yrs)[1])
  hi <- ceiling(range(yrs)[2])
  cands <- seq(low, hi, step) 
  
  yrs_new<- sapply(yrs, function(x){cands[which.min(abs(cands-x))]})
  
  return(yrs_new)
}


alpha_plot_continuous <- function(s_toTest, x, y, gtitle){
  rank_cor <- function(x, y){
    corr_mod <- cor.test(x, y, method = "spearman")
    corr_p <- paste("p = ", round(corr_mod$p.value, digits=3), sep="")
    corr_r <- paste("r = ", round(corr_mod$estimate, digits = 3), sep="")
    paste(corr_r,"(", corr_p, ")")
  }
  
  label = as.character(rank_cor(s_toTest[[x]], s_toTest[[y]]))
  xx = 0 #median(s_toPlot[[x]])
  yy = median(s_toPlot[[y]])
  
  gtitle <- paste(gtitle, "\n(", label,")")
  s_toTest %>%
    ggplot(aes(x=get(x), y=get(y))) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +
    labs(y=y, x=x) +
    theme_bw() + 
    ggtitle(gtitle) #+geom_text(aes(x = xx, y = yy, label=label), parse = FALSE)
}

pcoa_plot_continuous <- function(uu, s_toPlot, color_by, shape_by, gtitle="") {
  uu_pcoa <- pcoa(uu)
  uu_df <- merge(s_toPlot, uu_pcoa$vectors[, 1:5], by.x="SampleID", by.y="row.names")
  uu_pct <- round(uu_pcoa$values$Relative_eig * 100)
  
  nshape <- nlevels(s[,shape_by])
  uu_df %>%
    ggplot(aes(x=Axis.1, y=Axis.2,color=get(color_by), shape=get(shape_by))) +
    coord_equal() +
    theme_bw() +
    xlab(paste0("PCoA axis 1 (", uu_pct[1], "%)"))+
    ylab(paste0("PCoA axis 2 (", uu_pct[2], "%)")) +
    #scale_colour_gradient_tableau(name=color_by) + 
    scale_colour_gradient(low = "#D73027", high = "#4575B4", name=color_by) +
    scale_shape_manual(values = 1:nshape, name=shape_by) +
    #scale_color_brewer(palette = "Paired") + 
    ggtitle(gtitle) + 
    geom_point()
}

tidy_rank_test <- function(genus_cts_top_df, grp = "vdi_binary", color_by=NA){
  
  if (is.na(color_by))
    color_by <- grp
  
  expr <- lazyeval::interp(quote(! is.na(x)), x = as.name(grp)) 
  genus_cts_top_df %<>% filter_(expr) %>% droplevels()
  
  formula <- paste0("Proportion ~ ", grp)
  comparison <- paste0(levels(genus_cts_top_df[[grp]]) %>% as.character(), collapse = "-")
  
  ret1 <- genus_cts_top_df %>%
    group_by(Taxon) %>%
    do(tidy(wilcox.test(as.formula(formula), data=.))) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    mutate(comparision = comparison) %>%
    filter(p.value <= 0.05) %>%
    dplyr::select(-method) %>%
    dplyr::select(-alternative)
  
  if (dim(ret1)[1] > 0){
    taxa1 <- ret1$Taxon %>% droplevels() %>% as.character()
    nrow = ceiling(length(taxa1) / 4)
    
    fig <- genus_cts_top_df %>%
      filter(Taxon %in% taxa1) %>%
      droplevels() %>%
      ggplot(aes(x=get(grp), y=Proportion, color=get(color_by))) + 
      geom_boxplot(coef = 10) +
      geom_quasirandom(group = 1) +
      scale_color_brewer(palette = "Set2", name=color_by) +
      scale_y_log10() + 
      facet_wrap(~ Taxon, scales = "free_y", nrow = nrow) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw()
    
    print(fig)
  }
  ret1
}

tidy_rank_test_otu <- function(genus_cts_top_df, grp = "vdi_binary", pvalue=0.1, color_by=NA){
  
  if (is.na(color_by))
    color_by <- grp
  
  expr <- lazyeval::interp(quote(! is.na(x)), x = as.name(grp)) 
  genus_cts_top_df %<>% filter_(expr) %>% droplevels()
  
  formula <- paste0("Proportion ~ ", grp)
  comparison <- paste0(levels(genus_cts_top_df[[grp]]) %>% as.character(), collapse = "-")
  
  ret1 <- genus_cts_top_df %>%
    group_by(Taxon) %>%
    do(tidy(wilcox.test(as.formula(formula), data=.))) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    mutate(comparision = comparison) %>%
    filter(p.value <= pvalue) %>%
    dplyr::select(-method) %>%
    dplyr::select(-alternative)
  
  if (dim(ret1)[1] > 0){
    taxa1 <- ret1$Taxon %>% droplevels() %>% as.character()
    nrow = ceiling(length(taxa1) / 4)
    
    fig <- genus_cts_top_df %>%
      filter(Taxon %in% taxa1) %>%
      droplevels() %>%
      ggplot(aes(x=get(grp), y=Proportion, color=get(color_by))) + 
      geom_boxplot(coef = 10) +
      geom_quasirandom(group = 1) +
      scale_color_brewer(palette = "Set2", name=color_by) +
      scale_y_log10() + 
      facet_wrap(~ Taxon, scales = "free_y", nrow = nrow) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw()
    
    print(fig)
  }
  ret1
}

tidy_rank_test_multiple <- function(genus_cts_top_df, grp, color_by){
  #* i love this *
  combs <- combn(unlist(unique(genus_cts_top_df[, grp])), 2)
  num_tests <- dim(combs)[2]
  multi_comps <- lapply(1:num_tests, function(x){
    stemp <- dplyr::tbl_df(genus_cts_top_df)
    grps_toTest <- combs[,x] %>% as.character()
    expr <- lazyeval::interp(quote(x %in% y), x = as.name(grp), y = grps_toTest) 
    stemp <- stemp %>% filter_(expr) %>% droplevels()
    
    tidy_rank_test(stemp, grp,color_by)})
  
  do.call(rbind,multi_comps)
  
}

tidy_lm <- function(genus_cts_top_df, grp){
  fomula <- paste0("LogProp ~ ", grp)
  expr <- lazyeval::interp(quote(! is.na(x)), x = as.name(grp)) 
  genus_cts_top_df %<>%  filter_(expr) %>% droplevels()
  
  lm_models <- genus_cts_top_df %>%
    group_by(Taxon) %>%
    do(mod = lm(as.formula(fomula), data=.)) %>%
    ungroup()
  
  #lapply(lm_models$mod, summary)
  
  summaries <- lapply(1:length(lm_models$mod), function(x) data.frame(tidy(lm_models$mod[[x]]), taxa=lm_models$Taxon[[x]]))
  
  summaries_df <- do.call(rbind, summaries) %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value <= 0.05) %>%
    dplyr::select(taxa, everything())
  
  if (dim(summaries_df)[1] > 0){
    taxa <- summaries_df$taxa %>% as.character()
    nrow = ceiling(length(taxa) / 4)
    
    fig <- genus_cts_top_df %>%
      filter(Taxon %in% taxa) %>%
      droplevels() %>%
      ggplot(aes(x=get(grp), y=Proportion, color=get(grp))) + 
      geom_boxplot(coef = 10) +
      geom_quasirandom(group = 1) +
      scale_color_brewer(palette = "Set2", name=grp) +
      scale_y_log10() + 
      facet_wrap(~ TaxonLabel, scales = "free_y", nrow = nrow) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw() +
      ggtitle(grp) +
      theme(plot.title = element_text(hjust = 0.5))
    
    plot(fig)
  }
  summaries_df
}
tidy_lm_pvalue <- function(genus_cts_top_df, grp, pvalue=0.1){
  fomula <- paste0("LogProp ~ ", grp)
  expr <- lazyeval::interp(quote(! is.na(x)), x = as.name(grp)) 
  genus_cts_top_df %<>%  filter_(expr) %>% droplevels()
  
  lm_models <- genus_cts_top_df %>%
    group_by(Taxon) %>%
    do(mod = lm(as.formula(fomula), data=.)) %>%
    ungroup()
  
  #lapply(lm_models$mod, summary)
  
  summaries <- lapply(1:length(lm_models$mod), function(x) data.frame(tidy(lm_models$mod[[x]]), taxa=lm_models$Taxon[[x]]))
  
  summaries_df <- do.call(rbind, summaries) %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value <= pvalue) %>%
    dplyr::select(taxa, everything())
  
  if (dim(summaries_df)[1] > 0){
    taxa <- summaries_df$taxa %>% as.character()
    nrow = ceiling(length(taxa) / 4)
    
    fig <- genus_cts_top_df %>%
      filter(Taxon %in% taxa) %>%
      droplevels() %>%
      ggplot(aes(x=get(grp), y=Proportion, color=get(grp))) + 
      geom_boxplot(coef = 10) +
      geom_jitter() +
      scale_color_brewer(palette = "Set2", name=grp) +
      scale_y_log10() + 
      facet_wrap(~ Taxon, scales = "free_y", nrow = nrow) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw() +
      ggtitle(grp) +
      theme(plot.title = element_text(hjust = 0.5))
    
    plot(fig)
  }
  summaries_df
}
tidy_lm_continuous <- function(stest, grp){
  fomula <- paste0("LogProp ~ ", grp)
  expr <- lazyeval::interp(quote(! is.na(x)), x = as.name(grp)) 
  
  lm_models <- stest %>%
    filter_(expr) %>% 
    droplevels() %>%
    group_by(Taxon) %>%
    do(mod = lm(as.formula(fomula), data=.)) %>%
    ungroup()
  
  #lapply(lm_models$mod, summary)
  
  summaries <- lapply(1:length(lm_models$mod), function(x) data.frame(tidy(lm_models$mod[[x]]), taxa=lm_models$Taxon[[x]]))
  
  summaries_df <- do.call(rbind, summaries) %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value <= 0.05) %>%
    dplyr::select(taxa, everything())
  
  rank_cor <- function(x, y){
    corr_mod <- cor.test(x, y, method = "spearman")
    corr_p <- paste("p = ", round(corr_mod$p.value, digits=3), sep="")
    corr_r <- paste("r = ", round(corr_mod$estimate, digits = 3), sep="")
    paste(corr_r,"(", corr_p, ")")
  }
  
  #label = as.character(rank_cor(s_toTest[[x]], s_toPlot[[y]]))
  
  if (dim(summaries_df)[1] > 0){
    taxa <- summaries_df$taxa %>% as.character()
    nrow = ceiling(length(taxa) / 4)
    fig <- stest %>%
      filter(Taxon %in% taxa) %>%
      droplevels() %>%
      ggplot(aes(x=get(grp), y=Proportion, color=get(grp))) + 
      geom_point() +
      geom_smooth(method="lm", se=FALSE) +
      scale_colour_gradient(low = "#D73027", high = "#4575B4", name=grp) + 
      scale_y_log10() + 
      facet_wrap(~ Taxon, scales = "free_y", nrow = nrow) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw() +
      ggtitle(grp)
    plot(fig)
  }
  summaries_df
}
tidy_lm_continuous_pvalue <- function(stest, grp, pvalue=0.1){
  fomula <- paste0("LogProp ~ ", grp)
  expr <- lazyeval::interp(quote(! is.na(x)), x = as.name(grp)) 
  
  lm_models <- stest %>%
    filter_(expr) %>% 
    droplevels() %>%
    group_by(Taxon) %>%
    do(mod = lm(as.formula(fomula), data=.)) %>%
    ungroup()
  
  #lapply(lm_models$mod, summary)
  
  summaries <- lapply(1:length(lm_models$mod), function(x) data.frame(tidy(lm_models$mod[[x]]), taxa=lm_models$Taxon[[x]]))
  
  summaries_df <- do.call(rbind, summaries) %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value <= pvalue) %>%
    dplyr::select(taxa, everything())
  
  rank_cor <- function(x, y){
    corr_mod <- cor.test(x, y, method = "spearman")
    corr_p <- paste("p = ", round(corr_mod$p.value, digits=3), sep="")
    corr_r <- paste("r = ", round(corr_mod$estimate, digits = 3), sep="")
    paste(corr_r,"(", corr_p, ")")
  }
  
  #label = as.character(rank_cor(s_toTest[[x]], s_toPlot[[y]]))
  
  if (dim(summaries_df)[1] > 0){
    taxa <- summaries_df$taxa %>% as.character()
    nrow = ceiling(length(taxa) / 4)
    fig <- stest %>%
      filter(Taxon %in% taxa) %>%
      droplevels() %>%
      ggplot(aes(x=get(grp), y=Proportion, color=get(grp))) + 
      geom_point() +
      geom_smooth(method="lm", se=FALSE) +
      scale_colour_gradient(low = "#D73027", high = "#4575B4", name=grp) + 
      scale_y_log10() + 
      facet_wrap(~ Taxon, scales = "free_y", nrow = nrow) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw() +
      ggtitle(grp)
    plot(fig)
  }
  summaries_df
}

tidy_glmer <- function(genus_cts_top_df, grp, color_by=NA){
  if (is.na(color_by))
    color_by <- grp
  
  fomula <- paste0("cbind(Prop, PropNeg) ~ ", grp, " + (1|household_id)")
  expr <- lazyeval::interp(quote(! is.na(x)), x = as.name(grp)) 
  
  genus_models <- genus_cts_top_df %>%
    filter_(expr) %>% 
    droplevels() %>%
    group_by(Taxon) %>%
    do(mod = glmer(as.formula(fomula), data=., family="binomial")) %>%
    ungroup()
  
  #lapply(genus_models$mod, summary)
  
  summaries <- lapply(1:length(genus_models$mod), function(x) data.frame(tidy(genus_models$mod[[x]]), taxa=genus_models$Taxon[[x]],conv=ifelse(is.null(genus_models$mod[[x]]@optinfo$conv$lme4$code), 1, genus_models$mod[[x]]@optinfo$conv$lme4$code)))
  
  summaries_df <- do.call(rbind, summaries) %>%
    filter(group == "fixed") %>%
    filter(conv == 1) %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value<0.01) %>%
    dplyr::select(which(!(colnames(.) %in% c("conv", "group")))) %>%
    dplyr::select(taxa, everything())
  
  if (dim(summaries_df)[1] > 0){
    taxa <- summaries_df$taxa %>% as.character()
    nrow = ceiling(length(taxa) / 4)
    fig <- genus_cts_top_df %>%
      filter(Taxon %in% taxa) %>%
      droplevels() %>%
      ggplot(aes(x=get(grp), y=Proportion, color=get(color_by))) + 
      geom_boxplot(coef = 10) +
      geom_quasirandom(group = 1) +
      scale_color_brewer(palette = "Set2", name=color_by) +
      scale_y_log10() + 
      facet_wrap(~ TaxonLabel, scales = "free_y", nrow = nrow) +
      labs(x = "", y = "Relative Abundance") + 
      theme_bw() +
      ggtitle(grp)
    plot(fig)
  }
  summaries_df
}
