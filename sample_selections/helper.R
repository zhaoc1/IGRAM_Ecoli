library(magrittr)
library(dplyr)
library(pheatmap)

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
    arrange(both_kept) %>% #arrange_(.dots = top_taxa[1])
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
    geom_bar(stat="identity", width = ) +
    scale_fill_manual(values = c(fill, "#CCCCCC")) + #, breaks = NULL, drop=TRUE
    scale_y_continuous(limits = c(0,1.01), expand=c(0,0.1)) + 
    #facet_wrap( ~ factor(get(grp1)), scales = "free_x", nrow=nrow) + #factor(get(grp1))
    facet_grid( ~ factor(get(grp1)), scales = "free_x", space="free_x") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
  #theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
}

taxonomic_barchart_grid <- function(cts, props,s_toTest,fill=taxonomic_colors, grp=c("Site_Status"), repid="Patient", option=1){
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
  
  #rows_to_keep <- filter_low_coverage(props, perc_cutoff=0.5)
  #rows_to_keep <- apply(props,1,max) >= 0.8
  
  ## very nice: wide to long
  ## sort the SampleID by NumReads
  props_df <- s_toTest %>%
    merge(t(top_props), by.x = "SampleID", by.y = "row.names") %>%
    arrange(both_kept) %>% #arrange_(.dots = top_taxa[1])
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
    ggplot(aes(x=factor(get(repid)), y=Proportion, fill=Taxon, width=0.9)) +
    geom_bar(stat = "identity") + 
    theme_bw() +
    scale_fill_manual(values = c(fill, "#CCCCCC")) + #, breaks = NULL, drop=TRUE
    scale_y_continuous(limits = c(0,1.01), expand=c(0,0.1)) + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),strip.text.y = element_text(angle=0)) +
    labs(x = repid, y = "Relative Abundance") +
    facet_grid(factor(get(grp1)) ~ factor(get(grp2)), scales = "free_x",space="free") 
}

heatmap_grouped <- function(genus_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8,option=1){
  
  color = saturated_rainbow(101, 0.6)
  breaks = c(0, 1e-10, seq(2, 100, length.out = 100))
  
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
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, filename = fname, 
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, 
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
}

heatmap_clustered <- function(genus_props, heatmap_s, grps=c("SampleType"), fname=NULL, thre=0.1, option=1){
  
  color = saturated_rainbow(101, 0.6)
  #breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  breaks = c(0, 1e-10, seq(2, 100, length.out = 100))
  
  # clustered heatmap is mostly used to check spill over
  heatmap_s <- heatmap_s %>% arrange_(.dots=grps)
  heatmap_props <- genus_props[, heatmap_s$SampleID]
  
  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps
  
  if (option == 1)
    rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) #
  else
    rows_to_keep <- apply(heatmap_props,1,max) >= thre #
  heatmap_props <- heatmap_props[rows_to_keep,]
  
  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks,
             filename = fname, fontsize_col = 8, fontsize_row = 8, cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks,
             fontsize_col = 8, fontsize_row = 8, cellheight = 8, cellwidth = 8)
}

