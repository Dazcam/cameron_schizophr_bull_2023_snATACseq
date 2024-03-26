#--------------------------------------------------------------------------------------
#
#    snATACseq LDSR barcharts
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for ATAC figures

##  Load Packages  --------------------------------------------------------------------
library(tidyverse) 
library(viridis)
library(cowplot)

##  Set Variables  --------------------------------------------------------------------
LDSR_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/06LDSR/'
DATA_DIR <- paste0(LDSR_DIR, 'part_herit/baseline_v1.2/')
COND_DIR <- paste0(LDSR_DIR, 'CONDITIONAL_PEAKS/part_herit/baseline_v1.2/')
UNIQUE_DIR <- paste0(LDSR_DIR, 'unique_peaks/part_herit/baseline_v1.2/')

GWASs <- c('SCZ', 'BPD', 'ASD', 'MDD', 'ADHD', 'HEIGHT')
  
## MAIN L1 ATAC FIGS  -----------------------------------------------------------------
# Fig 5 and Fig S17
for (GWAS in GWASs) {

  SUMSTATS <- read_tsv(paste0(DATA_DIR, 'snATACseq_LDSR_', GWAS, '_baseline.v1.2_summary.tsv')) %>%
    mutate(GWAS = rep(GWAS, nrow(.)))
  assign(paste0(GWAS, '_ldsr'), SUMSTATS)
  
}

ldsr_grp_df <- rbind(SCZ_ldsr, BPD_ldsr, ASD_ldsr, MDD_ldsr, ADHD_ldsr, HEIGHT_ldsr) %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
##  separate_wider_delim(Category, '.', names = c('cell_type', 'suffix')) %>%
#  select(-suffix) %>%
  dplyr::rename(cell_type = Category) %>%
  dplyr::mutate(cell_type = str_replace(cell_type, 'progenitor', 'Progenitor')) %>%
  mutate(across('cell_type', str_replace, 'GE', 'GE-N'))

for (GWAS in GWASs) {
  
  ldsr_grp_df_subset <- ldsr_grp_df %>% filter(GWAS == (!!GWAS))
  

  ldsr_grp_df_subset$cell_type <- factor(ldsr_grp_df_subset$cell_type, 
                                  levels=c("CGE-N", "LGE-N", "MGE-N", "Progenitor"))
  
  LDSR_PLOT <- ggplot(data = ldsr_grp_df_subset, aes(x = LDSR, y = factor(cell_type, rev(levels(factor(cell_type)))), 
                             fill = cell_type)) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    geom_vline(xintercept=-log10(0.05/4), linetype = "dashed", color = "black") +
    geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    ggtitle(GWAS) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          legend.position = "none") +
    xlab(expression(-log[10](P))) +
    ylab('Cell type') +
    xlim(0, 8) +
    scale_fill_manual(values = c('#6098ab','#f18e2a', '#e1575a', '#58a14e', '#75b7b2', 
                                 '#edc949', '#b07aa1', '#ff9ca7', '#9c755f', '#bab0ab'))
  
  assign(paste0('ldsr_', GWAS, '_plot'), LDSR_PLOT, .GlobalEnv)

}

##  Fig S16 - Unique peaks  ------------
Fig_S16 <- read_tsv(paste0(UNIQUE_DIR, 'snATACseq_LDSR_SCZ_baseline.v1.2_summary.tsv')) %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
  dplyr::select(Category, LDSR) %>%
  mutate(across('Category', str_replace_all, '$', '-N')) %>%
  ggplot(aes(x = LDSR, y = factor(Category, rev(levels(factor(Category)))), 
             fill = Category)) +
  geom_bar(stat = "identity", color = 'black', position = "dodge") +
  geom_vline(xintercept=-log10(0.05/3), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  ggtitle(NULL) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.position = "none") +
  xlab(expression(-log[10](P))) +
  ylab('Cell type') +
  xlim(0, 8) +
  scale_fill_manual(values = c('#6098ab','#f18e2a', '#e1575a', '#58a14e', '#75b7b2', 
                               '#edc949', '#b07aa1', '#ff9ca7', '#9c755f', '#bab0ab'))

##  LDSR conditional  ----------
for (GWAS in 'SCZ') {
  
  SUMSTATS <- read_tsv(paste0(COND_DIR, 'snATACseq_LDSR_', GWAS, '_baseline.v1.2_summary.tsv')) %>%
    mutate(GWAS = rep(GWAS, nrow(.))) 
  assign(paste0(GWAS, '_ldsr_cond'), SUMSTATS)
  
}

## Fig S18 - L1 - Ziffra Conditional -------
ldsr_ziffra_cond_df <- rbind(SCZ_ldsr_cond) %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
  filter(!grepl("union_neurons", Category)) %>%
  filter(!grepl("progenitor", Category)) %>%
  separate_wider_delim(Category, '_ziffra_', names = c('cell_type', 'suffix')) %>%
  select(-suffix) %>%
  dplyr::mutate(across(c('cell_type'), str_replace, 'progenitor', 'Progenitor')) %>%
  mutate(across('cell_type', str_replace, 'GE', 'GE-N')) %>%
  mutate(across('cell_type', str_replace_all, '_', '-')) %>%
  mutate(across('cell_type', str_replace, '-m', '')) %>%
  mutate(across('cell_type', str_replace, '-vs-', ' cond. '))

fig_S18 <- ggplot(data = ldsr_ziffra_cond_df, aes(x = LDSR, y = factor(cell_type, rev(levels(factor(cell_type)))), 
                                                                fill = cell_type)) +
  geom_bar(stat = "identity", color = 'black', position = "dodge") +
  geom_vline(xintercept=-log10(0.05/4), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  ggtitle(NULL) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.position = "none") +
  xlab(expression(-log[10](P))) +
  ylab('Cell type') +
  xlim(0, 6) +
  scale_fill_manual(values = c('#6098ab', '#6098ab', '#6098ab',
                               '#f18e2a', '#f18e2a', '#f18e2a', 
                               '#e1575a', '#e1575a', '#e1575a'))

### PLOTS -----
# Fig 5 - SCZ only
Fig_5 <- ldsr_SCZ_plot + ggtitle(NULL)

# Fig S17 - All GWAS
Fig_S17 <- plot_grid(ldsr_SCZ_plot + ggtitle('Schizophrenia'), 
                     ldsr_ADHD_plot + ggtitle('ADHD'), 
                     ldsr_ASD_plot + ggtitle('Autism'), 
                     ldsr_BPD_plot + ggtitle('Bipolar Disorder'), 
                     ldsr_MDD_plot + ggtitle('Major Depressive Disorder'), 
                     ldsr_HEIGHT_plot + ggtitle('Height'))


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


