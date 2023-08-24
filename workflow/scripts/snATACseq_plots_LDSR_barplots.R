#--------------------------------------------------------------------------------------
#
#    snATACseq LDSR barcharts
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse) 
library(viridis)
library(cowplot)

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/06LDSR/part_herit/baseline_v1.2/'
COND_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/06LDSR/ZIFFRA_PEAKS/part_herit/baseline_v1.2/'

GWASs <- c('SCZ', 'HEIGHT')
  
##  Basic tests  ----------------------------------------------------------------------
for (GWAS in GWASs) {

  SUMSTATS <- read_tsv(paste0(DATA_DIR, 'snATACseq_LDSR_', GWAS, '_baseline.v1.2_summary.tsv')) %>%
    mutate(GWAS = rep(GWAS, nrow(.))) %>%
    mutate(REGION = ifelse(str_detect(Category, 'FC')  == TRUE, 'Frontal Cortex', 'Ganglionic eminence'))
  assign(paste0(GWAS, '_ldsr'), SUMSTATS)
  
}

ldsr_grp_df <- rbind(SCZ_ldsr, HEIGHT_ldsr) %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
##  separate_wider_delim(Category, '.', names = c('cell_type', 'suffix')) %>%
#  select(-suffix) %>%
  rename(cell_type = Category) %>%
  dplyr::mutate(across(c('cell_type'), str_replace, 'progenitor', 'Progenitor')) %>%
  mutate(across('cell_type', str_replace, 'GE', 'GE-N'))

for (GWAS in GWASs) {
  
  ldsr_grp_df_subset <- ldsr_grp_df %>% filter(GWAS == (!!GWAS))
  

  ldsr_grp_df_subset$cell_type <- factor(ldsr_grp_df_subset$cell_type, 
                                  levels=c("CGE-N", "MGE-N", "LGE-N", "Progenitor"))
  
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
  
  assign(paste0('ldsr_', GWAS, '_plot'), LDSR_PLOT, .GlobalEnv)

}
  
plot_grid(ldsr_SCZ_plot, ldsr_HEIGHT_plot, labels = c('AUTO'), label_size = 20)


# Color scheme - to match UMAPs
# scale_fill_manual(values = c("#1F78B4", "#FFD92F", "#E31A1C", "#33A02C", "#FF7F00",
#                              "#009F75", "#88C6ED", "#EA6A47", "#394BA0", "#D54799")

##  Ziffra conditional  ---------------------------------------------------------------
for (GWAS in GWASs) {
  
  SUMSTATS <- read_tsv(paste0(COND_DIR, 'snATACseq_LDSR_', GWAS, '_baseline.v1.2_summary.tsv')) %>%
    mutate(GWAS = rep(GWAS, nrow(.))) %>%
    mutate(REGION = ifelse(str_detect(Category, 'FC')  == TRUE, 'Frontal Cortex', 'Ganglionic eminence'))
  assign(paste0(GWAS, '_ldsr_cond'), SUMSTATS)
  
}

ldsr_grp_cond_df <- rbind(SCZ_ldsr_cond, HEIGHT_ldsr_cond) %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
  separate_wider_delim(Category, '_ziffra_', names = c('cell_type', 'suffix')) %>%
  select(-suffix) %>%
  dplyr::mutate(across(c('cell_type'), str_replace, 'progenitor', 'Progenitor')) %>%
  mutate(across('cell_type', str_replace, 'GE', 'GE-N')) %>%
  mutate(across('cell_type', str_replace_all, '_', '-')) %>%
  filter(GWAS == 'SCZ') %>%
  mutate(across('cell_type', str_replace, '-m', '')) %>%
  mutate(across('cell_type', str_replace, '-vs-', ' vs. '))

  
#  select(-suffix) %>%
#  dplyr::mutate(across(c('cell_type'), str_replace, 'progenitor', 'Progenitor')) %>%
#  mutate(across('cell_type', str_replace, 'E', 'E-InN'))

ldsr_cond_plot <- ggplot(data = ldsr_grp_cond_df, aes(x = LDSR, y = factor(cell_type, rev(levels(factor(cell_type)))), 
                                                        fill = cell_type)) +
  geom_bar(stat = "identity", color = 'black', position = "dodge") +
  geom_vline(xintercept=-log10(0.05/12), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  ggtitle('SCZ') +
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
  scale_fill_manual(values = c('#6098ab', '#6098ab', '#6098ab',
                               '#f18e2a', '#f18e2a', '#f18e2a', 
                               '#e1575a', '#e1575a', '#e1575a',  
                               '#58a14e', '#58a14e', '#58a14e'))
                    
                    

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


