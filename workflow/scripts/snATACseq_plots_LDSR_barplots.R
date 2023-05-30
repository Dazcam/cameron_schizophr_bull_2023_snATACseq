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
  separate_wider_delim(Category, '.', names = c('cell_type', 'suffix')) %>%
  select(-suffix) %>%
  dplyr::mutate(across(c('cell_type'), str_replace, 'progenitor', 'Progenitor')) %>%
  mutate(across('cell_type', str_replace, 'E', 'E-InN'))

ldsr_grp_df$cell_type <- factor(ldsr_grp_df$cell_type, 
                                levels=c("CGE-InN", "MGE-InN", "LGE-InN", "Progenitor"))
ldsr_grp_df$GWAS <- factor(ldsr_grp_df$GWAS, levels=c('SCZ', 'HEIGHT'))

ldsr_plot <- ggplot(ldsr_grp_df, aes(fill = cell_type, y = LDSR, x = GWAS)) + 
  geom_bar(position = position_dodge2(reverse = TRUE), stat = "identity") +
  theme_bw() +
  coord_flip() +
  facet_grid(cols = vars(REGION), rows = vars(GWAS), scales = 'free_y') +
  ylab(expression(-log[10](P))) +
  xlab(NULL) +
  geom_hline(yintercept=-log10(0.05/8), linetype = "dashed", color = "black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        panel.spacing.x = unit(2, "lines")) +
  scale_fill_manual(values = c('#6098ab','#f18e2a', '#e1575a', '#58a14e', '#75b7b2', 
                               '#edc949', '#b07aa1', '#ff9ca7', '#9c755f', '#bab0ab'))


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
  mutate(first_letter = substr(suffix, 1, 1)) %>%
  mutate(cell_type = paste(cell_type, first_letter, sep = '_')) %>%
  dplyr::mutate(across(c('cell_type'), str_replace, 'progenitor', 'Progenitor')) %>%
  mutate(across('cell_type', str_replace, 'GE', 'GE-InN'))
  
#  select(-suffix) %>%
#  dplyr::mutate(across(c('cell_type'), str_replace, 'progenitor', 'Progenitor')) %>%
#  mutate(across('cell_type', str_replace, 'E', 'E-InN'))

ldsr_cond_plot <- ggplot(ldsr_grp_cond_df, aes(fill = cell_type, y = LDSR, x = GWAS)) + 
  geom_bar(position = position_dodge2(reverse = TRUE), stat = "identity") +
  theme_bw() +
  coord_flip() +
  facet_grid(cols = vars(REGION), rows = vars(GWAS), scales = 'free_y') +
  ylab(expression(-log[10](P))) +
  xlab(NULL) +
  geom_hline(yintercept=-log10(0.05/8), linetype = "dashed", color = "black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        panel.spacing.x = unit(2, "lines")) +
  scale_fill_manual(values = c('#6098ab', '#6098ab', '#6098ab', '#6098ab', '#6098ab', '#6098ab',
                               '#f18e2a', '#f18e2a', '#f18e2a', '#f18e2a', '#f18e2a', '#f18e2a',
                               '#e1575a', '#e1575a', '#e1575a', '#e1575a', '#e1575a', '#e1575a',  
                               '#58a14e', '#58a14e', '#58a14e', '#58a14e', '#58a14e', '#58a14e'))
                    
                    

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


