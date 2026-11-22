# ==================== NULISA: DIFFERENTIAL ABUNDANCE (NULISAseqR) ====================
# Input : <de_workdir>/Raw Data.csv  (NULISAseq long-format export: SampleName, Target, NPQ)
# Output: <de_workdir>/  DE tables and volcano plots per comparison
# ---- Configuration: edit config/config.yaml (run from the repository root) ----
source("utils/load_config.R")
setwd(cfg_path("de_workdir", create = TRUE))

# Samples excluded from some comparisons (real sample names live only in config/config.yaml)
excluded_sample_1 <- cfg$exclude_samples$sample_1   # removed from the MutB comparisons ("wo_excluded_sample1")
excluded_sample_2 <- cfg$exclude_samples$sample_2   # removed from the "wo_noMut_pat" comparisons

library(NULISAseqR)
library(tidyverse)
library(ComplexHeatmap)
library(readxl)
library(ggplot2)
library(ggrepel)

Data <- read.csv("Raw Data.csv", stringsAsFactors = FALSE)

library(dplyr)  # already loaded via tidyverse

proteins_to_keep <- c(
  "ACHE","AGRN","ANXA5","APOE","APOE4","Aβ38","Aβ40","Aβ42","BACE1","BD-MAPT",
  "BD-pTau-181","BD-pTau-217","BD-pTau-231","CALB2","CCL11","CCL13","CCL17",
  "CCL2","CCL22","CCL26","CCL3","CCL4","CD40LG","CD63","CHI3L1","CHIT1",
  "CNTN2","CRH","CRP","CSF2","CST3","CX3CL1","CXCL1","CXCL10","CXCL8","DDC",
  "ENO2","FABP3","FCN2","FGF2","FLT1","FOLR1","GDF15","GDI1","GFAP","GOT1",
  "HBA1","HTT","ICAM1","IFNG","IGF1R","IGFBP7","IL10","IL12p70","IL13","IL15",
  "IL16","IL18","IL33","IL5","IL6","IL6R","IL7","IL9","KDR","KLK6","MAPT",
  "MDH1","MSLN","NEFH","NEFL","NGF","NPTX1","NPTX2","NPTXR","NPY","NRGN",
  "PARK7","PDGFRB","PDLIM5","PGF","PGK1","POSTN","PRDX6","PSEN1","pTau-181",
  "pTau-217","pTau-231","PTN","RUVBL2","S100A12","S100B","SAA1","SFRP1",
  "SFTPD","SLIT2","SMOC1","SNAP25","SNCA","SNCB","SOD1","SQSTM1","TAFA5",
  "TARDBP","TEK","TIMP3","TNF","TREM1","TREM2","UBB","UCHL1","VCAM1","VEGFA",
  "VEGFD","VGF","VSNL1","YWHAG"
)

Data <- Data %>%
  filter(Target %in% proteins_to_keep)


# Define your grouping based on SampleName — adjust pattern as needed
Data$Group <- case_when(
  grepl("Control", Data$SampleName) ~ "Control",
  grepl("NoMut", Data$SampleName) ~ "NoMut",
  grepl("MutC", Data$SampleName) ~ "MutC",
  grepl("MutBpos", Data$SampleName) ~ "MutBpos",
  grepl("DrugA", Data$SampleName) ~ "DrugA",
  grepl("DrugB-1month", Data$SampleName) ~ "DrugB-1month",
  grepl("DrugB-2month", Data$SampleName) ~ "DrugB-2month",
  grepl("DrugB-3month", Data$SampleName) ~ "DrugB-3month",
  grepl("DrugB-4month", Data$SampleName) ~ "DrugB-4month",
  grepl("DrugB-5month", Data$SampleName) ~ "DrugB-5month",
  grepl("MutB-Baseline", Data$SampleName) ~ "MutB-Baseline",
  TRUE ~ "Other"
)

# Filter for Control vs NoMut
data_cs <- Data %>% filter(Group %in% c("Control", "NoMut"))

# Reshape to wide format
data_cs_wide <- data_cs %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata
meta_cs <- data_cs %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)

# DE Analysis
lm_cs <- lmNULISAseq(
  data = as.matrix(data_cs_wide),
  sampleInfo = meta_cs,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)

NoMut_DE_proteins <- lm_cs$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )

write.csv(NoMut_DE_proteins, "NoMut_DE_proteins.csv", row.names = FALSE)

# Threshold for significance (you can adjust this)
pthresh <- 0.05

# Annotate significance
NoMut_DE_proteins$Significance <- with(NoMut_DE_proteins, 
                                          ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                 ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

NoMut_FDR_volcano_plot <- ggplot(NoMut_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  # Label only significant points for clarity
  ggrepel::geom_text_repel(data = subset(NoMut_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Save the plot as a PNG (adjust filename and size as desired)
ggsave("NoMut_FDR_volcano_plot.png", plot = NoMut_FDR_volcano_plot, width = 8, height = 6, dpi = 300)

print(NoMut_FDR_volcano_plot)

# Threshold for significance (you can adjust this)
pthresh <- 0.05

# Annotate significance
NoMut_DE_proteins$Significance <- with(NoMut_DE_proteins, 
                                 ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                        ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

NoMut_unadj_volcano_plot <- ggplot(NoMut_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  # Add labels only to significant points to reduce clutter
  ggrepel::geom_text_repel(data = subset(NoMut_DE_proteins, Significance != "Not Significant"), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Save the plot as a PNG (adjust filename and size as desired)
ggsave("NoMut_unadj_volcano_plot.png", plot = NoMut_unadj_volcano_plot, width = 8, height = 6, dpi = 300)

print(NoMut_unadj_volcano_plot)


# Filter for Control vs MutBpos
data_cc9 <- Data %>%
  filter(Group %in% c("Control", "MutBpos")) %>%
  filter(!grepl("DrugB", SampleName, ignore.case = TRUE)) %>%
  filter(!grepl("MutB-Baseline", SampleName, ignore.case = TRUE))

# Reshape to wide format
data_cc9_wide <- data_cc9 %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata
meta_cc9 <- data_cc9 %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)

# Now you can mutate the factor levels
meta_cc9 <- meta_cc9 %>%
  mutate(disease_type = factor(disease_type, levels = c("Control", "MutBpos")))


# DE Analysis
lm_cc9 <- lmNULISAseq(
  data = as.matrix(data_cc9_wide),
  sampleInfo = meta_cc9,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)


MutBpos_DE_proteins <- lm_cc9$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )

write.csv(MutBpos_DE_proteins, "MutBpos_DE_proteins.csv", row.names = FALSE)

# Threshold for significance (you can adjust this)
pthresh <- 0.05

# Annotate significance
MutBpos_DE_proteins$Significance <- with(MutBpos_DE_proteins, 
                                          ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                 ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

MutBpos_FDR_volcano_plot <- ggplot(MutBpos_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  # Label only significant points for clarity
  ggrepel::geom_text_repel(data = subset(MutBpos_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Save the plot as a PNG (adjust filename and size as desired)
ggsave("MutBpos_FDR_volcano_plot.png", plot = MutBpos_FDR_volcano_plot, width = 8, height = 6, dpi = 300)

print(MutBpos_FDR_volcano_plot)

# Threshold for significance (you can adjust this)
pthresh <- 0.05


# Annotate significance
MutBpos_DE_proteins$Significance <- with(MutBpos_DE_proteins, 
                                          ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                 ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

MutBpos_unadj_volcano_plot <- ggplot(MutBpos_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  # Add labels only to significant points to reduce clutter
  ggrepel::geom_text_repel(data = subset(MutBpos_DE_proteins, Significance != "Not Significant"), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )


# Save the plot as a PNG (adjust filename and size as desired)
ggsave("MutBpos_unadj_volcano_plot.png", plot = MutBpos_unadj_volcano_plot, width = 8, height = 6, dpi = 300)


print(MutBpos_unadj_volcano_plot)


# Filter for NoMut vs MutBpos
data_NoMut_mutB <- Data %>%
  filter(Group %in% c("NoMut", "MutBpos")) %>% 
  filter(!grepl("DrugB", SampleName, ignore.case = TRUE)) %>% 
  filter(!grepl("MutB-Baseline", SampleName, ignore.case = TRUE))

# Reshape to wide format
data_NoMut_mutB_wide <- data_NoMut_mutB %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata
meta_NoMut_mutB <- data_NoMut_mutB %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)

# Now you can mutate the factor levels
meta_NoMut_mutB <- meta_NoMut_mutB %>%
  mutate(disease_type = factor(disease_type, levels = c("NoMut", "MutBpos")))


# DE Analysis
lm_NoMut_mutB <- lmNULISAseq(
  data = as.matrix(data_NoMut_mutB_wide),
  sampleInfo = meta_NoMut_mutB,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)


NoMut_MutBpos_DE_proteins <- lm_NoMut_mutB$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )

write.csv(NoMut_MutBpos_DE_proteins, "NoMut_MutBpos_DE_proteins.csv", row.names = FALSE)

# Threshold for significance (you can adjust this)
pthresh <- 0.05

# Annotate significance
NoMut_MutBpos_DE_proteins$Significance <- with(NoMut_MutBpos_DE_proteins, 
                                         ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

NoMut_MutBpos_FDR_volcano_plot <- ggplot(NoMut_MutBpos_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  # Label only significant points for clarity
  ggrepel::geom_text_repel(data = subset(NoMut_MutBpos_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Save the plot as a PNG (adjust filename and size as desired)
ggsave("NoMut_MutBpos_FDR_volcano_plot.png", plot = NoMut_MutBpos_FDR_volcano_plot, width = 8, height = 6, dpi = 300)

print(NoMut_MutBpos_FDR_volcano_plot)

# Threshold for significance (you can adjust this)
pthresh <- 0.05

# Annotate significance
NoMut_MutBpos_DE_proteins$Significance <- with(NoMut_MutBpos_DE_proteins, 
                                         ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

NoMut_MutBpos_unadj_volcano_plot <- ggplot(NoMut_MutBpos_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  # Add labels only to significant points to reduce clutter
  ggrepel::geom_text_repel(data = subset(NoMut_MutBpos_DE_proteins, Significance != "Not Significant"), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Save the plot as a PNG (adjust filename and size as desired)
ggsave("NoMut_MutBpos_unadj_volcano_plot.png", plot = NoMut_MutBpos_unadj_volcano_plot, width = 8, height = 6, dpi = 300)

print(NoMut_MutBpos_unadj_volcano_plot)


# Filter for MutB-DrugB-1month vs MutB-Baseline groups
data_DrugB_baseline <- Data %>%
  filter(Group %in% c("DrugB-1month", "MutB-Baseline"))%>% 
  filter(!grepl("NoMut", SampleName, ignore.case = TRUE)) %>% 
  filter(!grepl("MutBpos", SampleName, ignore.case = TRUE))

# Reshape to wide format
data_DrugB_baseline_wide <- data_DrugB_baseline %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata
meta_DrugB_baseline <- data_DrugB_baseline %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)

# Mutate factor levels for group comparison
meta_DrugB_baseline <- meta_DrugB_baseline %>%
  mutate(disease_type = factor(disease_type, levels = c("MutB-Baseline", "DrugB-1month")))

# Differential Expression Analysis
lm_DrugB_baseline <- lmNULISAseq(
  data = as.matrix(data_DrugB_baseline_wide),
  sampleInfo = meta_DrugB_baseline,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)

# Select and rename relevant statistics
DrugB_Baseline_DE_proteins <- lm_DrugB_baseline$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )

# Save differential expression results to CSV
write.csv(DrugB_Baseline_DE_proteins, "DrugB_Baseline_DE_proteins.csv", row.names = FALSE)

# Threshold for significance
pthresh <- 0.05

# Annotate significance based on FDR-adjusted p-values
DrugB_Baseline_DE_proteins$Significance <- with(DrugB_Baseline_DE_proteins,
                                                    ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                           ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

# Volcano plot for FDR-adjusted p-values
DrugB_Baseline_FDR_volcano_plot <- ggplot(DrugB_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB_Baseline_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes (MutB-DrugB-1month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Save the FDR volcano plot
ggsave("DrugB_Baseline_FDR_volcano_plot.png", plot = DrugB_Baseline_FDR_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB_Baseline_FDR_volcano_plot)

# Annotate significance based on unadjusted p-values
DrugB_Baseline_DE_proteins$Significance <- with(DrugB_Baseline_DE_proteins,
                                                    ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                           ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)


# Volcano plot for unadjusted p-values
DrugB_Baseline_unadj_volcano_plot <- ggplot(DrugB_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB_Baseline_DE_proteins, Significance != "Not Significant"), 
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes (MutB-DrugB-1month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

# Save the unadjusted volcano plot
ggsave("DrugB_Baseline_unadj_volcano_plot.png", plot = DrugB_Baseline_unadj_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB_Baseline_unadj_volcano_plot)


# Filter for MutB-DrugB-2month vs MutB-Baseline groups
data_DrugB2mo_baseline <- Data %>%
  filter(Group %in% c("DrugB-2month", "MutB-Baseline")) %>% 
  filter(!grepl("NoMut", SampleName, ignore.case = TRUE)) %>% 
  filter(!grepl("MutBpos", SampleName, ignore.case = TRUE))
# Reshape to wide format
data_DrugB2mo_baseline_wide <- data_DrugB2mo_baseline %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")
# Metadata
meta_DrugB2mo_baseline <- data_DrugB2mo_baseline %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)
# Mutate factor levels for group comparison
meta_DrugB2mo_baseline <- meta_DrugB2mo_baseline %>%
  mutate(disease_type = factor(disease_type, levels = c("MutB-Baseline", "DrugB-2month")))
# Differential Expression Analysis
lm_DrugB2mo_baseline <- lmNULISAseq(
  data = as.matrix(data_DrugB2mo_baseline_wide),
  sampleInfo = meta_DrugB2mo_baseline,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)
# Select and rename relevant statistics
DrugB2mo_Baseline_DE_proteins <- lm_DrugB2mo_baseline$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )
# Save differential expression results to CSV
write.csv(DrugB2mo_Baseline_DE_proteins, "DrugB2mo_Baseline_DE_proteins.csv", row.names = FALSE)
# Threshold for significance
pthresh <- 0.05
# Annotate significance based on FDR-adjusted p-values
DrugB2mo_Baseline_DE_proteins$Significance <- with(DrugB2mo_Baseline_DE_proteins,
                                                       ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                              ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)
# Volcano plot for FDR-adjusted p-values
DrugB2mo_Baseline_FDR_volcano_plot <- ggplot(DrugB2mo_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB2mo_Baseline_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes (MutB-DrugB-2month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )
# Save the FDR volcano plot
ggsave("DrugB2mo_Baseline_FDR_volcano_plot.png", plot = DrugB2mo_Baseline_FDR_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB2mo_Baseline_FDR_volcano_plot)
# Annotate significance based on unadjusted p-values
DrugB2mo_Baseline_DE_proteins$Significance <- with(DrugB2mo_Baseline_DE_proteins,
                                                       ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                              ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)
# Volcano plot for unadjusted p-values
DrugB2mo_Baseline_unadj_volcano_plot <- ggplot(DrugB2mo_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB2mo_Baseline_DE_proteins, Significance != "Not Significant"), 
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes (MutB-DrugB-2month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )
# Save the unadjusted volcano plot
ggsave("DrugB2mo_Baseline_unadj_volcano_plot.png", plot = DrugB2mo_Baseline_unadj_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB2mo_Baseline_unadj_volcano_plot)


# Filter for MutB-DrugB-3month vs MutB-Baseline groups
data_DrugB3mo_baseline <- Data %>%
  filter(Group %in% c("DrugB-3month", "MutB-Baseline")) %>% 
  filter(!grepl("NoMut", SampleName, ignore.case = TRUE)) %>% 
  filter(!grepl("MutBpos", SampleName, ignore.case = TRUE))
# Reshape to wide format
data_DrugB3mo_baseline_wide <- data_DrugB3mo_baseline %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")
# Metadata
meta_DrugB3mo_baseline <- data_DrugB3mo_baseline %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)
# Mutate factor levels for group comparison
meta_DrugB3mo_baseline <- meta_DrugB3mo_baseline %>%
  mutate(disease_type = factor(disease_type, levels = c("MutB-Baseline", "DrugB-3month")))
# Differential Expression Analysis
lm_DrugB3mo_baseline <- lmNULISAseq(
  data = as.matrix(data_DrugB3mo_baseline_wide),
  sampleInfo = meta_DrugB3mo_baseline,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)
# Select and rename relevant statistics
DrugB3mo_Baseline_DE_proteins <- lm_DrugB3mo_baseline$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )
# Save differential expression results to CSV
write.csv(DrugB3mo_Baseline_DE_proteins, "DrugB3mo_Baseline_DE_proteins.csv", row.names = FALSE)
# Threshold for significance
pthresh <- 0.05
# Annotate significance based on FDR-adjusted p-values
DrugB3mo_Baseline_DE_proteins$Significance <- with(DrugB3mo_Baseline_DE_proteins,
                                                       ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                              ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)
# Volcano plot for FDR-adjusted p-values
DrugB3mo_Baseline_FDR_volcano_plot <- ggplot(DrugB3mo_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB3mo_Baseline_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes (MutB-DrugB-3month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )
# Save the FDR volcano plot
ggsave("DrugB3mo_Baseline_FDR_volcano_plot.png", plot = DrugB3mo_Baseline_FDR_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB3mo_Baseline_FDR_volcano_plot)
# Annotate significance based on unadjusted p-values
DrugB3mo_Baseline_DE_proteins$Significance <- with(DrugB3mo_Baseline_DE_proteins,
                                                       ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                              ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)
# Volcano plot for unadjusted p-values
DrugB3mo_Baseline_unadj_volcano_plot <- ggplot(DrugB3mo_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB3mo_Baseline_DE_proteins, Significance != "Not Significant"), 
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes (MutB-DrugB-3month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )
# Save the unadjusted volcano plot
ggsave("DrugB3mo_Baseline_unadj_volcano_plot.png", plot = DrugB3mo_Baseline_unadj_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB3mo_Baseline_unadj_volcano_plot)


# Filter for MutB-DrugB-4month vs MutB-Baseline groups
data_DrugB4mo_baseline <- Data %>%
  filter(Group %in% c("DrugB-4month", "MutB-Baseline")) %>% 
  filter(!grepl("NoMut", SampleName, ignore.case = TRUE)) %>% 
  filter(!grepl("MutBpos", SampleName, ignore.case = TRUE))
# Reshape to wide format
data_DrugB4mo_baseline_wide <- data_DrugB4mo_baseline %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")
# Metadata
meta_DrugB4mo_baseline <- data_DrugB4mo_baseline %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)
# Mutate factor levels for group comparison
meta_DrugB4mo_baseline <- meta_DrugB4mo_baseline %>%
  mutate(disease_type = factor(disease_type, levels = c("MutB-Baseline", "DrugB-4month")))
# Differential Expression Analysis
lm_DrugB4mo_baseline <- lmNULISAseq(
  data = as.matrix(data_DrugB4mo_baseline_wide),
  sampleInfo = meta_DrugB4mo_baseline,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)
# Select and rename relevant statistics
DrugB4mo_Baseline_DE_proteins <- lm_DrugB4mo_baseline$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )
# Save differential expression results to CSV
write.csv(DrugB4mo_Baseline_DE_proteins, "DrugB4mo_Baseline_DE_proteins.csv", row.names = FALSE)
# Threshold for significance
pthresh <- 0.05
# Annotate significance based on FDR-adjusted p-values
DrugB4mo_Baseline_DE_proteins$Significance <- with(DrugB4mo_Baseline_DE_proteins,
                                                       ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                              ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)
# Volcano plot for FDR-adjusted p-values
DrugB4mo_Baseline_FDR_volcano_plot <- ggplot(DrugB4mo_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB4mo_Baseline_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes (MutB-DrugB-4month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )
# Save the FDR volcano plot
ggsave("DrugB4mo_Baseline_FDR_volcano_plot.png", plot = DrugB4mo_Baseline_FDR_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB4mo_Baseline_FDR_volcano_plot)
# Annotate significance based on unadjusted p-values
DrugB4mo_Baseline_DE_proteins$Significance <- with(DrugB4mo_Baseline_DE_proteins,
                                                       ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                              ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)
# Volcano plot for unadjusted p-values
DrugB4mo_Baseline_unadj_volcano_plot <- ggplot(DrugB4mo_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB4mo_Baseline_DE_proteins, Significance != "Not Significant"), 
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes (MutB-DrugB-4month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )
# Save the unadjusted volcano plot
ggsave("DrugB4mo_Baseline_unadj_volcano_plot.png", plot = DrugB4mo_Baseline_unadj_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB4mo_Baseline_unadj_volcano_plot)


# Filter for MutB-DrugB-5month vs MutB-Baseline groups
data_DrugB5mo_baseline <- Data %>%
  filter(Group %in% c("DrugB-5month", "MutB-Baseline")) %>% 
  filter(!grepl("NoMut", SampleName, ignore.case = TRUE)) %>% 
  filter(!grepl("MutBpos", SampleName, ignore.case = TRUE))
# Reshape to wide format
data_DrugB5mo_baseline_wide <- data_DrugB5mo_baseline %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")
# Metadata
meta_DrugB5mo_baseline <- data_DrugB5mo_baseline %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)
# Mutate factor levels for group comparison
meta_DrugB5mo_baseline <- meta_DrugB5mo_baseline %>%
  mutate(disease_type = factor(disease_type, levels = c("MutB-Baseline", "DrugB-5month")))
# Differential Expression Analysis
lm_DrugB5mo_baseline <- lmNULISAseq(
  data = as.matrix(data_DrugB5mo_baseline_wide),
  sampleInfo = meta_DrugB5mo_baseline,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)
# Select and rename relevant statistics
DrugB5mo_Baseline_DE_proteins <- lm_DrugB5mo_baseline$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )
# Save differential expression results to CSV
write.csv(DrugB5mo_Baseline_DE_proteins, "DrugB5mo_Baseline_DE_proteins.csv", row.names = FALSE)
# Threshold for significance
pthresh <- 0.05
# Annotate significance based on FDR-adjusted p-values
DrugB5mo_Baseline_DE_proteins$Significance <- with(DrugB5mo_Baseline_DE_proteins,
                                                       ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                              ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)
# Volcano plot for FDR-adjusted p-values
DrugB5mo_Baseline_FDR_volcano_plot <- ggplot(DrugB5mo_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB5mo_Baseline_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes (MutB-DrugB-5month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )
# Save the FDR volcano plot
ggsave("DrugB5mo_Baseline_FDR_volcano_plot.png", plot = DrugB5mo_Baseline_FDR_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB5mo_Baseline_FDR_volcano_plot)
# Annotate significance based on unadjusted p-values
DrugB5mo_Baseline_DE_proteins$Significance <- with(DrugB5mo_Baseline_DE_proteins,
                                                       ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                              ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)
# Volcano plot for unadjusted p-values
DrugB5mo_Baseline_unadj_volcano_plot <- ggplot(DrugB5mo_Baseline_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(DrugB5mo_Baseline_DE_proteins, Significance != "Not Significant"), 
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes (MutB-DrugB-5month vs MutB-Baseline)",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )
# Save the unadjusted volcano plot
ggsave("DrugB5mo_Baseline_unadj_volcano_plot.png", plot = DrugB5mo_Baseline_unadj_volcano_plot, width = 8, height = 6, dpi = 300)
print(DrugB5mo_Baseline_unadj_volcano_plot)



# Filter for Control vs MutBpos, excluding excluded_sample_1
data_cc9 <- Data %>%
  filter(Group %in% c("Control", "MutBpos")) %>%
  filter(SampleName != excluded_sample_1) %>%
  filter(!grepl("DrugB", SampleName, ignore.case = TRUE)) %>%
  filter(!grepl("MutB-Baseline", SampleName, ignore.case = TRUE))

# Reshape to wide format
data_cc9_wide <- data_cc9 %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata
meta_cc9 <- data_cc9 %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)

# Now you can mutate the factor levels
meta_cc9 <- meta_cc9 %>%
  mutate(disease_type = factor(disease_type, levels = c("Control", "MutBpos")))

# DE Analysis
lm_cc9 <- lmNULISAseq(
  data = as.matrix(data_cc9_wide),
  sampleInfo = meta_cc9,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)

MutBpos_DE_proteins <- lm_cc9$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )

write.csv(MutBpos_DE_proteins, "MutBpos_DE_proteins_wo_excluded_sample1.csv", row.names = FALSE)

pthresh <- 0.05

# Annotate significance (FDR)
MutBpos_DE_proteins$Significance <- with(MutBpos_DE_proteins, 
                                         ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

MutBpos_FDR_volcano_plot <- ggplot(MutBpos_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(MutBpos_DE_proteins, Significance != "Not Significant"), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("MutBpos_FDR_volcano_plot_wo_excluded_sample1.png", plot = MutBpos_FDR_volcano_plot, width = 8, height = 6, dpi = 300)
print(MutBpos_FDR_volcano_plot)

# Annotate significance (unadjusted)
MutBpos_DE_proteins$Significance <- with(MutBpos_DE_proteins, 
                                         ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

MutBpos_unadj_volcano_plot <- ggplot(MutBpos_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(MutBpos_DE_proteins, Significance != "Not Significant"), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("MutBpos_unadj_volcano_plot_wo_excluded_sample1.png", plot = MutBpos_unadj_volcano_plot, width = 8, height = 6, dpi = 300)
print(MutBpos_unadj_volcano_plot)

# Filter for NoMut vs MutBpos, excluding excluded_sample_1
data_NoMut_mutB <- Data %>%
  filter(Group %in% c("NoMut", "MutBpos")) %>%
  filter(SampleName != excluded_sample_1) %>%
  filter(!grepl("DrugB", SampleName, ignore.case = TRUE)) %>%
  filter(!grepl("MutB-Baseline", SampleName, ignore.case = TRUE))

# Reshape to wide format
data_NoMut_mutB_wide <- data_NoMut_mutB %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata
meta_NoMut_mutB <- data_NoMut_mutB %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group)

# Now you can mutate the factor levels
meta_NoMut_mutB <- meta_NoMut_mutB %>%
  mutate(disease_type = factor(disease_type, levels = c("NoMut", "MutBpos")))

# DE Analysis
lm_NoMut_mutB <- lmNULISAseq(
  data = as.matrix(data_NoMut_mutB_wide),
  sampleInfo = meta_NoMut_mutB,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)

NoMut_MutBpos_DE_proteins <- lm_NoMut_mutB$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )

write.csv(NoMut_MutBpos_DE_proteins, "NoMut_MutBpos_DE_proteins_wo_excluded_sample1.csv", row.names = FALSE)

pthresh <- 0.05

# Annotate significance (FDR)
NoMut_MutBpos_DE_proteins$Significance <- with(NoMut_MutBpos_DE_proteins, 
                                                  ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                         ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

NoMut_MutBpos_FDR_volcano_plot <- ggplot(NoMut_MutBpos_DE_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(NoMut_MutBpos_DE_proteins, Significance != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("NoMut_MutBpos_FDR_volcano_plot_wo_excluded_sample1.png", plot = NoMut_MutBpos_FDR_volcano_plot, width = 8, height = 6, dpi = 300)
print(NoMut_MutBpos_FDR_volcano_plot)

# Annotate significance (unadjusted)
NoMut_MutBpos_DE_proteins$Significance <- with(NoMut_MutBpos_DE_proteins, 
                                                  ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
                                                         ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

NoMut_MutBpos_unadj_volcano_plot <- ggplot(NoMut_MutBpos_DE_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(NoMut_MutBpos_DE_proteins, Significance != "Not Significant"), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("NoMut_MutBpos_unadj_volcano_plot_wo_excluded_sample1.png", plot = NoMut_MutBpos_unadj_volcano_plot, width = 8, height = 6, dpi = 300)
print(NoMut_MutBpos_unadj_volcano_plot)

# Filter for Control vs MutBpos including all MutB-Baseline samples
data_cc9_with_baseline <- Data %>%
  filter(Group %in% c("Control", "MutBpos","MutB-Baseline")) %>%
  filter(SampleName != excluded_sample_1) %>%
  filter(!grepl("DrugB", SampleName, ignore.case = TRUE))
# Removed filter excluding MutB-Baseline samples

# Reshape to wide format
data_cc9_wide_with_baseline <- data_cc9_with_baseline %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata
# Correct grouping: treat 'MutB-Baseline' as 'MutBpos'
meta_cc9_with_baseline <- data_cc9_with_baseline %>%
  select(SampleName, Group) %>%
  distinct() %>%
  mutate(
    disease_type = ifelse(Group == "MutB-Baseline", "MutBpos", Group)
  ) %>%
  rename(sampleName_long = SampleName) %>%
  mutate(disease_type = factor(disease_type, levels = c("Control", "MutBpos")))


# Mutate factor levels
meta_cc9_with_baseline <- meta_cc9_with_baseline %>%
  mutate(disease_type = factor(disease_type, levels = c("Control", "MutBpos")))

# DE Analysis
lm_cc9_with_baseline <- lmNULISAseq(
  data = as.matrix(data_cc9_wide_with_baseline),
  sampleInfo = meta_cc9_with_baseline,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)

MutBpos_DE_with_basleine_proteins <- lm_cc9_with_baseline$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )

write.csv(MutBpos_DE_with_basleine_proteins, "MutBpos_DE_proteins_with_MutBBaseline.csv", row.names = FALSE)

pthresh <- 0.05

# Annotate significance (FDR)
MutBpos_DE_with_basleine_proteins$Significance <- with(MutBpos_DE_with_basleine_proteins, 
                                         ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
                                                ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

MutBpos_FDR_with_basleine_volcano_plot <- ggplot(MutBpos_DE_with_basleine_proteins, aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(MutBpos_DE_with_basleine_proteins, Significance != "Not Significant"), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("MutBpos_FDR_volcano_plot_with_MutBBaseline.png", plot = MutBpos_FDR_with_basleine_volcano_plot, width = 8, height = 6, dpi = 300)
print(MutBpos_FDR_with_basleine_volcano_plot)

# Annotate significance for unadjusted p-values
MutBpos_DE_with_basleine_proteins$Significance_unadj <- with(
  MutBpos_DE_with_basleine_proteins,
  ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
         ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

# Plot using unadjusted p-value annotations
MutBpos_with_basleine_unadj_volcano_plot <- ggplot(MutBpos_DE_with_basleine_proteins, aes(x = Coef, y = -log10(Pval_unadj), color = Significance_unadj, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(MutBpos_DE_with_basleine_proteins, Significance_unadj != "Not Significant"), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("MutBpos_unadj_volcano_plot_with_MutBBaseline.png", plot = MutBpos_with_basleine_unadj_volcano_plot, width = 8, height = 6, dpi = 300)
print(MutBpos_with_basleine_unadj_volcano_plot)


# Filter including MutB-Baseline
data_NoMut_mutB_with_baseline <- Data %>%
  filter(Group %in% c("NoMut", "MutBpos", "MutB-Baseline")) %>%
  filter(SampleName != excluded_sample_1) %>%
  filter(!grepl("DrugB", SampleName, ignore.case = TRUE))

# Reshape to wide format
data_NoMut_mutB_wide_with_baseline <- data_NoMut_mutB_with_baseline %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata: treat MutB-Baseline as MutBpos
meta_NoMut_mutB_with_baseline <- data_NoMut_mutB_with_baseline %>%
  select(SampleName, Group) %>%
  distinct() %>%
  mutate(disease_type = ifelse(Group == "MutB-Baseline", "MutBpos", Group)) %>%
  rename(sampleName_long = SampleName) %>%
  mutate(disease_type = factor(disease_type, levels = c("NoMut", "MutBpos")))

# DE Analysis
lm_NoMut_mutB_with_baseline <- lmNULISAseq(
  data = as.matrix(data_NoMut_mutB_wide_with_baseline),
  sampleInfo = meta_NoMut_mutB_with_baseline,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)

NoMut_MutBpos_DE_proteins_with_baseline <- lm_NoMut_mutB_with_baseline$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target = target,
    Coef = X_coef,
    t_stat = X_tstat,
    Pval_FDR = X_pval_FDR,
    Pval_unadj = X_pval_unadj,
    Pval_bonf = X_pval_bonf
  )

write.csv(NoMut_MutBpos_DE_proteins_with_baseline, "NoMut_MutBpos_DE_proteins_with_MutBBaseline.csv", row.names = FALSE)

pthresh <- 0.05

# Annotate significance (FDR)
NoMut_MutBpos_DE_proteins_with_baseline$Significance_FDR <- with(
  NoMut_MutBpos_DE_proteins_with_baseline, 
  ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
         ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

NoMut_MutBpos_FDR_volcano_plot_with_baseline <- ggplot(NoMut_MutBpos_DE_proteins_with_baseline, aes(x = Coef, y = -log10(Pval_FDR), color = Significance_FDR, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(NoMut_MutBpos_DE_proteins_with_baseline, Significance_FDR != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("NoMut_MutBpos_FDR_volcano_plot_with_MutBBaseline.png", plot = NoMut_MutBpos_FDR_volcano_plot_with_baseline, width = 8, height = 6, dpi = 300)
print(NoMut_MutBpos_FDR_volcano_plot_with_baseline)

# Annotate significance (unadjusted)
NoMut_MutBpos_DE_proteins_with_baseline$Significance_unadj <- with(
  NoMut_MutBpos_DE_proteins_with_baseline, 
  ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
         ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

NoMut_MutBpos_unadj_volcano_plot_with_baseline <- ggplot(NoMut_MutBpos_DE_proteins_with_baseline, aes(x = Coef, y = -log10(Pval_unadj), color = Significance_unadj, label = Target)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "coral2", "Downregulated" = "#4169E1", "Not Significant" = "gray70")) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(data = subset(NoMut_MutBpos_DE_proteins_with_baseline, Significance_unadj != "Not Significant"),
                           size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC for DE Genes",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("NoMut_MutBpos_unadj_volcano_plot_with_MutBBaseline.png", plot = NoMut_MutBpos_unadj_volcano_plot_with_baseline, width = 8, height = 6, dpi = 300)
print(NoMut_MutBpos_unadj_volcano_plot_with_baseline)


library(dplyr)
Data$Group <- case_when(
  grepl("Control", Data$SampleName) ~ "Control",
  grepl("NoMut", Data$SampleName) ~ "Case",
  grepl("MutBpos", Data$SampleName) ~ "Case",
  grepl("MutB-Baseline", Data$SampleName) ~ "Case",
  TRUE ~ "Other"
)

# Filter for Control vs NoMut
data_Case <- Data %>% filter(Group %in% c("Control", "Case"))

# Reshape to wide format
data_Case_wide <- data_Case %>%
  select(SampleName, Target, NPQ) %>%
  pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  column_to_rownames("Target")

# Metadata
meta_case <- data_Case %>%
  select(SampleName, Group) %>%
  distinct() %>%
  rename(sampleName_long = SampleName, disease_type = Group) %>%
  mutate(disease_type = factor(disease_type, levels = c("Control", "Case")))

# DE Analysis
lm_case <- lmNULISAseq(
  data = as.matrix(data_Case_wide),
  sampleInfo = meta_case,
  sampleName_var = "sampleName_long",
  modelFormula = "disease_type"
)

Case_DE_proteins <- lm_case$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target      = target,
    Coef        = X_coef,
    t_stat      = X_tstat,
    Pval_FDR    = X_pval_FDR,
    Pval_unadj  = X_pval_unadj,
    Pval_bonf   = X_pval_bonf
  )

write.csv(Case_DE_proteins, "Case_DE_proteins.csv", row.names = FALSE)

# Threshold for significance
pthresh <- 0.05

##### Volcano plot – FDR-adjusted p-values #####
Case_DE_proteins$Significance <- with(
  Case_DE_proteins,
  ifelse(Pval_FDR < pthresh & Coef > 0, "Upregulated",
         ifelse(Pval_FDR < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

Case_FDR_volcano_plot <- ggplot(
  Case_DE_proteins,
  aes(x = Coef, y = -log10(Pval_FDR), color = Significance, label = Target)
) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated"   = "coral2",
    "Downregulated" = "#4169E1",
    "Not Significant" = "gray70"
  )) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed",
             color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(
    data = subset(Case_DE_proteins, Significance != "Not Significant"),
    size = 3, max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC (Case vs Control)",
    x = "Log2 Fold Change",
    y = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("Case_FDR_volcano_plot.png",
       plot = Case_FDR_volcano_plot, width = 8, height = 6, dpi = 300)

print(Case_FDR_volcano_plot)

##### Volcano plot – unadjusted p-values #####
Case_DE_proteins$Significance <- with(
  Case_DE_proteins,
  ifelse(Pval_unadj < pthresh & Coef > 0, "Upregulated",
         ifelse(Pval_unadj < pthresh & Coef < 0, "Downregulated", "Not Significant"))
)

Case_unadj_volcano_plot <- ggplot(
  Case_DE_proteins,
  aes(x = Coef, y = -log10(Pval_unadj), color = Significance, label = Target)
) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated"   = "coral2",
    "Downregulated" = "#4169E1",
    "Not Significant" = "gray70"
  )) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed",
             color = "gray28", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray28", size = 0.5) +
  ggrepel::geom_text_repel(
    data = subset(Case_DE_proteins, Significance != "Not Significant"),
    size = 3, max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC (Case vs Control)",
    x = "Log2 Fold Change",
    y = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5)
  )

ggsave("Case_unadj_volcano_plot.png",
       plot = Case_unadj_volcano_plot, width = 8, height = 6, dpi = 300)

print(Case_unadj_volcano_plot)


# Remove excluded_sample_2 (see config: exclude_samples)
Data_wo_noMut_pat <- subset(Data, SampleName != excluded_sample_2)

# Grouping
Data_wo_noMut_pat$Group_wo_noMut_pat <- case_when(
  grepl("Control",  Data_wo_noMut_pat$SampleName) ~ "Control",
  grepl("NoMut", Data_wo_noMut_pat$SampleName) ~ "Case",
  grepl("MutBpos",  Data_wo_noMut_pat$SampleName) ~ "Case",
  grepl("MutB-Baseline", Data_wo_noMut_pat$SampleName) ~ "Case",
  TRUE ~ "Other"
)

# Filter for Control vs Case
data_Case_wo_noMut_pat <- Data_wo_noMut_pat %>%
  dplyr::filter(Group_wo_noMut_pat %in% c("Control", "Case"))

# Reshape to wide format
data_Case_wide_wo_noMut_pat <- data_Case_wo_noMut_pat %>%
  dplyr::select(SampleName, Target, NPQ) %>%
  tidyr::pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  tibble::column_to_rownames("Target")

# Metadata
meta_case_wo_noMut_pat <- data_Case_wo_noMut_pat %>%
  dplyr::select(SampleName, Group_wo_noMut_pat) %>%
  dplyr::distinct() %>%
  dplyr::rename(
    sampleName_long_wo_noMut_pat = SampleName,
    disease_type_wo_noMut_pat    = Group_wo_noMut_pat
  ) %>%
  dplyr::mutate(
    disease_type_wo_noMut_pat = factor(
      disease_type_wo_noMut_pat,
      levels = c("Control", "Case")
    )
  )

# DE Analysis
lm_case_wo_noMut_pat <- lmNULISAseq(
  data          = as.matrix(data_Case_wide_wo_noMut_pat),
  sampleInfo    = meta_case_wo_noMut_pat,
  sampleName_var = "sampleName_long_wo_noMut_pat",
  modelFormula   = "disease_type_wo_noMut_pat"
)

Case_DE_proteins_wo_noMut_pat <- lm_case_wo_noMut_pat$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target_wo_noMut_pat     = target,
    Coef_wo_noMut_pat       = X_coef,
    t_stat_wo_noMut_pat     = X_tstat,
    Pval_FDR_wo_noMut_pat   = X_pval_FDR,
    Pval_unadj_wo_noMut_pat = X_pval_unadj,
    Pval_bonf_wo_noMut_pat  = X_pval_bonf
  )

write.csv(
  Case_DE_proteins_wo_noMut_pat,
  "Case_DE_proteins_wo_noMut_pat.csv",
  row.names = FALSE
)

# Threshold for significance
pthresh_wo_noMut_pat <- 0.05

##### Volcano plot – FDR-adjusted p-values #####
Case_DE_proteins_wo_noMut_pat$Significance_wo_noMut_pat <- with(
  Case_DE_proteins_wo_noMut_pat,
  ifelse(
    Pval_FDR_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat > 0,
    "Upregulated",
    ifelse(
      Pval_FDR_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat < 0,
      "Downregulated",
      "Not Significant"
    )
  )
)

Case_FDR_volcano_plot_wo_noMut_pat <- ggplot(
  Case_DE_proteins_wo_noMut_pat,
  aes(
    x     = Coef_wo_noMut_pat,
    y     = -log10(Pval_FDR_wo_noMut_pat),
    color = Significance_wo_noMut_pat,
    label = Target_wo_noMut_pat
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated"    = "coral2",
    "Downregulated"  = "#4169E1",
    "Not Significant" = "gray70"
  )) +
  geom_hline(
    yintercept = -log10(pthresh_wo_noMut_pat),
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = subset(
      Case_DE_proteins_wo_noMut_pat,
      Significance_wo_noMut_pat != "Not Significant"
    ),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot (FDR): Case vs Control (wo noMut patient)",
    x     = "Log2 Fold Change",
    y     = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black", size = 0.5)
  )

ggsave(
  "Case_FDR_volcano_plot_wo_noMut_pat.png",
  plot   = Case_FDR_volcano_plot_wo_noMut_pat,
  width  = 8,
  height = 6,
  dpi    = 300
)

print(Case_FDR_volcano_plot_wo_noMut_pat)

##### Volcano plot – unadjusted p-values #####
Case_DE_proteins_wo_noMut_pat$Significance_wo_noMut_pat <- with(
  Case_DE_proteins_wo_noMut_pat,
  ifelse(
    Pval_unadj_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat > 0,
    "Upregulated",
    ifelse(
      Pval_unadj_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat < 0,
      "Downregulated",
      "Not Significant"
    )
  )
)

Case_unadj_volcano_plot_wo_noMut_pat <- ggplot(
  Case_DE_proteins_wo_noMut_pat,
  aes(
    x     = Coef_wo_noMut_pat,
    y     = -log10(Pval_unadj_wo_noMut_pat),
    color = Significance_wo_noMut_pat,
    label = Target_wo_noMut_pat
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated"    = "coral2",
    "Downregulated"  = "#4169E1",
    "Not Significant" = "gray70"
  )) +
  geom_hline(
    yintercept = -log10(pthresh_wo_noMut_pat),
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = subset(
      Case_DE_proteins_wo_noMut_pat,
      Significance_wo_noMut_pat != "Not Significant"
    ),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot (unadjusted): Case vs Control (wo noMut patient)",
    x     = "Log2 Fold Change",
    y     = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black", size = 0.5)
  )

ggsave(
  "Case_unadj_volcano_plot_wo_noMut_pat.png",
  plot   = Case_unadj_volcano_plot_wo_noMut_pat,
  width  = 8,
  height = 6,
  dpi    = 300
)

print(Case_unadj_volcano_plot_wo_noMut_pat)

# Define your grouping based on SampleName — adjust pattern as needed
Data_wo_noMut_pat$Group_wo_noMut_pat <- case_when(
  grepl("Control",         Data_wo_noMut_pat$SampleName) ~ "Control",
  grepl("NoMut",        Data_wo_noMut_pat$SampleName) ~ "NoMut",
  grepl("MutC",             Data_wo_noMut_pat$SampleName) ~ "MutC",
  grepl("MutBpos",         Data_wo_noMut_pat$SampleName) ~ "MutBpos",
  grepl("DrugA",        Data_wo_noMut_pat$SampleName) ~ "DrugA",
  grepl("DrugB-1month", Data_wo_noMut_pat$SampleName) ~ "DrugB-1month",
  grepl("DrugB-2month", Data_wo_noMut_pat$SampleName) ~ "DrugB-2month",
  grepl("DrugB-3month", Data_wo_noMut_pat$SampleName) ~ "DrugB-3month",
  grepl("DrugB-4month", Data_wo_noMut_pat$SampleName) ~ "DrugB-4month",
  grepl("DrugB-5month", Data_wo_noMut_pat$SampleName) ~ "DrugB-5month",
  grepl("MutB-Baseline",     Data_wo_noMut_pat$SampleName) ~ "MutB-Baseline",
  TRUE ~ "Other"
)

# Filter for Control vs NoMut
data_cs_wo_noMut_pat <- Data_wo_noMut_pat %>%
  dplyr::filter(Group_wo_noMut_pat %in% c("Control", "NoMut"))

# Reshape to wide format
data_cs_wide_wo_noMut_pat <- data_cs_wo_noMut_pat %>%
  dplyr::select(SampleName, Target, NPQ) %>%
  tidyr::pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  tibble::column_to_rownames("Target")

# Metadata
meta_cs_wo_noMut_pat <- data_cs_wo_noMut_pat %>%
  dplyr::select(SampleName, Group_wo_noMut_pat) %>%
  dplyr::distinct() %>%
  dplyr::rename(
    sampleName_long_wo_noMut_pat = SampleName,
    disease_type_wo_noMut_pat    = Group_wo_noMut_pat
  )

# DE Analysis
lm_cs_wo_noMut_pat <- lmNULISAseq(
  data           = as.matrix(data_cs_wide_wo_noMut_pat),
  sampleInfo     = meta_cs_wo_noMut_pat,
  sampleName_var = "sampleName_long_wo_noMut_pat",
  modelFormula   = "disease_type_wo_noMut_pat"
)

NoMut_DE_proteins_wo_noMut_pat <- lm_cs_wo_noMut_pat$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target_wo_noMut_pat     = target,
    Coef_wo_noMut_pat       = X_coef,
    t_stat_wo_noMut_pat     = X_tstat,
    Pval_FDR_wo_noMut_pat   = X_pval_FDR,
    Pval_unadj_wo_noMut_pat = X_pval_unadj,
    Pval_bonf_wo_noMut_pat  = X_pval_bonf
  )

write.csv(
  NoMut_DE_proteins_wo_noMut_pat,
  "NoMut_DE_proteins_wo_noMut_pat.csv",
  row.names = FALSE
)

# Threshold for significance (you can adjust this)
pthresh_wo_noMut_pat <- 0.05

# Annotate significance (FDR)
NoMut_DE_proteins_wo_noMut_pat$Significance_wo_noMut_pat <- with(
  NoMut_DE_proteins_wo_noMut_pat,
  ifelse(
    Pval_FDR_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat > 0,
    "Upregulated",
    ifelse(
      Pval_FDR_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat < 0,
      "Downregulated",
      "Not Significant"
    )
  )
)

NoMut_FDR_volcano_plot_wo_noMut_pat <- ggplot(
  NoMut_DE_proteins_wo_noMut_pat,
  aes(
    x     = Coef_wo_noMut_pat,
    y     = -log10(Pval_FDR_wo_noMut_pat),
    color = Significance_wo_noMut_pat,
    label = Target_wo_noMut_pat
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated"    = "coral2",
    "Downregulated"  = "#4169E1",
    "Not Significant" = "gray70"
  )) +
  geom_hline(
    yintercept = -log10(pthresh_wo_noMut_pat),
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = subset(
      NoMut_DE_proteins_wo_noMut_pat,
      Significance_wo_noMut_pat != "Not Significant"
    ),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: FDR-Adjusted P vs logFC (Control vs NoMut, wo noMut patient)",
    x     = "Log2 Fold Change",
    y     = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black", size = 0.5)
  )

ggsave(
  "NoMut_FDR_volcano_plot_wo_noMut_pat.png",
  plot   = NoMut_FDR_volcano_plot_wo_noMut_pat,
  width  = 8,
  height = 6,
  dpi    = 300
)

print(NoMut_FDR_volcano_plot_wo_noMut_pat)

# Annotate significance (unadjusted)
NoMut_DE_proteins_wo_noMut_pat$Significance_wo_noMut_pat <- with(
  NoMut_DE_proteins_wo_noMut_pat,
  ifelse(
    Pval_unadj_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat > 0,
    "Upregulated",
    ifelse(
      Pval_unadj_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat < 0,
      "Downregulated",
      "Not Significant"
    )
  )
)

NoMut_unadj_volcano_plot_wo_noMut_pat <- ggplot(
  NoMut_DE_proteins_wo_noMut_pat,
  aes(
    x     = Coef_wo_noMut_pat,
    y     = -log10(Pval_unadj_wo_noMut_pat),
    color = Significance_wo_noMut_pat,
    label = Target_wo_noMut_pat
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated"    = "coral2",
    "Downregulated"  = "#4169E1",
    "Not Significant" = "gray70"
  )) +
  geom_hline(
    yintercept = -log10(pthresh_wo_noMut_pat),
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = subset(
      NoMut_DE_proteins_wo_noMut_pat,
      Significance_wo_noMut_pat != "Not Significant"
    ),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Unadjusted P vs logFC (Control vs NoMut, wo noMut patient)",
    x     = "Log2 Fold Change",
    y     = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black", size = 0.5)
  )

ggsave(
  "NoMut_unadj_volcano_plot_wo_noMut_pat.png",
  plot   = NoMut_unadj_volcano_plot_wo_noMut_pat,
  width  = 8,
  height = 6,
  dpi    = 300
)

print(NoMut_unadj_volcano_plot_wo_noMut_pat)


# Filter including MutB-Baseline, using Data_wo_noMut_pat
data_NoMut_mutB_with_baseline_wo_noMut_pat <- Data_wo_noMut_pat %>%
  dplyr::filter(Group_wo_noMut_pat %in% c("NoMut", "MutBpos", "MutB-Baseline")) %>%
  dplyr::filter(SampleName != excluded_sample_1) %>%
  dplyr::filter(!grepl("DrugB", SampleName, ignore.case = TRUE))

# Reshape to wide format
data_NoMut_mutB_wide_with_baseline_wo_noMut_pat <- data_NoMut_mutB_with_baseline_wo_noMut_pat %>%
  dplyr::select(SampleName, Target, NPQ) %>%
  tidyr::pivot_wider(names_from = SampleName, values_from = NPQ) %>%
  tibble::column_to_rownames("Target")

# Metadata: treat MutB-Baseline as MutBpos
meta_NoMut_mutB_with_baseline_wo_noMut_pat <- data_NoMut_mutB_with_baseline_wo_noMut_pat %>%
  dplyr::select(SampleName, Group_wo_noMut_pat) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    disease_type_wo_noMut_pat = ifelse(
      Group_wo_noMut_pat == "MutB-Baseline", "MutBpos", Group_wo_noMut_pat
    )
  ) %>%
  dplyr::rename(sampleName_long_wo_noMut_pat = SampleName) %>%
  dplyr::mutate(
    disease_type_wo_noMut_pat = factor(
      disease_type_wo_noMut_pat,
      levels = c("NoMut", "MutBpos")
    )
  )

# DE Analysis
lm_NoMut_mutB_with_baseline_wo_noMut_pat <- lmNULISAseq(
  data           = as.matrix(data_NoMut_mutB_wide_with_baseline_wo_noMut_pat),
  sampleInfo     = meta_NoMut_mutB_with_baseline_wo_noMut_pat,
  sampleName_var = "sampleName_long_wo_noMut_pat",
  modelFormula   = "disease_type_wo_noMut_pat"
)

NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat <-
  lm_NoMut_mutB_with_baseline_wo_noMut_pat$modelStats %>%
  dplyr::select(target, X_coef, X_tstat, X_pval_FDR, X_pval_unadj, X_pval_bonf) %>%
  dplyr::rename(
    Target_wo_noMut_pat     = target,
    Coef_wo_noMut_pat       = X_coef,
    t_stat_wo_noMut_pat     = X_tstat,
    Pval_FDR_wo_noMut_pat   = X_pval_FDR,
    Pval_unadj_wo_noMut_pat = X_pval_unadj,
    Pval_bonf_wo_noMut_pat  = X_pval_bonf
  )

write.csv(
  NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat,
  "NoMut_MutBpos_DE_proteins_with_MutBBaseline_wo_noMut_pat.csv",
  row.names = FALSE
)

pthresh_wo_noMut_pat <- 0.05

# Annotate significance (FDR)
NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat$Significance_FDR_wo_noMut_pat <- with(
  NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat,
  ifelse(
    Pval_FDR_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat > 0,
    "Upregulated",
    ifelse(
      Pval_FDR_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat < 0,
      "Downregulated",
      "Not Significant"
    )
  )
)

NoMut_MutBpos_FDR_volcano_plot_with_baseline_wo_noMut_pat <- ggplot(
  NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat,
  aes(
    x     = Coef_wo_noMut_pat,
    y     = -log10(Pval_FDR_wo_noMut_pat),
    color = Significance_FDR_wo_noMut_pat,
    label = Target_wo_noMut_pat
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated"    = "coral2",
    "Downregulated"  = "#4169E1",
    "Not Significant" = "gray70"
  )) +
  geom_hline(
    yintercept = -log10(pthresh_wo_noMut_pat),
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = subset(
      NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat,
      Significance_FDR_wo_noMut_pat != "Not Significant"
    ),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot (FDR): NoMut vs MutBpos (with MutB baseline, wo noMut pat)",
    x     = "Log2 Fold Change",
    y     = "-Log10 (FDR-Adjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black", size = 0.5)
  )

ggsave(
  "NoMut_MutBpos_FDR_volcano_plot_with_MutBBaseline_wo_noMut_pat.png",
  plot   = NoMut_MutBpos_FDR_volcano_plot_with_baseline_wo_noMut_pat,
  width  = 8,
  height = 6,
  dpi    = 300
)
print(NoMut_MutBpos_FDR_volcano_plot_with_baseline_wo_noMut_pat)

# Annotate significance (unadjusted)
NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat$Significance_unadj_wo_noMut_pat <- with(
  NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat,
  ifelse(
    Pval_unadj_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat > 0,
    "Upregulated",
    ifelse(
      Pval_unadj_wo_noMut_pat < pthresh_wo_noMut_pat & Coef_wo_noMut_pat < 0,
      "Downregulated",
      "Not Significant"
    )
  )
)

NoMut_MutBpos_unadj_volcano_plot_with_baseline_wo_noMut_pat <- ggplot(
  NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat,
  aes(
    x     = Coef_wo_noMut_pat,
    y     = -log10(Pval_unadj_wo_noMut_pat),
    color = Significance_unadj_wo_noMut_pat,
    label = Target_wo_noMut_pat
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated"    = "coral2",
    "Downregulated"  = "#4169E1",
    "Not Significant" = "gray70"
  )) +
  geom_hline(
    yintercept = -log10(pthresh_wo_noMut_pat),
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "gray28",
    size       = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = subset(
      NoMut_MutBpos_DE_proteins_with_baseline_wo_noMut_pat,
      Significance_unadj_wo_noMut_pat != "Not Significant"
    ),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot (unadj): NoMut vs MutBpos (with MutB baseline, wo noMut pat)",
    x     = "Log2 Fold Change",
    y     = "-Log10 (Unadjusted p-value)"
  ) +
  theme(
    legend.title = element_blank(),
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black", size = 0.5)
  )

ggsave(
  "NoMut_MutBpos_unadj_volcano_plot_with_MutBBaseline_wo_noMut_pat.png",
  plot   = NoMut_MutBpos_unadj_volcano_plot_with_baseline_wo_noMut_pat,
  width  = 8,
  height = 6,
  dpi    = 300
)
print(NoMut_MutBpos_unadj_volcano_plot_with_baseline_wo_noMut_pat)
