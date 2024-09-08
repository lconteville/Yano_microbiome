library("phyloseq")
library("ggplot2")
library("reshape2")
library("ggpubr")
library("cowplot")
library("stringr")

setwd("D:/Github/Yano_microbiome/")

#### SET colors and theme ####
group.colors <- c("Yano_BR" = "brown2", "Yano_VE" = "#E7B800",
                  "Matses" ="#53B400", "Tunapuco" = "#4E84C4",
                  "US" = "darkorchid1")

groups.order= c("Yano_BR","Yano_VE","Matses","Tunapuco","US")

my_comparisons <- list( c("Yano_BR", "Yano_VE"), c("Yano_BR", "Matses"), c("Yano_BR", "Tunapuco"), c("Yano_BR", "US") )

theme_set(theme_classic()+ theme(
  axis.text = element_text(size=8, colour="black"),
  axis.title = element_text(size=8, colour="black"),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  strip.text = element_text(size = 8, colour="black"),
  plot.background = element_blank(),
  plot.title = element_text(size=8,colour="black")
))

#### INPUT DATA LEVEL 3 ####
otulevel3 <- read.csv("data/level3_relab.txt", header=TRUE, row.names = 1, sep="\t")
taxlevel3 <- read.csv("data/level3_tax.txt", header=TRUE, sep="\t")
metadata <- read.csv("data/metadata.tsv", header=TRUE, row.names=1, sep="\t")

rownames(taxlevel3) <- taxlevel3$Level_3
taxlevel3matrix =  as.matrix(taxlevel3)

# Ensure OTU and tax tables have matching row names
identical(rownames(otulevel3), rownames(taxlevel3))

# Create phyloseq object for Level3
OTU3 = otu_table(otulevel3, taxa_are_rows = TRUE)
TAX3 = tax_table(taxlevel3matrix)
sampledata = sample_data(metadata)
func_level3 = phyloseq(OTU3, TAX3, sampledata)

sample_data(func_level3)$Group <- factor((sample_data(func_level3)$Group),
                                         levels=groups.order)

#### Ordination Plot ####

# Perform Principal Coordinate Analysis (PCoA) using Bray-Curtis dissimilarity
PCoA_level3 <- ordinate(func_level3, "PCoA", "bray")

# Invert the first axis (PC1) for visualization
coords <- PCoA_level3$vectors
coords[, 1] <- -coords[, 1]
PCoA_level3$vectors <- coords

# Plot PCoA
fig4a <- plot_ordination(func_level3, PCoA_level3, color = "Group")+ 
  geom_point(size = 2.8, alpha = 0.9) + 
  scale_color_manual(values = group.colors) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.28, "cm"),
        legend.background = element_rect(color = "black", fill = "white",
                                         linewidth = 0.25, linetype = "solid"))

fig4a

#### INPUT DATA LEVEL 1 ####
otulevel1 <- read.csv("data/level1_relab.tsv", header=TRUE, row.names=1, sep="\t")
taxlevel1 <- data.frame(row.names = rownames(otulevel1), Level1 = rownames(otulevel1))
taxlevel1matrix =  as.matrix(taxlevel1)

# Create phyloseq object for Level3
OTU1 = otu_table(otulevel1, taxa_are_rows = TRUE)
TAX1 = tax_table(taxlevel1matrix)
func_level1 = phyloseq(OTU1, TAX1, sampledata)

sample_data(func_level1)$Group <- factor((sample_data(func_level1)$Group),
                                         levels=groups.order)

#### BOXPLOT MAIN FUNCS LEVEL 1 ####
main.funcs = c("Carbohydrates","Protein_Metabolism","Membrane_Transport")
mainfunc_level1 = prune_taxa((tax_table(func_level1)[, "Level1"] %in% main.funcs), func_level1)
mainfunc_level1df <- psmelt(mainfunc_level1)

mainfunc_level1df$Level1 <- factor(mainfunc_level1df$Level1,
                                   levels=c("Carbohydrates","Protein_Metabolism","Membrane_Transport"))

# Create boxplot of relative abundance for main functions
fig4b <- ggplot(mainfunc_level1df,
                aes(x = Group, y = Abundance, fill = Level1)) + 
  geom_boxplot(width = 0.9, outlier.shape = NA, 
               alpha = 1, position = position_dodge(1)) +
  ylab("Relative Abundance (%)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        legend.direction = "horizontal", 
        legend.position = "top",
        legend.justification = 'center',
        legend.title = element_blank(),
        legend.spacing.x = unit(.255, 'cm'),
        legend.background = element_rect(color = "black",
                                         fill = "white",
                                         linewidth = 0.25,
                                         linetype = "solid")) +
  scale_fill_manual(values = c("#FF8847", "deepskyblue3", "chartreuse3" ),
                    labels = c("Carbohydrates\nMetabolism", "Protein\nMetabolism", "Membrane\nTransport")) 

fig4b

#### INPUT DATA LEVEL 2 ####
otulevel2 <- read.csv("data/level2_relab.tsv", header=TRUE, row.names=1, sep="\t")
taxlevel2 <- read.csv("data/level2_tax.txt", header=TRUE, sep="\t")
rownames(taxlevel2) <- taxlevel2$Level2
taxlevel2matrix =  as.matrix(taxlevel2)

# Create phyloseq object for Level2
OTU2 = otu_table(otulevel2, taxa_are_rows = TRUE)
TAX2 = tax_table(taxlevel2matrix)
func_level2 = phyloseq(OTU2, TAX2, sampledata)

sample_data(func_level2)$Group <- factor((sample_data(func_level2)$Group),
                                         levels=groups.order)

#### BARPLOT CARBOHYDRATES METABOLISM ####
func_level2_carbo <- subset_taxa(func_level2, Level1=="Carbohydrates")

unique(tax_table(func_level2_carbo)[,"Level2"])

# Merge samples by group and calculate relative abundances
mergedlevel2_carbo = merge_samples(func_level2_carbo, "Group")
mergedlevel2_carbo2 = transform_sample_counts(mergedlevel2_carbo,
                                              function(x) x / sum(x) * 100)

# Collapse low abundance taxa
taxa_abund <- taxa_sums(mergedlevel2_carbo2) / sum(taxa_sums(mergedlevel2_carbo2))
low_abund_taxa <- names(taxa_abund[taxa_abund < 0.03])
tax_table(mergedlevel2_carbo2)[low_abund_taxa, "Level2"] <- "Others"

unique(tax_table(mergedlevel2_carbo2)[,"Level2"])

# Collapse to Level2 
physeq_collapsed <- tax_glom(mergedlevel2_carbo2, "Level2")
df_mergedlevel2_carbo <- psmelt(physeq_collapsed)

# Clean up function labels
df_mergedlevel2_carbo$Level2 <- str_replace_all(df_mergedlevel2_carbo$Level2,"_metabolism","")
df_mergedlevel2_carbo$Level2 <- str_replace_all(df_mergedlevel2_carbo$Level2,"_Metabolism","")
df_mergedlevel2_carbo$Level2 <- str_replace_all(df_mergedlevel2_carbo$Level2,"-_and_","/")
df_mergedlevel2_carbo$Level2 <- str_replace_all(df_mergedlevel2_carbo$Level2,"_"," ")

df_mergedlevel2_carbo$Level2 <- factor(df_mergedlevel2_carbo$Level2, levels=c("Others", "Aminosugars", "One-carbon",
                                                                              "Nucleotide sugars", "Fermentation", "Polysaccharides",
                                                                              "Di/oligosaccharides","Monosaccharides",
                                                                              "Central carbohydrate"))

df_mergedlevel2_carbo$Group = factor(df_mergedlevel2_carbo$Sample,
                                 levels = groups.order)

# Create stacked barplot with Level2 of Carbohydrates Metabolism
fig4c <- ggplot(df_mergedlevel2_carbo, aes(x=Group, y=Abundance, fill=Level2)) +
  geom_bar(stat="identity", position= "stack", alpha = 1, width = 0.7) +
  ylab("Relative Abundance (%)") +
  scale_fill_manual(values = c("#7DA1B9", "#80C2EC", "#D2EBF2", "#D7D4D4", "#FFD2AD", "#F2B380", 
                               "#DB717C", "#BB605D")) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.spacing.y = unit(0.04, "cm"),
        legend.key.spacing.x = unit(0.15, "cm"),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(margin = margin(l=2))
        ) +
  guides(fill = guide_legend(reverse = TRUE, nrow=3)) +
  scale_x_discrete(position = "top", expand = c(.1,0))  + 
  scale_y_continuous(expand = c(0,0)) + 
  geom_text(aes(label = round(Abundance,digits=1)), size = 3,
            position = position_stack(vjust = 0.5))  

fig4c

#### BOXPLOT LEVEL 3 FUNCS ####

func_level3_new <- prune_taxa(taxa_sums(func_level3) > 0, func_level3) 

func_level3_new <- subset_taxa(func_level3_new, Level_3 ==  "Pyruvate_metabolism_I:_anaplerotic_reactions_PEP" | Level_3 == "Glycolysis_and_Gluconeogenesis" | 
                          Level_3 == "D-Galacturonate_and_D-Glucuronate_Utilization" | Level_3 == "Xylose_utilization" | Level_3 == "Mannose_Metabolism"  |
                          Level_3 == "Lactose_utilization")

df_func_level3 <- psmelt(func_level3_new)

labs <- c("Pyruvate\n","Glycolysis/\nGluconeogenesis", "D-Galacturonate/\nD-Glucuronate",
          "Xylose\n","Mannose\n","Lactose\n")
names(labs) <- c("Pyruvate_metabolism_I:_anaplerotic_reactions_PEP","Glycolysis_and_Gluconeogenesis", "D-Galacturonate_and_D-Glucuronate_Utilization",
                 "Xylose_utilization","Mannose_Metabolism", "Lactose_utilization")
df_func_level3$Level_3 <- factor(df_func_level3$Level_3, levels = names(labs))

# Create boxplot with Level3 functions
fig4d <- ggplot(df_func_level3, aes(x=V2, y=Abundance, fill=V2)) + 
  facet_wrap( . ~ Level_3, ncol = 3, labeller = labeller(Level_3 = labs)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + 
  scale_fill_manual(values = c("green", "darkorchid1")) +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black",  linewidth = .55),
        strip.background = element_rect(linewidth = .55),
        panel.grid.major = element_line(),
        panel.grid.minor = element_line()
        ) +
  scale_x_discrete(labels=c("T", "U")) +
  scale_y_continuous(limits = c(0,1.7)) +
  ylab("Relative Abundance (%)")

fig4d

#### BOXPLOT YANO_BR UNIQUE FUNCS ####

# Convert phyloseq to data frame for plotting
df_func_level1 <- psmelt(func_level1)

level1_list <- unique(df_func_level1$Level1)

# Loop through each unique function and create boxplots
plot_list = list()
for (funcao in level1_list) {
  p <- ggplot(df_func_level1[df_func_level1$Level1 == funcao, ],
              aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    scale_fill_manual(values = group.colors) +
    ylab("Relative Abundance (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  plot_list[[funcao]] = p
}

Regulation <- plot_list[["Regulation_and_Cell_signaling"]] +
              ggtitle("Regulation/\nCell Signaling") +
              theme(plot.title = element_text(hjust=0.5))

Motility <-  plot_list[["Motility_and_Chemotaxis"]] +
             ggtitle("Motility/\nChemotaxis") +
             theme(plot.title = element_text(hjust=0.5),
                   axis.title.y=element_blank())

Virulence <- plot_list[["Virulence"]] +
             ggtitle("Virulence\n") +
             theme(plot.title = element_text(hjust=0.5),
                   axis.title.y=element_blank())

# Arrange the three plots in a single row
fig4e = ggarrange(Regulation, Motility, Virulence, ncol = 3, nrow = 1, align = "h")
fig4e

#### GENERATE PANEL WITH THE FIVE PLOTS ####

fig4 <- ggdraw() +
  draw_plot(fig4a, x = .0061, y = 0.661, width = 0.455, height = 0.33) +
  draw_plot(fig4b, x = 0.49, y = 0.639, width = 0.49, height = 0.35) +
  draw_plot(fig4c, x = .027, y = 0, width = 0.43, height = 0.65) +
  draw_plot(fig4d, x = 0.49, y = 0.275, width = 0.49, height = 0.37) +
  draw_plot(fig4e, x = 0.49, y = 0, width = 0.5, height = 0.27) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 11,
                  x = c(0, .5, 0, .5,.5), y = c(1, 1, .64, .64, .27))+
  theme(panel.background = element_rect(fill = "white",colour = "white"))

fig4

ggsave("Yano_fig4/figure4.png", fig4, width = 2300, height = 2300, dpi = 300, units = "px")


