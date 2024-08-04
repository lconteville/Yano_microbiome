library("phyloseq")
library("ggplot2")
library("ggpubr")
library("tidyverse")
library("cowplot")

setwd("D:/Github/Yano_microbiome/")

#### INPUT DATA ####
otumat <- read.table("data/kraken_bact_phyla.tsv", header=TRUE, row.names=1, sep="\t")
taxmat <- read.table("data/kraken_bact_phyla_tax.tsv", header=TRUE, sep="\t")
metadata <- read.table("data/metadata.tsv", header=TRUE, row.names=1, sep="\t")

# Ensure row names of OTU table match Phyla in taxmat
identical(rownames(otumat), taxmat$Phylum)

#### SET colors and theme ####
group.colors <- c("Yano_BR" = "brown2", "Yano_VE" = "#E7B800",
                  "Matses" ="#53B400", "Tunapuco" = "#4E84C4",
                  "US" = "darkorchid1")

phyla.colors <- c("Actinobacteria" = "#FF8847", "Bacteroidetes" = "#4FC79E",
                  "Firmicutes" = "#FFD651", "Proteobacteria" = "#1EA0CF",
                  "Spirochaetes" = "#B81A30", "Verrucomicrobia" = "#8491B4B2")

my_comparisons <- list( c("Yano_BR", "Yano_VE"), c("Yano_BR", "Matses"), c("Yano_BR", "Tunapuco"), c("Yano_BR", "US"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

theme_set(theme_classic()+ theme(
  axis.text = element_text(size=8, colour="black"),
  axis.title = element_text(size=8, colour="black"),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  strip.text.x = element_text(size = 8),
  plot.background = element_blank()
))

#### CREATE PHYLOSEQ OBJECT ####
rownames(taxmat) <- taxmat$Phylum
taxmatrix =  as.matrix(taxmat)

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmatrix)
sampledata = sample_data(metadata)

physeq = merge_phyloseq(OTU, TAX, sampledata)

sample_data(physeq)$Group <- factor((sample_data(physeq)$Group), levels=c("Yano_BR","Yano_VE","Matses","Tunapuco","US"))

#### GET THE TOP 6 PHYLA ####
phylum.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Phylum"], sum, na.rm=TRUE)
top6phyla = names(sort(phylum.sum, TRUE))[1:6]
physeq_top6 = prune_taxa((tax_table(physeq)[, "Phylum"] %in% top6phyla), physeq)

#### BARPLOT TOP 6 PHYLA ####
fig3A <- plot_bar(physeq_top6, fill = "Phylum") +
  facet_grid(. ~Group, drop=TRUE,scales = "free_x",  space = "free_x") +
  ylab("Relative Abundance (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        #legend.key.size = unit(.9, 'lines'),
        strip.background = element_blank()) +
  scale_fill_manual(values = phyla.colors) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0))
fig3A


#### BOXPLOT FIRMICUTES/BACTEROIDETES ####
df_top6 <- psmelt(physeq_top6)
bac_fir <- df_top6[df_top6$OTU %in% c("Bacteroidetes", "Firmicutes"),]

fig3B <- ggplot(data=bac_fir) + 
  geom_boxplot(aes(x=Group, y=Abundance, fill=Phylum), outlier.shape=NA, alpha=1) +
  ylab("Relative Abundance (%)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.position = "inside",
        legend.title = element_blank(),
        legend.key.size = unit(.8, 'lines'),
        legend.position.inside = c(0.25, 0.9),
        legend.spacing.y = unit(.7, 'cm'), 
        legend.background = element_rect(color = "black", fill = "white",
                                         linewidth = 0.25, linetype = "solid")) +
  scale_fill_manual(values = phyla.colors)

fig3B

#### BARPLOT LEFSE ####
lefse <- read.csv("data/lefse_bact_gen.tsv", sep="\t")
lefse$Group <- factor(lefse$Group, levels=c("Yano_BR","Yano_VE","Matses","Tunapuco","US"))
genera.order <- lefse$Genera[order(factor(lefse$Group), decreasing =TRUE)]
lefse$Genera <- factor(lefse$Genera, levels=genera.order)

fig3C <- ggplot(lefse, aes(x=LDAScore, y=Genera, fill=Group)) + 
  geom_bar(stat = "identity") + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 8.4),
        plot.margin = unit(c(0,1.2,0,0),"cm"),
        legend.position = "none") +
  scale_fill_manual(values= group.colors) +
  expand_limits(x = c(0,6)) +
  coord_cartesian(clip = 'off')+
  annotate("text", y = c(20.7,15.15,10.55,5.9,1.7), x = 6.2,  size = 3,
           label = c("Yano_BR","Yano_VE","Matses","Tunapuco","US    "),
           hjust = 0, parse = TRUE) +
  annotate("rect", ymin = c(.5,3.55,8.5,12.7,17.7), ymax = c(3.45,8.4,12.6,17.6,23.7), xmin = 5.9, xmax = 6,
           alpha = 1, fill = c( "darkorchid1","#4E84C4","#53B400", "#E7B800","brown2")) 

fig3C

#### BARPLOT ARCHAEA ####
otu_arch <- read.csv("data/kraken_arch_gen.tsv", header=TRUE, row.names=1, sep="\t")
tax_arch <- read.csv("data/kraken_arch_tax.tsv", header=TRUE, sep="\t")
metadata <- read.csv("data/metadata.tsv", header=TRUE, row.names=1, sep="\t")

### CREATE PHYLOSEQ OBJECT ###
rownames(tax_arch) <- tax_arch$ID
taxmatrix =  as.matrix(tax_arch)

OTU = otu_table(otu_arch, taxa_are_rows = TRUE)
TAX = tax_table(taxmatrix)
sampledata = sample_data(metadata)

physeq = phyloseq(OTU, TAX, sampledata)

sample_data(physeq)$Group <- factor((sample_data(physeq)$Group), levels=c("Yano_BR","Yano_VE","Matses","Tunapuco","CA","Norman","US"))

#### GET THE TOP 3 ARCHAEA ####
arch.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "ID"], sum, na.rm=TRUE)
top3arch = names(sort(arch.sum, TRUE))[1:3]
physeq_top3arch = prune_taxa((tax_table(physeq)[, "ID"] %in% top3arch), physeq)

df_top3arch <- psmelt(physeq_top3arch)

samples.order <- df_top3arch$Sample[order(df_top3arch$Group, df_top3arch$Sample, decreasing = FALSE)]
df_top3arch$Sample <- factor(df_top3arch$Sample, levels=unique(samples.order))

fig3D <- ggplot(df_top3arch, aes(x=Sample, y=Abundance, color=Group, fill=Group)) +
  geom_bar(stat="identity", position="stack") +
  facet_wrap(~ID) +
  ylab("Relative Abundance (%)") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = .3),
        axis.line.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black",  linewidth = .5),
        strip.background = element_rect(fill = "white", colour = "black")) + 
  grids(axis = "y", linetype = "dashed", color = "grey85") + 
  scale_y_continuous(limits = c(0,100), expand = c(0,0))

fig3D

#### GENERATE PANEL WITH THE FOUR PLOTS ####
fig3 <- ggdraw() +
  draw_plot(fig3A, x = .001, y = .68, width = 1, height = .33) +
  draw_plot(fig3B, x = .001, y = .33, width = .41, height = .35) +
  draw_plot(fig3C, x = .43, y = .356, width = .56, height = .32) +
  draw_plot(fig3D, x = .001, y = 0, width = .98, height = .33) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 11,
                  x = c(0, 0, .4, 0), y = c(1, .685, .685, .33))+
  theme(panel.background = element_rect(fill = "white",colour = "white"))

fig3

ggsave("Yano_fig3/figure3.png", fig3, width = 2300, height = 2700, dpi = 300, units = "px")

