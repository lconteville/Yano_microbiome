library("phyloseq")
library("ggplot2")
library("ggpubr")
library("cowplot")
library("reshape2")

setwd("D:/Github/Yano_microbiome/")

#### SET colors and theme ####
group.colors <- c("Yano_BR" = "brown2", "Yano_VE" = "#E7B800",
                  "Matses" ="#53B400", "Tunapuco" = "#4E84C4",
                  "US" = "darkorchid1")

my_comparisons <- list( c("Yano_BR", "Yano_VE"), c("Yano_BR", "Matses"), c("Yano_BR", "Tunapuco"), c("Yano_BR", "US"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

theme_set(theme_classic()+ theme(
  axis.text = element_text(size=9,colour="black"),
  axis.title = element_text(size=9,colour="black"),
  legend.text = element_text(size = 9)
))

#### INPUT DATA ####
otumat <- read.table("data/kraken_bact_gen.tsv", header=TRUE, row.names=1, sep="\t")
taxmat <- read.table("data/kraken_bact_gen_tax.tsv", header=TRUE, sep="\t")
metadata <- read.table("data/metadata.tsv", header=TRUE, row.names=1, sep="\t")

# Ensure row names of OTU table match Genus column in taxmat
identical(rownames(otumat), taxmat$Genus)

#### CREATE PHYLOSEQ OBJECT ####
rownames(taxmat) <- taxmat$Genus
taxmatrix =  as.matrix(taxmat)

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmatrix)
sampledata = sample_data(metadata)

physeq = merge_phyloseq(OTU, TAX, sampledata)

sample_data(physeq)$Group <- factor((sample_data(physeq)$Group), levels=c("Yano_BR","Yano_VE","Matses","Tunapuco","US"))

#### BOXPLOT OF ALPHA DIVERSITY ####
Shannon_gen <- estimate_richness(physeq, measures="Shannon")
Shannon_gen2 <- merge(Shannon_gen, sample_data(physeq), by='row.names')

fig2A <- ggplot(Shannon_gen2, aes(x=Group, y=Shannon, fill=Group)) + 
  geom_boxplot(alpha=1) + 
  scale_fill_manual(values=group.colors) +
  ylab("Shannon diversity") +
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif", symnum.args = symnum.args, size = 3.5)+
  theme(axis.title.x=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

fig2A 

#### DATA IN RELATIVE ABUNDANCE ####
Ab_rel_genero = transform_sample_counts(physeq, function(x) x / sum(x) * 100  )

#### BOXPLOT OF BETA DIVERSITY (Bray-Curtis) ####
# Calculate Bray-Curtis distance among samples
abrel_bray <- phyloseq::distance(Ab_rel_genero, method = "bray")

# Distance among samples of each group
abrel_bray_cadagrupo <- sapply(c("Y","SRR","SM","HC","SRS"),
                               function(letters) as.matrix(abrel_bray)[grep(letters, rownames(sample_data(Ab_rel_genero))), 
                                                                                 grep(letters, rownames(sample_data(Ab_rel_genero)))])
names(abrel_bray_cadagrupo) <- c("Yano_BR","Yano_VE","Matses", "Tunapuco","US")

# Extract the upper triangle of the distance matrices
extract_upper_tri <- function(mat) {
  return(mat[upper.tri(mat)])
}

upper_tri_vectors <- lapply(abrel_bray_cadagrupo, extract_upper_tri)
df.upper_tri_vectors = melt(upper_tri_vectors)

df.upper_tri_vectors$L1 <- factor(df.upper_tri_vectors$L1, levels=c("Yano_BR","Yano_VE","Matses","Tunapuco","US"))

fig2B <- ggplot(df.upper_tri_vectors, aes(x=L1, y=value, fill=L1)) + 
  geom_boxplot(alpha=1) + 
  scale_fill_manual(values=group.colors) +
  ylab("Bray-Curtis diversity") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", symnum.args = symnum.args, size = 3.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.position="none")

fig2B

#### Ordination Plot ####

pcoa_ord = ordinate(Ab_rel_genero, method = "PCoA", distance = abrel_bray)

fig2C <- plot_ordination(Ab_rel_genero, pcoa_ord, color = "Group") + 
  geom_point(size = 3, alpha=0.9) + 
  scale_color_manual(values=group.colors)+
  theme(legend.title = element_blank(), 
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.25),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.background = element_rect(color = "black", fill = "white", 
                                         linewidth = 0.25, linetype = "solid"))

fig2C


#### GENERATE PANEL WITH THE THREE PLOTS ####

fig2 <- ggdraw() +
  draw_plot(fig2A, x = .01, y = 0, width = .24, height = .98) +
  draw_plot(fig2B, x = .27, y = 0, width = .25, height = .98) +
  draw_plot(fig2C, x = .52, y = 0.055, width = .46, height = .925) +
  draw_plot_label(label = c("A", "B", "C"), size = 12,
                  x = c(0, .26, .51), y = c(1, 1, 1)) +
  theme(panel.background = element_rect(fill = "white",colour = "white"))

fig2

ggsave("Yano_fig2/figure2.png", fig2, width = 9, height = 3.5, dpi = 300, units = "in")
