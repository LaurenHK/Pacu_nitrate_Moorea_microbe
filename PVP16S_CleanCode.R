# Jordan Sims
# PVP 16S - Code for manuscript submission
# March 2021
# Condensing several R scripts with all analyses into one clean script with the analyses included in the manuscript text

library(phyloseq); packageVersion('phyloseq') #1.30.0
library(dplyr); packageVersion('dplyr') #1.0.3
library(vegan); packageVersion('vegan') #2.5.7
library(ggplot2); packageVersion('ggplot2') #3.3.3
library(tidyr); packageVersion('tidyr') #1.1.2

#---------Build and subset phyloseq object---------

# Read in the ASV table, taxonomy table, and sample metadata files
otu <- as.matrix(read.csv("ASVtable_PVP16S.csv", sep= ",", row.names = 1))
taxonomy <- as.matrix(read.csv("Taxonomy_PVP16S.csv", row.names = 1))
sampledat <- read.csv("SampleMetadata_PVP.csv", sep = ",", row.names = 1)

# Construct phyloseq object
PVP16S <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
                   sample_data(sampledat), 
                   tax_table(taxonomy))

sum(taxa_sums(PVP16S))
# 6683435 total reads

# Subset full phyloseq object into two based on sequencing facility
# GG includes parents and planulae sequenced at Georgia Genomics, OSU includes parents sequenced at OSU
GG <- subset_samples(PVP16S, Facility == "GG")
GG <- prune_taxa(taxa_sums(GG) > 0, GG)
OSU <- subset_samples(PVP16S, Facility == "OSU")
OSU <- prune_taxa(taxa_sums(OSU) > 0, OSU)

#---------OSU: Explore phyloseq object and filter reads---------

# Look at the distribution of sequence depth
hist(sample_sums(OSU), breaks=length(sample_sums(OSU)))

# Filter out unassigned reads and reads that identify as Mitochondria or Chloroplast DNA
OSU2 <- subset_taxa(OSU, (Kingdom != "NA"))
OSU3 <- subset_taxa(OSU2, (Family != "Mitochondria"))
OSU4 <- subset_taxa(OSU3, (Family != "Chloroplast"))

# Check the number of reads and ASVs remaining at each filtering step
sum(taxa_sums(OSU)); ncol(otu_table(OSU))
# 1903433; 1085 ASVs
sum(taxa_sums(OSU2)); ncol(otu_table(OSU2))
# 1903433 reads; 1085 ASVs
sum(taxa_sums(OSU3)); ncol(otu_table(OSU3))
# 1894562 reads; 876 ASVs
sum(taxa_sums(OSU4)); ncol(otu_table(OSU4))
# 1894562 reads; 876 ASVs

# Rename back to "OSU" for simplicity and remove transitional objects
OSU <- OSU4
rm(OSU2, OSU3, OSU4)

#---------OSU: Remove ASVs found in negative control---------
# Negative control was molecular water sent for sequencing, called "OSblank"

# Make a phyloseq object with just the negative sample and negative taxa
OSUblank <- subset_samples(OSU, sample_names(OSU) == "OSblank")
OSUblank <- prune_taxa(taxa_sums(OSUblank) > 0, OSUblank)

# Extract negative ASV list
OSUblank.ASVs <- rownames(data.frame(tax_table(OSUblank)))
# There are 210 unique sequences in OSblank

# Extract full list of ASVs from all OSU samples
OSU.ASVs <- rownames(data.frame(tax_table(OSU)))

# Make a list of the ASVs to keep by subtracting out the ones from the negative
OSUkeeptaxa <- OSU.ASVs[!(OSU.ASVs %in% OSUblank.ASVs)]

# Remove the taxa found in the negative control from the parent data
OSU <- prune_taxa(OSUkeeptaxa, OSU)
OSU <- prune_taxa(taxa_sums(OSU) > 0, OSU)

# Confirm there are no reads remaining for OSblank and remove from phyloseq object
sample_sums(OSU)
OSU <- subset_samples(OSU, SampleID != "Neg")

nrow(tax_table(OSU))
# 666 ASVs

#---------OSU: Compare diversity between pre- and post-planulation parents---------

# Subset phyloseq object into pre-planulation and post-planulation parents
OSUpre <- subset_samples(OSU, Stage == "ParentPre")
OSUpre <- prune_taxa(taxa_sums(OSUpre) > 0, OSUpre)
OSUpost <- subset_samples(OSU, Stage == "ParentPost")
OSUpost <- prune_taxa(taxa_sums(OSUpost) > 0, OSUpost)

# Compare ASV richness between pre-planulation and post-planulation parents
nrow(tax_table(OSUpre))
# 433 ASVs
nrow(tax_table(OSUpost))
# 263 ASVs

# Compare Family level richness between pre-planulation and post-planulation parents
# Overall Family level richness
OSUtax <- as.data.frame(tax_table(OSU))
nlevels(OSUtax$Family)
# 204 Families
OSUpretax <- as.data.frame(tax_table(OSUpre))
nlevels(OSUpretax$Family)
# 177 Families
OSUposttax <- as.data.frame(tax_table(OSUpost))
nlevels(OSUposttax$Family)
# 133 Families

#---------OSU: Determine relative abundance of Endozoicomonadaceae in parents---------

# Calculate relative abundance of each ASV, convert phyloseq information to data frame, and filter out all non-Endozoicomonadaceae ASVs
OSUrel <- transform_sample_counts(OSU, function(OTU) OTU/sum(OTU))
OSU.melted <- psmelt(OSUrel)
OSU.endo <- OSU.melted %>% filter(Family == "Endozoicomonadaceae") %>% filter(Abundance > 0)

# Reorganize data frame to contain relative abundance of Endozoicomonadaceae in each sample
OSU.endo <- OSU.endo %>% ungroup() %>% group_by(SampleID) %>% 
  summarize(SampleID = SampleID, EndoAbun = sum(Abundance), Colony = Colony, Treatment = Treatment, Stage = Stage, CloneID = CloneID) %>% distinct()

# Calculate mean and standard error of relative abundance of Endozoicomonadaceae present in all samples
mean(OSU.endo$EndoAbun)
# 0.9915409
sd(OSU.endo$EndoAbun)/sqrt(40)
# 0.002716017

#---------OSU: Statistical analyses of parent bacterial communities by treatment--------- 

# Function for reading in distance matrices
read_dist <- function(file) {
  pre.dist <- read.csv(file)
  dist <- pre.dist[2:ncol(pre.dist)]
  colnames(dist) <- t(pre.dist[1])
  rownames(dist) <- t(pre.dist[1])
  invisible(dist)
}

# Bacterial community composition of pre-planulation parents in control vs nitrate treatment
pre.w.dist <- as.dist(read_dist("Unifrac_pre_dist_PVP16S.csv"))
pre_df <- data.frame(sample_data(OSUpre))
pre_df <- tibble::rownames_to_column(pre_df, "Sample")

adonis(pre.w.dist ~ Treatment, data = pre_df)
#           Df  SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)  
# Treatment  1 0.00020808 0.00020808  1.9929 0.09968  0.086 .
# Residuals 18 0.00187940 0.00010441         0.90032         
# Total     19 0.00208748                    1.00000         

# Bacterial community composition of post-planulation parents in control vs nitrate treatment
post.w.dist <- as.dist(read_dist("Unifrac_post_dist_PVP16S.csv"))
post_df <- data.frame(sample_data(OSUpost))
post_df <- tibble::rownames_to_column(post_df, "Sample")

adonis(post.w.dist ~ Treatment, data = post_df)
#           Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)
# Treatment  1 0.0000890 8.9000e-05 0.94121 0.04969  0.671
# Residuals 18 0.0017021 9.4559e-05         0.95031       
# Total     19 0.0017910                    1.00000

# Alpha diversity of bacterial communities of pre-planulation parents in control vs nitrate treatment
prerich <- estimate_richness(OSUpre, measures = c("Shannon", "Simpson"))
samplepre <- row.names(prerich)
prerichdf <- data.frame(Sample = samplepre, prerich)
prerichdiv <- merge(prerichdf, pre_df, by = "Sample")

# Test for normal distribution of alpha diversity values
shapiro.test(prerichdiv$Shannon) # p = 0.0001265, not normal
shapiro.test(prerichdiv$Simpson) # p = 0.006236, not normal

# Neither set of diversity values was normal, so assess differences in alpha diversity using a Wilcoxon rank sum test
wilcox.test(Simpson ~ Treatment, data = prerichdiv) # p = 0.315
wilcox.test(Shannon ~ Treatment, data = prerichdiv) # p = 0.5787

# Alpha diversity of bacterial communities of post-planulation parents in control vs nitrate treatment
postrich <- estimate_richness(OSUpost, measures = c("Shannon", "Simpson"))
samplepost <- row.names(postrich)
postrichdf <- data.frame(Sample = samplepost, postrich)
postrichdiv <- merge(postrichdf, post_df, by = "Sample")

# Test for normal distribution of alpha diversity values
shapiro.test(postrichdiv$Shannon) # p = 8.509e-08, not normal
shapiro.test(postrichdiv$Simpson) # p = 4.043e-06, not normal

# Neither set of diversity values was normal, so assess differences in alpha diversity using a Wilcoxon rank sum test
wilcox.test(Simpson ~ Treatment, data = postrichdiv) # p = 0.1431
wilcox.test(Shannon ~ Treatment, data = postrichdiv) # p = 0.1431

# Bacterial community dispersion of pre-planulation parents in control vs nitrate treatment
permutest(betadisper(pre.w.dist, pre_df$Treatment))
#           Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0000501 5.0101e-05 0.5573    999  0.461
# Residuals 18 0.0016182 8.9897e-05 

# Bacterial community dispersion of post-planulation parents in control vs nitrate treatment
permutest(betadisper(post.w.dist, post_df$Treatment))
#           Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00006656 6.6564e-05 0.7195    999  0.637
# Residuals 18 0.00166520 9.2511e-05 

#---------OSU: Statistical analyses of parent bacterial communities by clone--------- 

# Bacterial community composition of pre-planulation parents by clone
adonis(pre.w.dist ~ CloneID, data = pre_df)
#           Df  SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)  
# CloneID    5 0.00141204 2.8241e-04  5.8536 0.67644  0.023 *
# Residuals 14 0.00067544 4.8245e-05         0.32356         
# Total     19 0.00208748                    1.00000

# Bacterial community composition of post-planulation parents by clone
adonis(post.w.dist ~ CloneID, data = post_df)
#           Df  SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)
# CloneID    5 0.00091179 1.8236e-04  2.9036 0.50908  0.106
# Residuals 14 0.00087926 6.2804e-05         0.49092       
# Total     19 0.00179105                    1.00000

#---------OSU: Statistical analyses comparing parents pre- and post-planulation--------- 

# Will only use samples from control treatment for these analyses, so need to subset phyloseq object
Control <- subset_samples(OSU, Treatment == "Control")
Control <- prune_taxa(taxa_sums(Control) > 0, Control)

# Bacterial community composition between pre- and post-planulation parents
control.w.dist <- as.dist(read_dist("Unifrac_control_dist_PVP16S.csv"))
control_df <- data.frame(sample_data(Control))
control_df <- tibble::rownames_to_column(control_df, "Sample")

adonis(control.w.dist ~ Stage, data = control_df)
#           Df  SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)  
# Stage      1 0.00004217 4.2173e-05  1.2059 0.06279  0.028 *
# Residuals 18 0.00062950 3.4972e-05         0.93721         
# Total     19 0.00067167                    1.00000 

# Alpha diversity of bacterial communities between pre- and post-planulation parents
controlrich <- estimate_richness(Control, measures = c("Shannon", "Simpson"))
samplecontrol <- row.names(controlrich)
controlrichdf <- data.frame(Sample = samplecontrol, controlrich)
controlrichdiv <- merge(controlrichdf, control_df, by = "Sample")

# Test for normal distribution of alpha diversity values
shapiro.test(controlrichdiv$Shannon) # p = 3.875e-06, not normal
shapiro.test(controlrichdiv$Simpson) # p = 3.405e-05, not normal

# Neither set of diversity values was normal, so assess differences in alpha diversity using a Wilcoxon rank sum test
wilcox.test(Simpson ~ Stage, data = controlrichdiv) # p = 0.1051
wilcox.test(Shannon ~ Stage, data = controlrichdiv) # p = 0.01854

# Bacterial community dispersion between pre- and post-planulation parents
permutest(betadisper(control.w.dist, control_df$Stage))
#           Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00004276 4.2758e-05 1.4589    999  0.288
# Residuals 18 0.00052756 2.9309e-05

#---------GG: Explore phyloseq object and filter reads from planulae only---------

# Look at the distribution of sequence depth
hist(sample_sums(GG), breaks=length(sample_sums(GG)))

# Identify samples will low sequence depth (<10,000 reads)
sums <- sample_sums(GG)
sums[sums<=10000]
# GG10N6 GG6C10  GG6C6  GG6C7  GG8C6  GGNeg 
#   6300   4936   9397   8562   9668     96

# This phyloseq object contains all individual planulae samples and a subset of parent samples. To assess the number of reads and ASVs in planluae alone, will first subset phyloseq object to only include planulae and negative control
planula <- subset_samples(GG, Stage == "Planula" | SampleID == "Neg")
planula <- prune_taxa(taxa_sums(planula) > 0, planula)

# Filter out samples with less than 10000 reads, unassigned reads, and reads that identify as Mitochondria or Chloroplast DNA
planula2 <- subset_samples(planula, SampleID != "10N6" & SampleID != "6C10" & SampleID != "6C6" & SampleID != "6C7" & SampleID != "8C6")
planula2 <- prune_taxa(taxa_sums(planula2) > 0, planula2)
planula3 <- subset_taxa(planula2, (Kingdom != "NA"))
planula4 <- subset_taxa(planula3, (Family != "Mitochondria"))
planula5 <- subset_taxa(planula4, (Family != "Chloroplast"))

# Check the number of reads and ASVs remaining at each filtering step
sum(taxa_sums(planula)); ncol(otu_table(planula))
# 2726708 reads; 2935 ASVs
sum(taxa_sums(planula2)); ncol(otu_table(planula2))
# 2687845 reads; 2901 ASVs
sum(taxa_sums(planula3)); ncol(otu_table(planula3))
# 2687845 reads; 2901 ASVs
sum(taxa_sums(planula4)); ncol(otu_table(planula4))
# 2325996 reads; 2410 ASVs
sum(taxa_sums(planula5)); ncol(otu_table(planula5))
# 2325996 reads; 2410 ASVs

# Rename back to "planula" for simplicity and remove transitional objects
planula <- planula5
rm(planula2, planula3, planula4, planula5)

#---------GG: Remove ASVs found in negative control from planulae only---------
# Negative control was an extraction from molecular water sent for sequencing, called "GGNeg"

# Make a phyloseq object with just the negative sample and negative taxa
GGNeg <- subset_samples(planula, sample_names(planula) == "GGNeg")
GGNeg <- prune_taxa(taxa_sums(GGNeg) > 0, GGNeg)

# Extract negative ASV list
GGNeg.ASVs <- rownames(data.frame(tax_table(GGNeg)))
# There are 5 unique sequences in GGNeg

# Extract full list of ASVs from all planula samples
planula.ASVs <- rownames(data.frame(tax_table(planula)))

# Make a list of the ASVs to keep by subtracting out the ones from the negative
planulakeeptaxa <- planula.ASVs[!(planula.ASVs %in% GGNeg.ASVs)]

# Remove the taxa found in the negative control from the planula data
planula <- prune_taxa(planulakeeptaxa, planula)
planula <- prune_taxa(taxa_sums(planula) > 0, planula)

# Confirm there are no reads remaining for GGNeg and remove from phyloseq object
sample_sums(planula)
planula <- subset_samples(planula, SampleID != "Neg")

nrow(tax_table(planula))
# 2405 ASVs

#---------GG: Filter reads from planulae + subset of parents---------
# Because I will do one analysis comparing the planulae to the parents, I will also filter and subset the whole GG phyloseq object

# Filter out samples with less than 10000 reads, unassigned reads, and reads that identify as Mitochondria or Chloroplast DNA
GG2 <- subset_samples(GG, SampleID != "10N6" & SampleID != "6C10" & SampleID != "6C6" & SampleID != "6C7" & SampleID != "8C6")
GG2 <- prune_taxa(taxa_sums(GG2) > 0, GG2)
GG3 <- subset_taxa(GG2, (Kingdom != "NA"))
GG4 <- subset_taxa(GG3, (Family != "Mitochondria"))
GG5 <- subset_taxa(GG4, (Family != "Chloroplast"))

# Rename back to "GG" for simplicity and remove transitional objects
GG <- GG5
rm(GG2, GG3, GG4, GG5)

#---------GG: Remove ASVs found in negative control from planulae + subset of parents---------
# Negative control was an extraction from molecular water sent for sequencing, called "GGNeg"

# Extract full list of ASVs from all GG samples
GG.ASVs <- rownames(data.frame(tax_table(GG)))

# Make a list of the ASVs to keep by subtracting out the ones from the negative
GGkeeptaxa <- GG.ASVs[!(GG.ASVs %in% GGNeg.ASVs)]

# Remove the taxa found in the negative control from the data
GG <- prune_taxa(GGkeeptaxa, GG)
GG <- prune_taxa(taxa_sums(GG) > 0, GG)

# Confirm there are no reads remaining for GGNeg and remove from phyloseq object
sample_sums(GG)
GG <- subset_samples(GG, SampleID != "Neg")

nrow(tax_table(GG))
# 3005 ASVs

#---------GG: Statistical analyses comparing parents to planulae ---------

# Bacterial community composition of all parents vs planulae
GG.w.dist <- as.dist(read_dist("Unifrac_GG_dist_PVP16S.csv"))
GG_df <- data.frame(sample_data(GG))
GG_df <- tibble::rownames_to_column(GG_df, "Sample")

# Group pre- and post-planulation parents into one group "Parent"
GG_df$Stage <- as.character(GG_df$Stage)
GG_df$Stage[GG_df$Stage == "ParentPre" | GG_df$Stage == "ParentPost"] <- "Parent"
GG_df$Stage <- as.factor(GG_df$Stage)

adonis(GG.w.dist ~ Stage, data = GG_df)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Stage      1    1.8333 1.83335  53.174 0.37401  0.001 ***
# Residuals 89    3.0686 0.03448         0.62599           
# Total     90    4.9019                 1.00000 

# Bacterial community dispersion of all parents vs planulae
permutest(betadisper(GG.w.dist, GG_df$Stage))
#           Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.43796 0.43796 123.05    999  0.001 ***
# Residuals 89 0.31678 0.00356   

#---------GG: Statistical analyses of planulae bacterial communities by treatment--------- 

# Bacterial community composition of planulae in control vs nitrate treatment
planula.w.dist <- as.dist(read_dist("Unifrac_planula_dist_PVP16S.csv"))
planula_df <- data.frame(sample_data(planula))
planula_df <- tibble::rownames_to_column(planula_df, "Sample")

adonis(planula.w.dist ~ Treatment, data = planula_df)
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Treatment  1   0.04826 0.048262  1.1676 0.01574  0.273
# Residuals 73   3.01732 0.041333         0.98426       
# Total     74   3.06558                  1.00000

# Alpha diversity of bacterial communities of planulae in control vs nitrate treatment
planrich <- estimate_richness(planula, measures = c("Shannon", "Simpson"))
sampleplan <- row.names(planrich)
planrichdf <- data.frame(Sample = sampleplan, planrich)
planrichdiv <- merge(planrichdf, planula_df, by = "Sample")

# Test for normal distribution of alpha diversity values
shapiro.test(planrichdiv$Shannon) # p = 0.5801, normal
shapiro.test(planrichdiv$Simpson) # p = 6.508e-09, not normal

# Shannon diversity values are normally distributed, so assess differences in Shannon diversity using a t-test. Simpson diversity values are not normally distributed, so assess differences in Simpson diversity using a Wilcoxon rank sum test
wilcox.test(Simpson ~ Treatment, data = planrichdiv) # p = 0.1233
t.test(Shannon ~ Treatment, data = planrichdiv) # p = 0.07106

# Bacterial community dispersion of planulae in control vs nitrate treatment
permutest(betadisper(planula.w.dist, planula_df$Treatment))
#           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00303 0.0030298 0.7324    999  0.388
# Residuals 73 0.30198 0.0041367

#---------GG: Statistical analyses of planulae bacterial communities by clone--------- 

# Bacterial community composition of planulae by parent clone
adonis(planula.w.dist ~ CloneID, data = planula_df)
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# CloneID    1   0.06939 0.069389  1.6906 0.02263  0.078 .
# Residuals 73   2.99619 0.041044         0.97737         
# Total     74   3.06558                  1.00000

#---------Bubble plot of parent bacterial communities---------

# Transform sample counts to relative abundance, filter out sequences with <0.1% abundance in a sample, and prepare data to be used in ggplot2
OSU <- tax_glom(OSU, rank_names(OSU)[5])
OSU.rel <- transform_sample_counts(OSU, function(OTU) OTU/sum(OTU))
OSU.rel <- prune_taxa(taxa_sums(OSU.rel) > 0.001, OSU.rel)
OSU.melted <- psmelt(OSU.rel)

# Add a new column to data frame with percent abundance (ex. 85% instead of 0.85)
perc <- OSU.melted$Abundance * 100
OSU.melted <- cbind(OSU.melted, perc)

# Order Family levels by the Phylum they belong to for plot aesthetics
OSU.melted$Family <- factor(
  OSU.melted$Family, 
  levels = unique(OSU.melted$Family[order(OSU.melted$Phylum)]), 
  ordered=TRUE)

# Add new column to remove "Pr" or "Po" from SampleID and order by CloneID for plot aesthetics
OSU.melted$newID <- substr(OSU.melted$SampleID, 1, nchar(as.character(OSU.melted$SampleID)) - 2)
OSU.melted$newID <- factor(
  OSU.melted$newID, 
  levels = unique(
    c("10C", "10N", "3C", "3N", "4C", "5C", "5N", "6C", "6N", "9C", "9N", "2C", "2N", "1C", "4N", "8C", "8N", "1N", "7C", "7N", "10C", "10N", "3C", "3N", "4C", "5C", "5N", "6C", "6N", "9C", "9N",  "2C", "2N", "1C", "4N", "8C", "8N",  "1N", "7C", "7N"), 
    ordered=TRUE))

# Rename levels of Stage for plot aesthetics
OSU.melted$Stage <- factor(OSU.melted$Stage, levels = c("ParentPre", "ParentPost"))
stage.labs <- c("Pre-Planulation", "Post-Planulation")
names(stage.labs) <- c("ParentPre", "ParentPost")

# Create the bubble plot
ParentBubble <- ggplot(OSU.melted, aes(x = newID, y = Family)) + 
  geom_vline(xintercept = c("1N", "2N", "3N", "4N", "5N", "6N", "7N", "8N", "9N", "10N"), 
             color = "grey60") +
  geom_point(aes(size = perc, fill = Phylum), shape = 21, alpha = 0.75, stroke = 0.1) +
  scale_size_continuous(limits = c(0.000001, 100), range = c(2, 15), breaks = c(1, 10, 50, 75)) + 
  facet_grid( ~ Stage, scales="free_x", labeller = labeller(Stage = stage.labs)) +
  labs(y = "Family", size = "Relative Abundance (%)", fill = "Phylum") +
  guides(fill = guide_legend(override.aes = list(size = 8))) +
  theme(legend.key = element_blank(),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 16, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 14), 
        legend.text = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 16), 
        strip.text.x = element_text(size = 18),
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +  
  scale_fill_manual(values = c("grey46", "mediumpurple1", "wheat3", "darkorange", "hotpink2", "darkolivegreen3", "firebrick", "darkslategray1", "goldenrod1", "royalblue3", "lightcoral", "darkcyan", "darkmagenta")) + 
  scale_y_discrete(limits = rev(levels(OSU.melted$Family))) 
ParentBubble

#---------Relative abundance stacked bar plot of planulae bacterial communities---------

# Transform sample counts to relative abundance and remove any OTUs with 0% abundance
plan.rel <- transform_sample_counts(planula, function(OTU) OTU/sum(OTU))
plan.rel <- prune_taxa(taxa_sums(plan.rel) > 0, plan.rel)

# Create a new column "PhyFam" that contains Phylum_Class_Family and prep data for ggplot2
tax.df <- as.data.frame(tax_table(plan.rel)) %>% 
  unite("PhyFam", c(Phylum, Class, Family), sep = "_", remove = FALSE)
tax_table(plan.rel) <- as.matrix(tax.df)
plan.melted <- psmelt(plan.rel)

# Bin all ASVs with <5% abundance into a category called "Other"
plan.melted <- plan.melted %>% filter(Abundance > 0)
plan.melted$PhyFam <- as.character(plan.melted$PhyFam)
plan.melted$PhyFam[plan.melted$Abundance < 0.05] <- "Other"

# Rename PhyFam to "Phylum_Class_Family"
colnames(plan.melted)[11] <- "Phylum_Class_Family"

# Order samples properly
six <- c("6C1", "6C2", "6C3", "6C4", "6C5", "6C6", "6C7", "6C8", "6C9", "6C10", "6N1", "6N2", "6N3", "6N4", "6N5", "6N6", "6N7", "6N8", "6N9", "6N10")
eight <- c("8C1", "8C2", "8C3", "8C4", "8C5", "8C6", "8C7", "8C8", "8C9", "8C10", "8N1", "8N2", "8N3", "8N4", "8N5", "8N6", "8N7", "8N8", "8N9", "8N10")
nine <- c("9C1", "9C2", "9C3", "9C4", "9C5", "9C6", "9C7", "9C8", "9C9", "9C10", "9N1", "9N2", "9N3", "9N4", "9N5", "9N6", "9N7", "9N8", "9N9", "9N10")
ten <- c("10C1", "10C2", "10C3", "10C4", "10C5", "10C6", "10C7", "10C8", "10C9", "10C10", "10N1", "10N2", "10N3", "10N4", "10N5", "10N6", "10N7", "10N8", "10N9", "10N10")
plan.melted$SampleID <- factor(plan.melted$SampleID, levels = c(six, eight, nine, ten))

#Assign color gradients
Actino <- colorRampPalette(c("lavender", "darkorchid4"))
Bactero <- colorRampPalette(c("peachpuff2", "coral4"))
Chloro <- colorRampPalette(c("cyan1", "cyan3"))
Cyano <- c("slategray1", "slategray3")
Deino <- c("khaki1", "gold")
Eury <- c("lightslateblue")
Firm <- colorRampPalette(c("pink", "palevioletred3"))
Fuso <- c("navyblue")
Gemma <- c("darkorchid1")
Lenti <- c("black")
Pates <- c("coral")
Plant <- colorRampPalette(c("lemonchiffon1", "lightgoldenrod3"))
#Proteo <- colorRampPalette(c("darkseagreen1", "darkgreen"))
#Proteo <- Proteo(26)
Proteo <- c("#C1FFC1", "#B9F8B9", "#B1F2B1", "#A9ECA9", "#A2E6A2", "#9AE09A", "#92D992", "#8AD38A", "#83CD83", "#7BC77B", "#73C173", "chartreuse", "#64B464", "#5CAE5C", "#54A854", "#4DA24D", "#459B45", "#3D953D", "#368F36", "#2E892E", "#268326", "#1E7C1E", "#177617", "#0F700F", "#076A07", "#006400")
Verru <- colorRampPalette(c("steelblue1", "steelblue4"))
Other <- c("grey81")
color <- c(Other, Actino(12), Bactero(8), Chloro(3), Cyano, Deino, Eury, Firm(13), Fuso, Gemma, Lenti, Pates, Plant(5), Proteo, Verru(3))

#Order bars so "Other" is on the bottom of each stack
PhyFam <- sort(unique(plan.melted$Phylum_Class_Family), decreasing = FALSE)
PhyFam <- PhyFam[-45]
PhyFam <- c("Other", PhyFam)
plan.melted$Phylum_Class_Family <- factor(plan.melted$Phylum_Class_Family, levels = PhyFam)
plan.melted$Phylum_Class_Family <- as.factor(plan.melted$Phylum_Class_Family)

#Create the bar plot
PlanBar <- ggplot(data = plan.melted, aes(x = SampleID, y = Abundance, fill = Phylum_Class_Family)) +
  geom_bar(aes(), stat = "identity", width = 0.8) +
  facet_grid( ~ CloneID, scales = "free_x", space = "free") +
  scale_fill_manual(values = color) +
  ylab("Relative Abundance of Bacterial Taxa (>5% abundance)") +
  theme(panel.border = element_rect(fill = NA, colour = "black",size = 0.7),
        panel.spacing = unit(0.5, "lines"), 
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(hjust = 1, size = 10, angle = 45),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank()) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6.75)) +
  theme(strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"))
PlanBar
