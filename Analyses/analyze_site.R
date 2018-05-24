##### analyze_site.R #####
# Kuan-lin Huang @ WashU 2017 July
# updated in 2017 Dec
# updated 2018 April
# determine if CPTAC sites are more enriched in clusters (hybrid, pho only) than dbPAF sites and simulated sites

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")

# MC3 driver score
mc3_score_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/somatic/Driver_BaileyCell2018/Mutation.CTAT.3D.Scores.txt.gz"
mc3_score = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(mc3_score_f))
colnames(mc3_score) = gsub("\\.","",colnames(mc3_score))

# CPTAC sites in clusters
clustered_cptac_site= cptac_site[paste(cptac_site$hgnc_symbol,cptac_site$amino_acid_residue) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene),]
cat("Total number of CPTAC sites clustered:",nrow(clustered_cptac_site),"\n")
rank_vectors(clustered_cptac_site$hgnc_symbol)
tn = "output/clusteredCPTACsites.txt"
write.table(clustered_cptac_site, quote=F, sep="\t", file = tn, row.names = F)


##### pairwise #####
# reference panel of phosphosites
PTMcosmo_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/input/PTMPortal_ptm_sites_dump.phospho_only.with_enst.tsv"
PTMcosmo = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = PTMcosmo_f)
PTMcosmo$Mutation_Gene = paste("p.",PTMcosmo$residue,PTMcosmo$p_coord,sep="")
PTMcosmo_map = PTMcosmo[,c("Ensembl.Transcript","p_coord","Mutation_Gene")]
cptac_site$Position = gsub("p.[A-Z]","",cptac_site$amino_acid_residue)
cptac_site_map = cptac_site[,c("ensembl_transcript_id","Position","amino_acid_residue")]
colnames(PTMcosmo_map) = c("Transcript2","Position","originalLabel")
colnames(cptac_site_map) = c("Transcript2","Position","originalLabel")
phosphosites_map = rbind(cptac_site_map,PTMcosmo_map)
phosphosites_map = phosphosites_map[!duplicated(paste(phosphosites_map$Transcript,phosphosites_map$Position)),]

# mapped mutation and phosphosites
site_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201805/3D_Proximity.musites.gz"
site = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(site_f))
site_uniq = site[!duplicated(paste(site$Gene1,site$Mutation1,site$Transcript2,site$TranscriptPosition2)),]

# examine the average distances of clustered sites
clustered_sites = site_uniq[paste(site_uniq$Gene1,site_uniq$Mutation1) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)
                       & paste(site_uniq$Transcript2,site_uniq$Site2) %in% paste(annotated_cluster$Transcript,annotated_cluster$Mutation_Gene),]
site_uniq$StructuralDistance = as.numeric(gsub(" .*","",site_uniq$DistanceInfo))
site_uniq$LinearDistance = as.numeric(site_uniq$LinearDistance)
p = ggplot(site_uniq,aes(x = LinearDistance, y=StructuralDistance))
p = p + geom_point(alpha=0.01,stroke=0, shape=16) + scale_x_log10(breaks = c(1,5,10,100,1000))#+ xlim(0,1000) 
p = p + ylim(0,20) + labs(title = "Distance between co-clustered phosphosites and mutations", 
                          x = "Linear Distance (# of amino acid residues)", y = "Distance on 3D Protein Structures (Angstrom, Å)")
p = p + theme_bw() + theme(legend.position = "bottom")
p
fn = paste("output/clustered_site_distances_linear_vs_3D.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

mutations = site[,c("Gene1","Mutation1")]
mutations = mutations[!duplicated(paste(mutations$Gene1,mutations$Mutation1)),]
mutations$Position = as.numeric(gsub("p.[A-Z](.*)[A-Z]","\\1",mutations$Mutation1))

phosphosites = site[,c("Gene2","Transcript2","Site2")]
phosphosites = phosphosites[!duplicated(paste(phosphosites$Transcript2,phosphosites$Site2)),]
phosphosites$Site2 = gsub("p. ","p.",phosphosites$Site2) # some residues had a gap...
phosphosites$Position = as.numeric(gsub("p.[A-Z]*([0-9]*)[A-Z]*","\\1",phosphosites$Site2))
#clean phosphosite residues here:
phosphosites_merge = merge(phosphosites,phosphosites_map,by=c("Transcript2","Position")) # some are still not mapped but ignoring those for now, none of those clustered...
cat("Number of residues with inconsistent residues from the original phosphosite files (True being inconsistent):\n")
table(phosphosites_merge$Site2 != phosphosites_merge$originalLabel)

phosphosites_merge$Residue = gsub("p.([A-Z])[0-9]+","\\1",phosphosites_merge$originalLabel)
background_rate = data.frame(table(phosphosites_merge$Residue))
background_rate = background_rate[background_rate$Var1 != "H",] #negletible
background_rate$Var1 = factor(background_rate$Var1,levels=c("S","T","Y"))
background_rate$Var2 = "Mapped"

p = ggplot(background_rate,aes(x = Var2, y=Freq, fill=Var1))
p = p + geom_bar(stat = "identity") + labs(x = "", y = "Count of phosphosites")
p = p + theme_bw() + theme(legend.position = "bottom")
p
fn = paste("output/Data_201805_mapped_sites_by_residue.pdf",sep="_")
ggsave(fn, w=1,h=4,useDingbat=F)

# compare rate of STY to background rate of STY "other" residues
annotated_clusterPTM = annotated_cluster[annotated_cluster$Alternate=="ptm",]
annotated_clusterPTM$ref = gsub("p.([A-Z])[0-9]+","\\1",annotated_clusterPTM$Mutation_Gene)
annotated_clusterPTM_h = annotated_clusterPTM[annotated_clusterPTM$Type=="Hybrid",]
clustered_rate = data.frame(table(annotated_clusterPTM_h$ref))
clustered_rate$Var2 = "InHybridClusters"

p = ggplot(clustered_rate,aes(x = Var2, y=Freq, fill=Var1))
p = p + geom_bar(stat = "identity") + labs(x = "", y = "Count of phosphosites")
p = p + theme_bw() + theme(legend.position = "bottom")
p
fn = paste("output/Data_201805_coclustered_sites_by_residue.pdf",sep="_")
ggsave(fn, w=1,h=4,useDingbat=F)

cat("Enrichment for T:\n")
fisher.test(alternative = "greater",matrix(ncol=2,nrow=2,c(sum(annotated_clusterPTM_h$ref=="T",na.rm = T),sum(phosphosites_merge$Residue=="T",na.rm = T)-sum(annotated_clusterPTM_h$ref=="T",na.rm = T),
            sum(annotated_clusterPTM_h$ref !="T" ,na.rm = T),sum(phosphosites_merge$Residue %in% c("S","Y"),na.rm = T)-sum(annotated_clusterPTM_h$ref !="T" ,na.rm = T))))
cat("Enrichment for S:\n")
fisher.test(alternative = "greater",matrix(ncol=2,nrow=2,c(sum(annotated_clusterPTM_h$ref=="S",na.rm = T),sum(phosphosites_merge$Residue=="S",na.rm = T)-sum(annotated_clusterPTM_h$ref=="S",na.rm = T),
                                   sum(annotated_clusterPTM_h$ref !="S" ,na.rm = T),sum(phosphosites_merge$Residue %in% c("S","Y"),na.rm = T)-sum(annotated_clusterPTM_h$ref !="S" ,na.rm = T))))
cat("Enrichment for Y:\n")
fisher.test(alternative = "greater",matrix(ncol=2,nrow=2,c(sum(annotated_clusterPTM_h$ref=="Y",na.rm = T),sum(phosphosites_merge$Residue=="Y",na.rm = T)-sum(annotated_clusterPTM_h$ref=="Y",na.rm = T),
                                   sum(annotated_clusterPTM_h$ref !="Y" ,na.rm = T),sum(phosphosites_merge$Residue %in% c("S","Y"),na.rm = T)-sum(annotated_clusterPTM_h$ref !="Y" ,na.rm = T))))

# clustered phosphosites compared to Direct/proximal
annotated_clusterMUT = annotated_cluster[annotated_cluster$Alternate != "ptm",]
annotated_clusterPTM_h$SiteType = "Clustered"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position-1)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position-2)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position+1)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position+2)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position)] = "Direct"
cat("Clustered phosphosites also showing up as direct or proximal:\n")
table(annotated_clusterPTM_h$SiteType)

# look at phosphosite ratio within each gene
background_by_gene = data.frame(table(phosphosites_merge$Gene2))
clustered_by_gene = data.frame(table(annotated_clusterPTM_h$Gene_Drug))
#all_genes = rbind(background_by_gene,clustered_by_gene)
colnames(background_by_gene) = c("Gene","Mapped")
colnames(clustered_by_gene) = c("Gene","Clustered")
all_by_gene = merge(background_by_gene, clustered_by_gene, by="Gene")
all_by_gene$NotClustered = all_by_gene$Mapped - all_by_gene$Clustered
all_by_gene_m = melt(all_by_gene,by="Gene")
all_by_gene$Rate = all_by_gene$Clustered/all_by_gene$Mapped
top_rate = all_by_gene$Gene[all_by_gene$Rate>0.5 & all_by_gene$Mapped > 5]
all_by_gene_m = all_by_gene_m[all_by_gene_m$variable!="Mapped",]
p = ggplot(all_by_gene_m[all_by_gene_m$Gene %in% top_rate,],aes(x = Gene, y= value, fill=variable))
p = p + geom_bar(stat = "identity") + labs(x = "", y = "Count of phosphosites") + coord_flip()
p = p + theme_bw() + theme(legend.position = "bottom")
p
fn = paste("output/Data_201805_coclustered_sites_by_gene_ratio.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

tn = "output/coclustered_phosphosite_by_gene_ratio.txt"
write.table(all_by_gene, quote=F, sep="\t", file = tn, row.names = F)

##### examine whether clustered mutations have higher mc3 score ######
mc3_score$Position = as.numeric(gsub("p.[A-Z](.*)[A-Z]","\\1",mc3_score$protein_change))
mc3_score = mc3_score[mc3_score$gene %in% annotated_cluster$Gene_Drug,]
mc3_score$Type = "None"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$protein_change) %in% paste(annotated_cluster$Transcript,annotated_cluster$Mutation_Gene)] = "Clustered"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(phosphosites$Transcript2,phosphosites$Position-1)] = "Proximal"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(phosphosites$Transcript2,phosphosites$Position-2)] = "Proximal"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(phosphosites$Transcript2,phosphosites$Position+1)] = "Proximal"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(phosphosites$Transcript2,phosphosites$Position+2)] = "Proximal"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(phosphosites$Transcript2,phosphosites$Position)] = "Direct"
table(mc3_score$Type)
# mc3_score$Type[paste(mc3_score$gene,mc3_score$protein_change) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)] = "Clustered"
# mc3_score$Type[paste(mc3_score$gene,mc3_score$Position) %in% paste(phosphosites$Gene2,phosphosites$Position-1)] = "Proximal"
# mc3_score$Type[paste(mc3_score$gene,mc3_score$Position) %in% paste(phosphosites$Gene2,phosphosites$Position-2)] = "Proximal"
# mc3_score$Type[paste(mc3_score$gene,mc3_score$Position) %in% paste(phosphosites$Gene2,phosphosites$Position+1)] = "Proximal"
# mc3_score$Type[paste(mc3_score$gene,mc3_score$Position) %in% paste(phosphosites$Gene2,phosphosites$Position+2)] = "Proximal"
# mc3_score$Type[paste(mc3_score$gene,mc3_score$Position) %in% paste(phosphosites$Gene2,phosphosites$Position)] = "Direct"
# table(mc3_score$Type)
# # see if proximal and direct being in the cluster have higher mutation score: not true
# mc3_score$Type[mc3_score$Type== "Proximal" & paste(mc3_score$gene,mc3_score$protein_change) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)] = "Clustered & Proximal"
# mc3_score$Type[mc3_score$Type== "Direct" & paste(mc3_score$gene,mc3_score$protein_change) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)] = "Clustered & Direct"

cat("Testing whether overlapping or clustered mutations have higher functional score than others","\n")
wilcox.test(mc3_score$eigenscorefunctional[mc3_score$Type=="Clustered"],mc3_score$eigenscorefunctional[mc3_score$Type=="None"])
wilcox.test(mc3_score$eigenscorefunctional[mc3_score$Type=="Clustered"],mc3_score$eigenscorefunctional[mc3_score$Type=="Proximal"])
wilcox.test(mc3_score$eigenscorefunctional[mc3_score$Type=="Clustered"],mc3_score$eigenscorefunctional[mc3_score$Type=="Direct"])

cat("Testing whether overlapping or clustered mutations have higher cancer-specific functional score than others","\n")
wilcox.test(mc3_score$eigenscorecancer[mc3_score$Type=="Clustered"],mc3_score$eigenscorecancer[mc3_score$Type=="None"])
wilcox.test(mc3_score$eigenscorecancer[mc3_score$Type=="Clustered"],mc3_score$eigenscorecancer[mc3_score$Type=="Proximal"])
wilcox.test(mc3_score$eigenscorecancer[mc3_score$Type=="Clustered"],mc3_score$eigenscorecancer[mc3_score$Type=="Direct"])

# fit = glm(data= mc3_score, eigenscorefunctional~Type+gene)
# summary(fit)
mc3_score$Type = factor(mc3_score$Type,levels = c("None","Clustered","Proximal","Direct"))

for (score in colnames(mc3_score)[5:17]){
  p = ggplot(mc3_score,aes_string(x = "Type", y = score, fill="Type"))
  p = p + geom_violin(alpha=0.8)
  p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8)
  p = p + getTypeFillScale()
  p = p + theme_bw()
  p = p + labs(y = "Score", title = score, x = "Interaction with Phosphosites")
  p
  fn = paste("output/mutation",score,"vs_interaction_wPhospho.pdf",sep="_")
  ggsave(fn,h=5, w = 5,useDingbat=F)
}

p = ggplot(mc3_score,aes(x = Type, y = eigenscorefunctional, fill=Type))
p = p + geom_violin(alpha=0.8)
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8)
p = p + getTypeFillScale()
p = p + theme_bw()
p = p + labs(y = "Mutation Eigen Functional Score", x = "Interaction with Phosphosites")
p
fn = paste("output/mutation_functionalscore_vs_interaction_wPhospho.pdf",sep="_")
ggsave(fn,h=5, w = 5,useDingbat=F)

mc3_score$Type = factor(mc3_score$Type,levels = c("None","Clustered","Proximal","Direct"))
p = ggplot(mc3_score,aes(x = Type, y = eigenscorecancer, fill=Type))
p = p + geom_violin(alpha=0.8)
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.8)
p = p + getTypeFillScale() + ylim(-5,10)
p = p + theme_bw()
p = p + labs(y = "Mutation Eigen Cancer Functional Score", x = "Interaction with Phosphosites")
p
fn = paste("output/mutation_cancerscore_vs_interaction_wPhospho.pdf",sep="_")
ggsave(fn,h=5, w = 5,useDingbat=F)

# look at gene distribution of mutations
annotated_clusterMUT_h = annotated_clusterMUT[annotated_clusterMUT$Type=="Hybrid",]
annotated_clusterMUT_h$SiteType = "Clustered"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_cluster_h$Transcript,annotated_clusterPTM_h$Position-1)] = "Proximal"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position-2)] = "Proximal"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position+1)] = "Proximal"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position+2)] = "Proximal"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position)] = "Direct"
cat("Clustered mutations showing up as direct or proximal:\n")
table(annotated_clusterMUT_h$SiteType)

top_genes_cluster = names(rank_vectors(annotated_clusterMUT_h$Gene_Drug[annotated_clusterMUT_h$SiteType == "Clustered"],n=15))
# top_genes_cluster = names(rank_vectors(mc3_score$gene[mc3_score$Type == "Clustered"],n=10))
# top_genes_proxi = names(rank_vectors(mc3_score$gene[mc3_score$Type == "Proximal"],n=10))
# top_genes_direct = names(rank_vectors(mc3_score$gene[mc3_score$Type =="Direct"],n=3))
# top_genes = c(top_genes_cluster,top_genes_proxi,top_genes_direct)
gene_dist = data.frame(table(annotated_clusterMUT_h$SiteType,annotated_clusterMUT_h$Gene_Drug))
colnames(gene_dist) = c("Type","Gene","Count")
gene_dist_g = gene_dist[gene_dist$Gene %in% top_genes_cluster,]
gene_dist_g$Count_plot = gene_dist_g$Count
gene_dist_g$Count_plot[gene_dist_g$Count_plot>50] = 50
getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
p = ggplot(gene_dist_g,aes(y = Type, x = Gene))
p = p + geom_tile(aes(fill=Count_plot), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
p = p + geom_text(aes(label = Count), color="black", size=3)
p = p + theme_bw()
p = p + labs(y = "Mutation interaction with phosphosites", x = "Gene")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=14),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = paste("output/top_genes_category_heatmap.pdf",sep="_")
ggsave(fn,h=4, w = 8,useDingbat=F)

p = ggplot(gene_dist_g[gene_dist_g$Type!="None",],aes(x = Gene, y = Count, fill=Type))
#p = p + facet_grid(Type~., space="free")
p = p + geom_bar(stat = "identity") #+ scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
#p = p + geom_text(aes(label = Count), color="black", size=3)
p = p + theme_bw()
p = p + getTypeFillScale()
p = p + labs(title="Co-clustered mutations")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=14),axis.ticks = element_blank())#element_text(colour="black", size=14))
p = p + theme(legend.position = "bottom")
p
fn = paste("output/top_genes_category_bar.pdf",sep="_")
ggsave(fn,h=4, w = 7,useDingbat=F)

#####
cat("##### examine whether clustered phosphosites have more citation ######\n")

PTMcosmo$clustered = F
PTMcosmo$clustered[paste(PTMcosmo$Ensembl.Transcript,PTMcosmo$Mutation_Gene) %in% paste(annotated_cluster$Transcript,annotated_cluster$Mutation_Gene)] = T
table(PTMcosmo$clustered)
wilcox.test(PTMcosmo$n_evidence[PTMcosmo$clustered],PTMcosmo$n_evidence[!PTMcosmo$clustered])

p = ggplot(PTMcosmo,aes(y = n_evidence, x = clustered, fill=clustered))
p = p + geom_violin(alpha=0.5)  + ylim(0,20)#+ scale_y_log10()
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw()
p = p + labs(y = "Number of Evidence for the Phosphosite")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=14),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = paste("output/n_evidence_for_phosphosites.pdf",sep="_")
ggsave(fn,useDingbat=F)

wilcox.test(PTMcosmo$n_literature[PTMcosmo$clustered],PTMcosmo$n_literature[!PTMcosmo$clustered])
p = ggplot(PTMcosmo,aes(y = n_literature, x = clustered, fill=clustered))
p = p + geom_violin(alpha=0.5)  + ylim(0,20)#+ scale_y_log10()
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw()
p = p + labs(y = "Number of Literature for the Phosphosite")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=14),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = paste("output/n_literature_for_phosphosites.pdf",sep="_")
ggsave(fn,useDingbat=F)

#  evaluate mutations
### note: may need to consider only recurrent mutations
mutations$Type = NA
mutations$Type[paste(mutations$Gene1,mutations$Mutation1) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)] = "Clustered"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites$Gene2,phosphosites$Position-1)] = "Proximal"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites$Gene2,phosphosites$Position-2)] = "Proximal"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites$Gene2,phosphosites$Position+1)] = "Proximal"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites$Gene2,phosphosites$Position+2)] = "Proximal"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites$Gene2,phosphosites$Position)] = "Direct"
mutation_type_count = data.frame(table(mutations$Type,mutations$Gene1))

# evaluate phosphosites
phosphosites$Type = NA
phosphosites$Type[paste(phosphosites$Gene2,phosphosites$Site2) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)] = "Clustered"
phosphosites$Type[paste(phosphosites$Gene2,phosphosites$Position) %in% paste(mutations$Gene1,mutations$Position-1)] = "Proximal"
phosphosites$Type[paste(phosphosites$Gene2,phosphosites$Position) %in% paste(mutations$Gene1,mutations$Position-2)] = "Proximal"
phosphosites$Type[paste(phosphosites$Gene2,phosphosites$Position) %in% paste(mutations$Gene1,mutations$Position+1)] = "Proximal"
phosphosites$Type[paste(phosphosites$Gene2,phosphosites$Position) %in% paste(mutations$Gene1,mutations$Position+2)] = "Proximal"
phosphosites$Type[paste(phosphosites$Gene2,phosphosites$Position) %in% paste(mutations$Gene1,mutations$Position)] = "Direct"
phosphosites_type_count = data.frame(table(phosphosites$Type,phosphosites$Gene2))

table(phosphosites$Type)
table(mutations$Type)
top_mut_genes = names(rank_vectors(mutations$Gene1[mutations$Type!="None"]))

### check something odd about the direct number being different ###
# # likely due to this:
# Gene1    Chromosome1    Start1    Stop1    Mutation1    Chain1    Position1    Feature1    COSMIC1    Gene2    Transcript2    TranscriptPosition2    Site2    Chain2    Position2    Feature2    COSMIC2    LinearDistance
# MKNK2    19    2042789    2042789    p.A192T    [A]    192    Protein    kinase.    {ECO:0000255|PROSITE    ENST00000250896    220    p.S220    [A]    220    Phosphorylation    N/A    28
#   MKNK2    19    2043146    2043146    p.L157Q    [A]    204    Protein    kinase.    {ECO:0000255|PROSITE    ENST00000250896    220    p.E220    [A]    267    Phosphorylation    N/A    63

mutation_type_count$type = "mutation"
phosphosites_type_count$type = "phosphosite"

type_count = rbind(mutation_type_count,phosphosites_type_count)
colnames(type_count) = c("Interaction","Gene","Count","Type")

type_count_g = type_count[type_count$Gene %in% top_mut_genes,]
type_count_g$Gene = factor(type_count_g$Gene,levels = top_mut_genes)
type_count_g$Interaction = factor(type_count_g$Interaction,levels = c("Clustered","Proximal","Direct"))

p = ggplot(type_count_g,aes(x = Gene, y = Count, fill=Interaction))
p = p + facet_grid(Type~., scale="free")
p = p + geom_bar(stat = "identity")
p = p + getTypeFillScale()
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5))
p
fn = paste("output/gene_mut_count.pdf",sep="_")
ggsave(fn,h=8,useDingbat=F)

##### cancer sites only #####

phosphosites = site[,c("Gene2","Site2")]
phosphosites = phosphosites[!duplicated(paste(phosphosites$Gene2,phosphosites$Site2)),]
phosphosites_CPTAC = phosphosites[paste(phosphosites$Gene2,phosphosites$Site2) %in% paste(cptac_site$hgnc_symbol,cptac_site$amino_acid_residue),]

annotated_cluster_in_cptac = annotated_cluster$Cluster[paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene) %in%
                                                         paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Site2)]
annotated_cluster_cptac = annotated_cluster[annotated_cluster$Cluster %in% annotated_cluster_in_cptac,]

mutations$Type = NA
mutations$Type[paste(mutations$Gene1,mutations$Mutation1) %in% paste(annotated_cluster_cptac$Gene_Drug,annotated_cluster_cptac$Mutation_Gene)] = "Clustered"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position-1)] = "Proximal"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position-2)] = "Proximal"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position+1)] = "Proximal"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position+2)] = "Proximal"
mutations$Type[paste(mutations$Gene1,mutations$Position) %in% paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position)] = "Direct"
mutation_type_count = data.frame(table(mutations$Type,mutations$Gene1))

# evaluate phosphosites_CPTAC
phosphosites_CPTAC$Type = NA
phosphosites_CPTAC$Type[paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Site2) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)] = "Clustered"
phosphosites_CPTAC$Type[paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position) %in% paste(mutations$Gene1,mutations$Position-1)] = "Proximal"
phosphosites_CPTAC$Type[paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position) %in% paste(mutations$Gene1,mutations$Position-2)] = "Proximal"
phosphosites_CPTAC$Type[paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position) %in% paste(mutations$Gene1,mutations$Position+1)] = "Proximal"
phosphosites_CPTAC$Type[paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position) %in% paste(mutations$Gene1,mutations$Position+2)] = "Proximal"
phosphosites_CPTAC$Type[paste(phosphosites_CPTAC$Gene2,phosphosites_CPTAC$Position) %in% paste(mutations$Gene1,mutations$Position)] = "Direct"
phosphosites_CPTAC_type_count = data.frame(table(phosphosites_CPTAC$Type,phosphosites_CPTAC$Gene2))

table(phosphosites_CPTAC$Type)
table(mutations$Type)
top_mut_genes = names(rank_vectors(mutations$Gene1[mutations$Type!="None"]))

### check something odd about the direct number being different ###
# # likely due to this:
# Gene1    Chromosome1    Start1    Stop1    Mutation1    Chain1    Position1    Feature1    COSMIC1    Gene2    Transcript2    TranscriptPosition2    Site2    Chain2    Position2    Feature2    COSMIC2    LinearDistance
# MKNK2    19    2042789    2042789    p.A192T    [A]    192    Protein    kinase.    {ECO:0000255|PROSITE    ENST00000250896    220    p.S220    [A]    220    Phosphorylation    N/A    28
#   MKNK2    19    2043146    2043146    p.L157Q    [A]    204    Protein    kinase.    {ECO:0000255|PROSITE    ENST00000250896    220    p.E220    [A]    267    Phosphorylation    N/A    63

mutation_type_count$type = "mutation"
phosphosites_CPTAC_type_count$type = "phosphosite"

type_count = rbind(mutation_type_count,phosphosites_CPTAC_type_count)
colnames(type_count) = c("Interaction","Gene","Count","Type")

type_count_g = type_count[type_count$Gene %in% top_mut_genes,]
type_count_g$Gene = factor(type_count_g$Gene,levels = top_mut_genes)
type_count_g$Interaction = factor(type_count_g$Interaction,levels = c("Clustered","Proximal","Direct"))

p = ggplot(type_count_g,aes(x = Gene, y = Count, fill=Interaction))
p = p + facet_grid(Type~., scale="free")
p = p + geom_bar(stat = "identity")
p = p + getTypeFillScale()
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5))
p
fn = paste("output/gene_mut_count_phoCPTAC.pdf",sep="_")
ggsave(fn,h=8,useDingbat=F)

### clustered sites may have higher driver scores ###
# p = ggplot(site_uniq,aes(x = firstDistance))
# p = p + geom_density(alpha=0.2,size=0.5)
# p = p + theme_bw()
# p
# fn = paste("output/pass.fast.MC3_site_distance_density.pdf",sep="_")
# ggsave(fn, useDingbat=F)
# 
# # check whether sites have enrichment in TCGA
# cptac_site_f = "/Users/khuang/Box\ Sync/hotpho_data/input/CPTAC.BRCA.OV.phosphoproteome.phosphosite.itraq.wgene.tsv"
# cptac_site = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = cptac_site_f)
# 
# site_uniq$inCPTAC = site_uniq$GeneSite2 %in% (paste(cptac_site$hgnc_symbol,cptac_site$ensembl_transcript_id,cptac_site$amino_acid_residue))
# 
# annotated_cluster_filt = annotated_cluster[(annotated_cluster$Type=="Hybrid"),]
# annotated_cluster_filt$GeneSite2 = paste(annotated_cluster_filt$Gene_Drug,annotated_cluster_filt$Transcript,annotated_cluster_filt$Mutation_Gene)
# annotated_cluster_filt_brief = annotated_cluster_filt[,c("Closeness_Centrality","Geodesic_From_Centroid","GeneSite2")]
# 
# site_uniq_c = merge(site_uniq,annotated_cluster_filt_brief,by="GeneSite2",all.x=T)
# site_uniq_c$inCluster = !is.na(site_uniq_c$Closeness_Centrality)
# fisher.test(table(site_uniq_c$inCluster,site_uniq_c$inCPTAC))
# 
# p = ggplot(site_uniq_c,aes(fill = inCluster,x=inCPTAC))
# p = p + facet_grid(inCluster~.,drop=T, scale="free")
# p = p + geom_bar()
# p = p + theme_bw()
# p = p + labs(x="In CPTAC cancer sample")
# p
# fn = paste("output/inCluster_vs_inCPTAC.pdf",sep="_")
# ggsave(fn, useDingbat=F)
# 
# # see if cancer sites have higher CC score
# wilcox.test(site_uniq_c$Closeness_Centrality[site_uniq_c$inCPTAC],site_uniq_c$Closeness_Centrality[!site_uniq_c$inCPTAC])
# t.test(site_uniq_c$Closeness_Centrality[site_uniq_c$inCPTAC],site_uniq_c$Closeness_Centrality[!site_uniq_c$inCPTAC])
# 
# p = ggplot(site_uniq_c,aes(x = Closeness_Centrality,fill=inCPTAC))
# p = p + facet_grid(inCPTAC~.)
# p = p + geom_density(alpha=0.5)
# p = p + theme_bw()
# p
# fn = paste("output/closeness_centrality_inHybrid_vs_inCPTAC.pdf",sep="_")
# ggsave(fn, useDingbat=F)
# 
# p = ggplot(site_uniq_c,aes(x = Closeness_Centrality,fill=inCPTAC))
# p = p + facet_grid(inCPTAC~.)
# p = p + geom_density(alpha=0.5)
# p = p + theme_bw() + xlim(0.5,10)
# p
# fn = paste("output/closeness_centrality_inHybrid_vs_inCPTAC_zoom.pdf",sep="_")
# ggsave(fn, useDingbat=F)

### check whether clustered sites have different residue distribution ###
# can check this in each domain as well


# # previous checks: applied to .sites file only
# # > sum(!(site$GeneSite1 %in% site$GeneSite2))
# # [1] 10342
# # > sum(!(site$GeneSite2 %in% site$GeneSite1))
# # [1] 15551
# # okay there are unique entries in each
# 
# site1 = site[,grep("1",colnames(site))]
# site2 = site[,grep("2",colnames(site))]
# colnames(site2) = colnames(site1)
# all_sites = rbind(site1,site2)
# dim(all_sites)
# all_sites_uniq = all_sites[!duplicated(all_sites$GeneSite1),]
# dim(all_sites_uniq)
# # [1] 29625     9
# #some positions still have multiple residues gotta check on this
# all_sites_uniq = all_sites_uniq[!duplicated(paste(all_sites_uniq$Transcript1,all_sites_uniq$TranscriptPosition1)),]
# dim(all_sites_uniq)
# # [1] 26141     9
# all_sites_uniq$residue = gsub("p.([A-Z]).*","\\1",all_sites_uniq$Site1)
# 
# # examine gene distribution
# gene_counts = table(all_sites_uniq$Gene1)[order(table(all_sites_uniq$Gene1),decreasing = T)]
# top_genes = names(gene_counts[1:20])
# all_sites_uniq$Gene1= factor(all_sites_uniq$Gene1,levels=names(gene_counts))
# 
# all_sites_uniq_g = all_sites_uniq[all_sites_uniq$Gene1 %in% top_genes,]
# plot_barplot(all_sites_uniq_g, x_string = "Gene1", fill_string = "residue",
#              fileName = "pass.fast.MC3_site_counts_top20gene.pdf")