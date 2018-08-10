##### plot_cluster.R #####
# Kuan-lin Huang @ WashU 2017 May
# updated 2017 Dec
# updated 2018 April

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")


# basic stats
cat("Number of unique hybrid clusters:",length(unique(annotated_cluster_h$Cluster)),"\n")
cat("Number of unique genes in hybrid clusters:",length(unique(annotated_cluster_h$Gene_Drug)),"\n")
cat("Number of phosphosites in hybrid clusters:",sum(annotated_cluster_h$Alternate=="ptm",na.rm=T),"\n")
cat("Number of mutations in hybrid clusters:",sum(annotated_cluster_h$Alternate!="ptm",na.rm=T),"\n")
cat("Number of unique mutations in hybrid clusters:",sum(annotated_cluster_h$Alternate[!duplicated(paste(annotated_cluster_h$Gene_Drug,annotated_cluster_h$Mutation_Gene))]!="ptm",na.rm=T),"\n")
annotated_cluster$Position = as.numeric(gsub("p.[A-Z]([0-9]+)[A-Z]*","\\1",annotated_cluster$Mutation_Gene))
annotated_clusterMUT = annotated_cluster[annotated_cluster$Alternate!="ptm",]
annotated_clusterPTM = annotated_cluster[annotated_cluster$Alternate=="ptm",]
annotated_cluster_centroids = annotated_cluster[annotated_cluster$Geodesic_From_Centroid==0,]
annotated_cluster_centroids_unique = annotated_cluster_centroids[!duplicated(annotated_cluster_centroids$Cluster),]

annotated_clusterPTM$ref = gsub("p.([A-Z])[0-9]+","\\1",annotated_clusterPTM$Mutation_Gene)
annotated_clusterPTM_h = annotated_clusterPTM[annotated_clusterPTM$Type=="Hybrid",]
annotated_clusterMUT_h = annotated_clusterMUT[annotated_clusterMUT$Type=="Hybrid",]
annotated_clusterPTM_h$residue = annotated_clusterPTM_h$ref
annotated_clusterPTM_h$residue[!(annotated_clusterPTM_h$residue %in% c("S","T","Y"))]= "Other"

# cluster sites also as proximal or directly mutated by mutations
annotated_clusterPTM_h$SiteType = "Clustered"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position-1)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position-2)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position+1)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position+2)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position)] = "Direct"
cat("Clustered phosphosites showing up as direct or proximal:\n")
table(annotated_clusterPTM_h$SiteType)

# # add transvar ones 
# transvar_f = "HotSpot3D/Data_201807/PTM_Site_Transvar.txt.gz"
# transvar = read.table(header=F, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(transvar_f))
# transvar_anno = transvar[,c(3,6,12)]
# colnames(transvar_anno) = c("Transcript","Position","GenomicPosition")
# annotated_cluster_h = merge(annotated_cluster_h, transvar_anno, by=c("Transcript","Position"), all.x=T)
# annotated_cluster_h$Start[annotated_cluster_h$Alternate == "ptm"] = gsub("chr.*g.([0-9]*)_.*","\\1",annotated_cluster_h$GenomicPosition[annotated_cluster_h$Alternate == "ptm"])
# annotated_cluster_h$Stop[annotated_cluster_h$Alternate == "ptm"] = gsub("chr.*g.([0-9]*)_([0-9]*)/.*","\\2",annotated_cluster_h$GenomicPosition[annotated_cluster_h$Alternate == "ptm"])
# # note for PTM these transvar annotation are HG38
# ### use transvar to translate these back 
# annotated_cluster_h_PTM_transvarIn = annotated_cluster_h[annotated_cluster_h$Alternate == "ptm",c("GenomicPosition")]
# annotated_cluster_h_PTM_transvarIn = paste(gsub("/.*","",annotated_cluster_h_PTM_transvarIn),"delinsAAA",sep="")
# tn = "output/annotated_cluster_h_PTM_transvarIn.txt"
# write.table(annotated_cluster_h_PTM_transvarIn, quote=F, sep="\t", file = tn, row.names = F)

annotated_clusterPTM_h$SiteType = "Clustered"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position-1)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position-2)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position+1)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position+2)] = "Proximal"
annotated_clusterPTM_h$SiteType[paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position) %in% paste(annotated_clusterMUT$Transcript,annotated_clusterMUT$Position)] = "Direct"
cat("Clustered phosphosites showing up as direct or proximal:\n")
table(annotated_clusterPTM_h$SiteType)

annotated_clusterMUT_h$SiteType = "Clustered"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position-1)] = "Proximal"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position-2)] = "Proximal"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position+1)] = "Proximal"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position+2)] = "Proximal"
annotated_clusterMUT_h$SiteType[paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position)] = "Direct"
cat("Clustered mutations showing up as direct or proximal:\n")
table(annotated_clusterMUT_h$SiteType)

##### plotting cluster #####
all_gene_type_count = data.frame(table(annotated_cluster_centroids_unique$Gene_Drug,annotated_cluster_centroids_unique$Type))
p = ggplot(all_gene_type_count,aes(x = Var2, y = Freq, color=Var2))
p = p + geom_point(alpha=0.3)
p = p + geom_text_repel(aes(label=ifelse((Var2=="Hybrid" & Freq > 2) | (Var2=="Mut_Only" & Freq > 8) | (Var2=="Site_Only" & Freq > 1), as.character(Var1), NA)))
p = p + theme_bw()+ scale_size_continuous(range = c(0,4))
p = p + labs(x="Cluster type", y = "Number of clusters")
p
fn = paste("output/Data_201807_p0.05_cluster_count_by_type_gene.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

# examine gene distribution
gene_counts_all = table(annotated_cluster_centroids_unique$Gene_Drug)[order(table(annotated_cluster_centroids_unique$Gene_Drug),decreasing = T)]
gene_counts = table(annotated_cluster_centroids_unique$Gene_Drug[annotated_cluster_centroids_unique$Type=="Hybrid"])[order(table(annotated_cluster_centroids_unique$Gene_Drug[annotated_cluster_centroids_unique$Type=="Hybrid"]),decreasing = T)]
top_genes = c(names(gene_counts[gene_counts>2]))
top_genes = top_genes[top_genes %in% annotated_cluster_centroids_unique$Gene_Drug[annotated_cluster_centroids_unique$Type=="Hybrid"]]
annotated_cluster_centroids_unique_g = annotated_cluster_centroids_unique[annotated_cluster_centroids_unique$Gene_Drug %in% top_genes,]
annotated_cluster_centroids_unique_g$Gene_Drug= factor(annotated_cluster_centroids_unique_g$Gene_Drug,levels=names(gene_counts))
annotated_cluster_centroids_unique_g$Type= factor(annotated_cluster_centroids_unique_g$Type,levels=c("Site_Only","Mut_Only","Hybrid"  ))
annotated_clusterPTM_h_g = annotated_clusterPTM_h[annotated_clusterPTM_h$Gene_Drug %in% top_genes,]
annotated_clusterPTM_h_g$Gene_Drug= factor(annotated_clusterPTM_h_g$Gene_Drug,levels=names(gene_counts))
annotated_clusterMUT_h_g = annotated_clusterMUT_h[annotated_clusterMUT_h$Gene_Drug %in% top_genes,]
annotated_clusterMUT_h_g$Gene_Drug= factor(annotated_clusterMUT_h_g$Gene_Drug,levels=names(gene_counts))

annotated_cluster_centroids_unique_g = annotated_cluster_centroids_unique_g[annotated_cluster_centroids_unique_g$Gene_Drug != "HIST1H4G",]
p = ggplot(annotated_cluster_centroids_unique_g,aes(x = Gene_Drug, fill=Type))
p = p + geom_bar()
p = p + theme_bw() + coord_flip()
p = p + labs(x="Gene", y = "Count of clusters")
p
fn = paste("output/Data_201807_cc.p0.05_gene_counts.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

# top cluster in each gene
cluster_summary = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = paste(cluster_f,"summary",sep="."), colClasses=c(rep('character',2),rep("numeric",13),rep("character",7)))
top_clusters = c()
for (gene in top_genes){
  cluster_summary_gene = cluster_summary[gsub(":.*","",cluster_summary$Centroid) == gene,]
  top_cluster = cluster_summary_gene$Cluster_ID[order(cluster_summary_gene$Cluster_Closeness, decreasing = T)][1]
  top_clusters = c(top_cluster,top_clusters)
}

annotated_cluster_h_top_clusters = annotated_cluster_h[annotated_cluster_h$Cluster %in% top_clusters & annotated_cluster_h$Gene_Drug %in% top_genes & annotated_cluster_h$Gene_Drug != "MGAM",]
annotated_cluster_h_top_clusters$Gene_Drug= factor(annotated_cluster_h_top_clusters$Gene_Drug,levels=names(gene_counts))
p = ggplot(annotated_cluster_h_top_clusters,aes(x = as.numeric(Position), y = paste(Gene_Drug,Cluster), color=Alternate=="ptm"))
p = p + geom_point(alpha=0.3,stroke=0,size=2)# + coord_flip()
p = p + geom_text_repel(aes(label=ifelse(Alternate=="ptm" & !duplicated(Mutation_Gene),gsub("p.","",Mutation_Gene),NA)))
p = p + theme_nogrid() + geom_line(color="black",alpha=0.5)
p = p + labs(y="Cluster", x = "Protein coordinate") #+ xlim(0,2000)
p
fn = paste("output/top_cluster_by_coordinate.pdf",sep="_")
ggsave(fn, w=10,h=6,useDingbat=F)

##### centroid size plot #####
p = ggplot(annotated_cluster_centroids_unique_g,aes(x = Position, y = Gene_Drug, color=Type, size=Total_count))
p = p + geom_point(alpha=0.3,stroke=0)# + coord_flip()
p = p + theme_bw()
p = p + labs(y="Gene", x = "Protein coordinate*") + xlim(0,2000)
p
fn = paste("output/cluster_by_coordinate.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

annotated_cluster_centroids_unique_g_hy = annotated_cluster_centroids_unique_g[annotated_cluster_centroids_unique_g$Type=="Hybrid",]
annotated_cluster_centroids_unique_g_hy_m = melt(annotated_cluster_centroids_unique_g_hy,id.vars = colnames(annotated_cluster_centroids_unique_g_hy)[!(colnames(annotated_cluster_centroids_unique_g_hy) %in% c("Mut_count", "Site_count"))])

feat_genes = annotated_cluster_centroids_unique_g_hy_m$Gene_Drug[annotated_cluster_centroids_unique_g_hy_m$Alternate=="ptm"]
annotated_cluster_centroids_unique_g_hy_m_g = annotated_cluster_centroids_unique_g_hy_m[annotated_cluster_centroids_unique_g_hy_m$Gene_Drug %in% feat_genes,]

p = ggplot(annotated_cluster_centroids_unique_g_hy_m_g,aes(y = Position, x = Gene_Drug, color=variable))
p = p + geom_point(alpha=0.3,stroke=0, aes(size=value)) + coord_flip()
p = p + geom_text_repel(aes(color="centroid",label=ifelse(Alternate=="ptm" & variable=="Site_count", gsub("p.","",Mutation_Gene), NA)))
p = p + theme_bw()#+ scale_size_continuous(range = c(0,4))
p = p + labs(x="Gene", y = "Protein coordinate*") #+ ylim(0,900)
p = p + theme(axis.title = element_text(size=14), axis.text.x = element_text(colour="black", size=12, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p
fn = paste("output/cluster_by_coordinate_hybrid_only.pdf",sep="_")
ggsave(fn, w=8,h=4,useDingbat=F)


##### STY distribution #####
annotated_clusterPTM_dis = data.frame(table(annotated_clusterPTM$ref,annotated_clusterPTM$Type))
annotated_clusterPTM_dis$Var1 = as.character(annotated_clusterPTM_dis$Var1)
annotated_clusterPTM_dis$Var1[!(annotated_clusterPTM_dis$Var1 %in% c("S","T","Y"))]= "Other"

p = ggplot(annotated_clusterPTM_h_g,aes(x = Gene_Drug, fill=residue))
p = p + geom_bar() + coord_flip() + labs(x = "Gene", y = "Count of co-clustered phosphosites")
p = p + theme_bw() #+ scale_y_log10() 
#p = p + geom_vline(xintercept = 1,alpha=0.5)
p
fn = paste("output/Data_201807_cc_filtered_clus_by_residueType_byGene.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

p = ggplot(annotated_clusterPTM_dis,aes(x = Var1, y = Freq, fill=Var2))
p = p + facet_grid(Var2~.)
p = p + geom_bar(stat="identity")
p = p + theme_bw() #+ scale_y_log10()
#p = p + geom_vline(xintercept = 1,alpha=0.5)
p
fn = paste("output/Data_201807_cc_filtered_clus_by_residueType.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

annotated_clusterPTM_dis_h = annotated_clusterPTM_dis[annotated_clusterPTM_dis$Var2=="Hybrid",]
annotated_clusterPTM_dis_h$Pass = "All"
PTM_dis_h = rbind(annotated_clusterPTM_dis_h,annotated_clusterPTM_dis_h)
p = ggplot(PTM_dis_h,aes(x = Var1, y = Freq, fill=Var1))
p = p + facet_grid(Pass~., scale="free")
p = p + geom_bar(stat="identity")
p = p + theme_bw()
p
fn = paste("output/Data_201807_cc_clus_by_residueType.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)


##### plot mutation-phosphosite relation #####

annotated_cluster_centroids_unique_g$Gene_Drug= factor(annotated_cluster_centroids_unique_g$Gene_Drug,levels=names(gene_counts))

annotated_clusterPTM_h_g$SiteType = factor(annotated_clusterPTM_h_g$SiteType,levels = c("None","Clustered","Proximal","Direct"))
p = ggplot(annotated_clusterPTM_h_g,aes(x = Gene_Drug, fill=SiteType))
p = p + geom_bar() + labs(x = "Gene", y = "Count of co-clustered phosphosites", title = "Co-clustered Phosphosites")
p = p + theme_bw() + getTypeFillScale()
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle = 90, vjust=0.5))
p
fn = paste("output/Data_201807_clustered_sites_byInteractionType_byGene.pdf",sep="_")
ggsave(fn, w=8,h=4,useDingbat=F)

annotated_clusterMUT_h_g$SiteType = factor(annotated_clusterMUT_h_g$SiteType,levels = c("None","Clustered","Proximal","Direct"))

annotated_clusterMUT_h_g = annotated_clusterMUT_h_g[annotated_clusterMUT_h_g$Gene_Drug != "HIST1H4G",]
p = ggplot(annotated_clusterMUT_h_g,aes(x = Gene_Drug, fill=SiteType))
p = p + geom_bar() + labs(x = "Gene", y = "Count of co-clustered phosphosites", title = "Co-clustered Mutations")
p = p + theme_bw() + getTypeFillScale() + scale_y_sqrt(breaks = c(0,1, 5, 10,50,100,500,1000))#+ scale_y_log10()
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle = 90, vjust=0.5))
p
fn = paste("output/Data_201807_clustered_mutations_byInteractionType_byGene.pdf",sep="_")
ggsave(fn, w=8,h=4,useDingbat=F)

# MC3 driver score
mc3_score_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/somatic/Driver_BaileyCell2018/Mutation.CTAT.3D.Scores.txt.gz"
mc3_score = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(mc3_score_f))
colnames(mc3_score) = gsub("\\.","",colnames(mc3_score))


##### examine whether clustered mutations have higher mc3 score ###### # move this to plot_cluster.R?
mc3_score$Position = as.numeric(gsub("p.[A-Z](.*)[A-Z]","\\1",mc3_score$protein_change))
mc3_score = mc3_score[mc3_score$gene %in% annotated_cluster$Gene_Drug,]
mc3_score$Type = "None"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$protein_change) %in% paste(annotated_clusterMUT_h$Transcript,annotated_clusterMUT_h$Mutation_Gene)] = "Clustered"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position-1)] = "Proximal"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position-2)] = "Proximal"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position+1)] = "Proximal"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position+2)] = "Proximal"
mc3_score$Type[paste(mc3_score$transcript,mc3_score$Position) %in% paste(annotated_clusterPTM_h$Transcript,annotated_clusterPTM_h$Position)] = "Direct"
cat("MC3 mutations and its relations to phosphosites:\n")
table(mc3_score$Type)

cat("Testing whether overlapping or clustered mutations have higher functional score than others","\n")
wilcox.test(mc3_score$eigenscorefunctional[mc3_score$Type=="Clustered"],mc3_score$eigenscorefunctional[mc3_score$Type=="None"])
wilcox.test(mc3_score$eigenscorefunctional[mc3_score$Type=="Clustered"],mc3_score$eigenscorefunctional[mc3_score$Type=="Proximal"])
wilcox.test(mc3_score$eigenscorefunctional[mc3_score$Type=="Clustered"],mc3_score$eigenscorefunctional[mc3_score$Type=="Direct"])

cat("Testing whether overlapping or clustered mutations have higher cancer-specific functional score than others","\n")
wilcox.test(mc3_score$eigenscorecancer[mc3_score$Type=="Clustered"],mc3_score$eigenscorecancer[mc3_score$Type=="None"])
wilcox.test(mc3_score$eigenscorecancer[mc3_score$Type=="Clustered"],mc3_score$eigenscorecancer[mc3_score$Type=="Proximal"])
wilcox.test(mc3_score$eigenscorecancer[mc3_score$Type=="Clustered"],mc3_score$eigenscorecancer[mc3_score$Type=="Direct"])

mc3_score$Type = factor(mc3_score$Type,levels = c("None","Clustered","Proximal","Direct"))

mc3_score_m = melt(mc3_score[,c(10:13,34,63)], id.var = c("Type"))
colnames(mc3_score_m) = c("InteractionWithPhosphosite","PredictionTool","Score")
mc3_score_m$PredictionTool = gsub("score","",mc3_score_m$PredictionTool)
p = ggplot(mc3_score_m,aes(x = InteractionWithPhosphosite, y = Score, fill=InteractionWithPhosphosite))
p = p + facet_grid(PredictionTool~., scale = "free_y")
p = p + geom_violin(alpha=0.8)
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8)
p = p + getTypeFillScale()
p = p + theme_bw() + theme(legend.position = "none")
p = p + labs(y = "Functional Score", x = "Interaction with Phosphosites")
p
fn = paste("output/all_mutation_functionalscore_vs_interaction_wPhospho.pdf",sep="_")
ggsave(fn,h=8, w = 4,useDingbat=F)
