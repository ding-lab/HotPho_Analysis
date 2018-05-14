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

# plotting cluster
annotated_cluster$coord = as.numeric(gsub("p.[A-Z]([0-9]+)[A-Z]*","\\1",annotated_cluster$Mutation_Gene))

annotated_cluster_centroids = annotated_cluster[annotated_cluster$Geodesic_From_Centroid==0,]
annotated_cluster_centroids_unique = annotated_cluster_centroids[!duplicated(annotated_cluster_centroids$Cluster),]

all_gene_type_count = data.frame(table(annotated_cluster_centroids_unique$Gene_Drug,annotated_cluster_centroids_unique$Type))
p = ggplot(all_gene_type_count,aes(x = Var2, y = Freq, color=Var2))
p = p + geom_point(alpha=0.3)
p = p + geom_text_repel(aes(label=ifelse((Var2=="Hybrid" & Freq > 2) | (Var2=="Mut_Only" & Freq > 8) | (Var2=="Site_Only" & Freq > 1), as.character(Var1), NA)))
p = p + theme_bw()+ scale_size_continuous(range = c(0,4))
p = p + labs(x="Cluster type", y = "Number of clusters")
p
fn = paste("output/Data_201803_p0.05_cluster_count_by_type_gene.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

# examine gene distribution
gene_counts_all = table(annotated_cluster_centroids_unique$Gene_Drug)[order(table(annotated_cluster_centroids_unique$Gene_Drug),decreasing = T)]
gene_counts = table(annotated_cluster_centroids_unique$Gene_Drug[annotated_cluster_centroids_unique$Type=="Hybrid"])[order(table(annotated_cluster_centroids_unique$Gene_Drug[annotated_cluster_centroids_unique$Type=="Hybrid"]),decreasing = T)]
top_genes = c(names(gene_counts[gene_counts>2]),intersect(names(gene_counts[gene_counts>1]),names(gene_counts_all[gene_counts_all>2])))
top_genes = top_genes[top_genes %in% annotated_cluster_centroids_unique$Gene_Drug[annotated_cluster_centroids_unique$Type=="Hybrid"]]
annotated_cluster_centroids_unique_g = annotated_cluster_centroids_unique[annotated_cluster_centroids_unique$Gene_Drug %in% top_genes,]
annotated_cluster_centroids_unique_g$Gene_Drug= factor(annotated_cluster_centroids_unique_g$Gene_Drug,levels=names(gene_counts))
annotated_cluster_centroids_unique_g$Type= factor(annotated_cluster_centroids_unique_g$Type,levels=c("Site_Only","Mut_Only","Hybrid"  ))

# plot_barplot(annotated_cluster_centroids_unique_g, x_string = "Gene_Drug", fill_string = "Type",
#              fileName = "Data_201803_cc.p0.05_gene_counts.pdf")
p = ggplot(annotated_cluster_centroids_unique_g,aes(x = Gene_Drug, fill=Type))
p = p + geom_bar()
p = p + theme_bw() + coord_flip()
p = p + labs(x="Gene", y = "Count of clusters")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p
fn = paste("output/Data_201803_cc.p0.05_gene_counts.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

gene_type_count = data.frame(table(annotated_cluster_centroids_unique_g$Type,annotated_cluster_centroids_unique_g$Gene_Drug))
gene_type_count$cancer_gene = gene_type_count$Var2 %in% cancer_gene
p = ggplot(gene_type_count,aes(x = Var1, y = Var2, color=Var1))
p = p + facet_grid(cancer_gene~., drop=T, scale="free", space="free")
p = p + geom_point(aes(size=Freq))
p = p + theme_bw()+ scale_size_continuous(range = c(0,4))
p = p + labs(x="Cluster type", y = "Gene")
p
fn = paste("output/Data_201803_cc.p0.05_dist_by_cluster_type_bubble.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

p = ggplot(annotated_cluster_centroids_unique_g,aes(x = Mut_count, y = Site_count, color=Type))
#p = p + facet_grid(cancer_gene~., drop=T, scale="free", space="free")
p = p + geom_point(alpha=0.2)
p = p + geom_text_repel(aes(label=ifelse((Type=="Hybrid" & Total_count > 17) & (Site_count > 1 | Total_count > 6), paste(as.character(Gene_Drug),Mutation_Gene), NA)))
p = p + theme_bw()  + ylim(0,7) + xlim(0,13)
p = p + labs(x="Number of mutations", y = "Number of phosphosites")
#p = p + coord_equal()
p
fn = paste("output/Data_201803_cc.p0.05_hybrid_cluster_siteCounts.pdf",sep="_")
ggsave(fn, w=6,h=5,useDingbat=F)

annotated_cluster_centroids_unique_g$Gene_Drug= factor(annotated_cluster_centroids_unique_g$Gene_Drug,levels=names(gene_counts))

annotated_clusterPTM = annotated_cluster[annotated_cluster$Alternate=="ptm",]
annotated_clusterPTM$ref = gsub("p.([A-Z])[0-9]+","\\1",annotated_clusterPTM$Mutation_Gene)
annotated_clusterPTM_g = annotated_clusterPTM[annotated_clusterPTM$Gene_Drug %in% top_genes,]
annotated_clusterPTM_g$Gene_Drug= factor(annotated_clusterPTM_g$Gene_Drug,levels=names(gene_counts))
annotated_clusterPTM_g_h = annotated_clusterPTM_g[annotated_clusterPTM_g$Type=="Hybrid",]
annotated_clusterPTM_g_h$residue = annotated_clusterPTM_g_h$ref
annotated_clusterPTM_g_h$residue[!(annotated_clusterPTM_g_h$residue %in% c("S","T","Y"))]= "Other"

annotated_clusterPTM_dis = data.frame(table(annotated_clusterPTM$ref,annotated_clusterPTM$Type))
annotated_clusterPTM_dis$Var1 = as.character(annotated_clusterPTM_dis$Var1)
annotated_clusterPTM_dis$Var1[!(annotated_clusterPTM_dis$Var1 %in% c("S","T","Y"))]= "Other"

p = ggplot(annotated_clusterPTM_g_h,aes(x = Gene_Drug, fill=residue))
#p = p + facet_grid(Var2~.)
p = p + geom_bar() + coord_flip() + labs(x = "Gene", y = "Count of co-clustered phosphosites")
p = p + theme_bw() #+ scale_y_log10() 
#p = p + geom_vline(xintercept = 1,alpha=0.5)
p
fn = paste("output/Data_201803_cc_filtered_clus_by_residueType_byGene.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

p = ggplot(annotated_clusterPTM_dis,aes(x = Var1, y = Freq, fill=Var2))
p = p + facet_grid(Var2~.)
p = p + geom_bar(stat="identity")
p = p + theme_bw() #+ scale_y_log10()
#p = p + geom_vline(xintercept = 1,alpha=0.5)
p
fn = paste("output/Data_201803_cc_filtered_clus_by_residueType.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

annotated_clusterPTM_dis_h = annotated_clusterPTM_dis[annotated_clusterPTM_dis$Var2=="Hybrid",]
annotated_clusterPTM_dis_h = annotated_clusterPTM_dis[annotated_clusterPTM_dis$Var2=="Hybrid",]
annotated_clusterPTM_dis_h$Pass = "Pass"
annotated_clusterPTM_dis_h$Pass = "All"
PTM_dis_h = rbind(annotated_clusterPTM_dis_h,annotated_clusterPTM_dis_h)
p = ggplot(PTM_dis_h,aes(x = Var1, y = Freq, fill=Var1))
p = p + facet_grid(Pass~., scale="free")
p = p + geom_bar(stat="identity")
p = p + theme_bw() #+ scale_y_log10()
#p = p + geom_vline(xintercept = 1,alpha=0.5)
p
fn = paste("output/Data_201803_cc_clus_by_residueType.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

### centroid size plot ###
p = ggplot(annotated_cluster_centroids_unique_g,aes(x = coord, y = Gene_Drug, color=Type, size=Total_count))
p = p + geom_point(alpha=0.3,stroke=0)# + coord_flip()
#p = p + geom_text_repel(aes(label=ifelse((Var2=="Hybrid" & Freq > 3) | (Var2=="Mut_Only" & Freq > 9) | (Var2=="Site_Only" & Freq > 1), as.character(Var1), NA)))
p = p + theme_bw()#+ scale_size_continuous(range = c(0,4))
p = p + labs(y="Gene", x = "Protein coordinate*") + xlim(0,2000)
p
fn = paste("output/cluster_by_coordinate.pdf",sep="_")
ggsave(fn, w=5,h=5,useDingbat=F)

annotated_cluster_centroids_unique_g_hy = annotated_cluster_centroids_unique_g[annotated_cluster_centroids_unique_g$Type=="Hybrid",]
annotated_cluster_centroids_unique_g_hy_m = melt(annotated_cluster_centroids_unique_g_hy,id.vars = colnames(annotated_cluster_centroids_unique_g_hy)[!(colnames(annotated_cluster_centroids_unique_g_hy) %in% c("Mut_count", "Site_count"))])

feat_genes = annotated_cluster_centroids_unique_g_hy_m$Gene_Drug[annotated_cluster_centroids_unique_g_hy_m$Alternate=="ptm"]
annotated_cluster_centroids_unique_g_hy_m_g = annotated_cluster_centroids_unique_g_hy_m[annotated_cluster_centroids_unique_g_hy_m$Gene_Drug %in% feat_genes,]

p = ggplot(annotated_cluster_centroids_unique_g_hy_m_g,aes(y = coord, x = Gene_Drug, color=variable))
p = p + geom_point(alpha=0.3,stroke=0, aes(size=value)) + coord_flip()
p = p + geom_text_repel(aes(color="centroid",label=ifelse(Alternate=="ptm" & variable=="Site_count", gsub("p.","",Mutation_Gene), NA)))
p = p + theme_bw()#+ scale_size_continuous(range = c(0,4))
p = p + labs(x="Gene", y = "Protein coordinate*") + ylim(0,900)
p = p + theme(axis.title = element_text(size=14), axis.text.x = element_text(colour="black", size=12, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p
fn = paste("output/cluster_by_coordinate_hybrid_only.pdf",sep="_")
ggsave(fn, w=8,h=5,useDingbat=F)

