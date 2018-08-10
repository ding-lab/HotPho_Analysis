##### analyze_PFAM.R #####
# Kuan-lin Huang @ WashU 2017 July
# updated 2017 Dec
# updated 2018 April

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")

# PFAM file
PFAM_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/data/pdb_pfam_mapping.txt.gz"

##### CLUSTERs #####

annotated_cluster$Gene_site = paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)

PFAM = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, file = gzfile(PFAM_f))
site_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201807/3D_Proximity_cleaned.musites.gz"
site = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(site_f))
site_uniq = site[!duplicated(paste(site$Gene1,site$Mutation1,site$Transcript2,site$TranscriptPosition2)),]
site_uniq$PDB_ID =gsub(".* (.*) .*","\\1",site_uniq$DistanceInfo)

### annotate pFAM name to musite file
site_uniq$PFAM_Name1 = NA
site_uniq$PFAM_Name2 = NA
site_uniq$PFAM_desc1 = NA
site_uniq$PFAM_desc2 = NA
for (i in 1:nrow(site_uniq)){
  PDB_ID = site_uniq$PDB_ID[i]
  pfam = PFAM[PFAM$PDB_ID == PDB_ID,]
  if (nrow(pfam) > 0){
    for (k in 1:nrow(pfam)){
      if (site_uniq$Position1[i] > pfam$PdbResNumStart[k] & site_uniq$Position1[i] < pfam$PdbResNumEnd[k]) {
        site_uniq$PFAM_Name1[i] = pfam$PFAM_Name[k]
        site_uniq$PFAM_desc1[i] = pfam$PFAM_desc[k]
      }
      if (site_uniq$Position2[i] > pfam$PdbResNumStart[k] & site_uniq$Position2[i] < pfam$PdbResNumEnd[k]) {
        site_uniq$PFAM_Name2[i] = pfam$PFAM_Name[k]
        site_uniq$PFAM_desc2[i] = pfam$PFAM_desc[k]
      }
    }
  }
}

mutPFAM1 = site_uniq[,which(colnames(site_uniq) %in% c("Gene1","Mutation1","PFAM_Name1","PFAM_desc1"))]
sitePFAM1 = site_uniq[,which(colnames(site_uniq) %in% c("Gene2","Site2","PFAM_Name2","PFAM_desc2"))]
mutPFAM1$Gene_site = paste(mutPFAM1$Gene1,mutPFAM1$Mutation1)
sitePFAM1$Gene_site = paste(sitePFAM1$Gene2,sitePFAM1$Site2)
colnames(mutPFAM1)=c("Gene","Mutation","PFAM_Name","PFAM_desc","Gene_site")
colnames(sitePFAM1)=c("Gene","Mutation","PFAM_Name","PFAM_desc","Gene_site")
musitePFAM = rbind(mutPFAM1,sitePFAM1) # 11477527  lines
musitePFAM$Gene_site_pfam = paste(musitePFAM$Gene,musitePFAM$Gene_site,musitePFAM$PFAM_Name)
musitePFAM = musitePFAM[!duplicated(musitePFAM$Gene_site_pfam),] # 88029 lines
musitePFAM = musitePFAM[!is.na(musitePFAM$PFAM_Name),] # 31361 lines

annotated_cluster$site_type = "Mutation"
annotated_cluster$site_type[annotated_cluster$Alternate=="ptm"]= "Phosphosite"
annotated_cluster$coord = as.numeric(gsub("p.[A-Z]([0-9]+)[A-Z]*","\\1",annotated_cluster$Mutation_Gene))
annotated_cluster$inCPTAC = FALSE
annotated_cluster$inCPTAC[paste(annotated_cluster$Transcript,annotated_cluster$Mutation_Gene) %in% paste(cptac_site$ensembl_transcript_id,cptac_site$amino_acid_residue)]=TRUE

annotated_cluster_feature = merge(annotated_cluster, musitePFAM, by="Gene_site", all.x=T)
annotated_cluster_feature$site_type = "Mutation"
annotated_cluster_feature$site_type[annotated_cluster_feature$Alternate=="ptm"]= "Phosphosite" 

annotated_cluster_feature_clusterOnly = annotated_cluster_feature[annotated_cluster_feature$Type == "Hybrid",c("Cluster","Gene_Drug","PFAM_Name","PFAM_desc")]
annotated_cluster_feature_clusterOnly_h_u = annotated_cluster_feature_clusterOnly[!is.na(annotated_cluster_feature_clusterOnly$PFAM_desc) & !duplicated(annotated_cluster_feature_clusterOnly),]
PFAMcounts = data.frame(table(annotated_cluster_feature_clusterOnly_h_u$PFAM_desc)[table(annotated_cluster_feature_clusterOnly_h_u$PFAM_desc) > 7])
colnames(PFAMcounts) = c("PFAM_desc","HybridClusterCounts")
PFAMcounts$PFAM_desc = factor(PFAMcounts$PFAM_desc,levels = PFAMcounts$PFAM_desc[order(PFAMcounts$HybridClusterCounts)])

p = ggplot(PFAMcounts,aes(x = PFAM_desc, y=HybridClusterCounts, fill=PFAM_desc))
p = p + geom_bar(stat="identity") + theme_bw()
#p = p + labs(x = "Feature", y="Counts")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + coord_flip() + theme(legend.position = "none")
p
ggsave(file="output/top_PFAM_features_in_hybrid_clusters.pdf", w = 8,useDingbats=FALSE)

annotated_cluster_feature_sele = annotated_cluster_feature[annotated_cluster_feature$PFAM_desc %in% PFAMcounts$PFAM_desc,]
annotated_cluster_feature_sele_mut = annotated_cluster_feature_sele[annotated_cluster_feature_sele$site_type=="Mutation",]
annotated_cluster_feature_sele_pho = annotated_cluster_feature_sele[annotated_cluster_feature_sele$site_type=="Phosphosite",]
annotated_cluster_feature_sele_MutCount = data.frame(table(annotated_cluster_feature_sele_mut$Gene_Drug,annotated_cluster_feature_sele_mut$PFAM_desc))
annotated_cluster_feature_sele_PhoCount = data.frame(table(annotated_cluster_feature_sele_pho$Gene_Drug,annotated_cluster_feature_sele_pho$PFAM_desc))
colnames(annotated_cluster_feature_sele_MutCount) = c("Gene","PFAM_desc","MutCount")
colnames(annotated_cluster_feature_sele_PhoCount) = c("Gene","PFAM_desc","PhoCount")
annotated_cluster_feature_sele_count = merge(annotated_cluster_feature_sele_MutCount,annotated_cluster_feature_sele_PhoCount,by=c("Gene","PFAM_desc"))
annotated_cluster_feature_sele_count = annotated_cluster_feature_sele_count[annotated_cluster_feature_sele_count$MutCount != 0 | annotated_cluster_feature_sele_count$PhoCount != 0,]

annotated_cluster_feature_sele_count$PFAM_desc = factor(annotated_cluster_feature_sele_count$PFAM_desc,levels = PFAMcounts$PFAM_desc[order(PFAMcounts$HybridClusterCounts)])

# annotated_cluster_feature_sele_count = data.frame(table(annotated_cluster_feature_sele$Gene_Drug,annotated_cluster_feature_sele$PFAM_desc,annotated_cluster_feature_sele$site_type))
# colnames(annotated_cluster_feature_sele_count) = c("Gene","PFAM_desc","SiteType","Count")
p = ggplot(annotated_cluster_feature_sele_count,aes(x = MutCount, y = PhoCount, color=PFAM_desc))#, fill=Var2))
p = p + geom_point(alpha=0.3) + theme_bw() + ylim(0,10)
#p = p + labs(x = "Feature", y="Counts") 
p = p + geom_text_repel(aes(label= ifelse(MutCount > 9 | PhoCount > 3,as.character(Gene),NA)))
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
#p = p + theme(legend.position = "bottom")
p
ggsave(file="output/top_PFAM_features_in_hybrid_clusters_gene_site_counts.pdf", w=10,useDingbats=FALSE)


##### look into some specific gene families #####
# histone
histones = read.table("input/histones.tsv",header=F,sep="\t")
histone_genes = histones[,1]
annotated_cluster_h = annotated_cluster[annotated_cluster$Gene_Drug %in% histone_genes,]
annotated_cluster_h$SubFamily = gsub("[A-Z]+$","",annotated_cluster_h$Gene_Drug) # subfamily for histones:https://en.wikipedia.org/wiki/Histone#Actively_transcribed_genes
annotated_cluster_h$SubFamily[annotated_cluster_h$SubFamily=="H2"] = "H2AF"  
annotated_cluster_h_hybrid = annotated_cluster_h[annotated_cluster_h$Type=="Hybrid",]

annotated_cluster_feature_h = annotated_cluster_feature[annotated_cluster_feature$Gene_Drug %in% histone_genes,]
annotated_cluster_feature_h$SubFamily = gsub("[A-Z]+$","",annotated_cluster_feature_h$Gene_Drug) # subfamily for histones:https://en.wikipedia.org/wiki/Histone#Actively_transcribed_genes
annotated_cluster_feature_h$SubFamily[annotated_cluster_feature_h$SubFamily=="H2"] = "H2AF"  
annotated_cluster_feature_h_hybrid = annotated_cluster_feature_h[annotated_cluster_feature_h$Type=="Hybrid",]
annotated_cluster_feature_h_hybrid_p = annotated_cluster_feature_h_hybrid[!is.na(annotated_cluster_feature_h_hybrid$PFAM_desc)
                                                                          & annotated_cluster_feature_h_hybrid$SubFamily %in% c("HIST1H3","HIST1H4"),]
hist_sele_domain = names(table(annotated_cluster_feature_h_hybrid_p$PFAM_desc)[table(annotated_cluster_feature_h_hybrid_p$PFAM_desc)>9])
annotated_cluster_feature_h_hybrid_p_d = annotated_cluster_feature_h_hybrid_p[annotated_cluster_feature_h_hybrid_p$PFAM_desc %in% hist_sele_domain,]
p = ggplot(annotated_cluster_feature_h_hybrid_p_d,aes(x = coord, y=Gene_Drug, shape=site_type, color = PFAM_desc, fill = PFAM_desc))
p = p + facet_grid(SubFamily~.,scale="free",space="free")
p = p + geom_jitter(alpha=0.3,width=0, height=0.2) + theme_bw()
p = p + geom_text_repel(aes(label=ifelse(inCPTAC & !duplicated(Mutation_Gene),Mutation_Gene,NA)))
p = p + labs(x = "Protein coordinate", y="Gene")
#p = p + xlim(0,1250)
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))
p
fn = paste("output/histone_domain_sites_mut_site_linear.plot.pdf",sep=".")
ggsave(file=fn, w=10, useDingbats=FALSE)

# kinase
annotated_cluster_feature_k_cluster = annotated_cluster_feature$Cluster[!is.na(annotated_cluster_feature$PFAM_desc) & 
                                     annotated_cluster_feature$PFAM_desc=="Protein tyrosine kinase" | annotated_cluster_feature$PFAM_desc=="Protein kinase domain"]
annotated_cluster_feature_k = annotated_cluster_feature[annotated_cluster_feature$Cluster %in% annotated_cluster_feature_k_cluster,]
annotated_cluster_feature_k_hybrid = annotated_cluster_feature_k[annotated_cluster_feature_k$Type=="Hybrid" &
                                                                   annotated_cluster_feature_k$Gene_Drug %in% cancer_gene,]
colnames(manning_kinome_wgene_map)[4] = "Gene_Drug"
annotated_cluster_feature_k_hybrid_m = merge(annotated_cluster_feature_k_hybrid,manning_kinome_wgene_map,by="Gene_Drug",all.x=T)
annotated_cluster_feature_k_hybrid_m$PFAM_desc[!(annotated_cluster_feature_k_hybrid_m$PFAM_desc %in% c("Protein tyrosine kinase","Protein kinase domain"))] = "Other"
p = ggplot(annotated_cluster_feature_k_hybrid_m,aes(x = coord, y=Gene_Drug, shape=site_type, color = PFAM_desc))
#p = p + facet_grid(Family~.,scale="free",space="free")
p = p + facet_grid(GroupName~.,scale="free",space="free")
p = p + geom_point(alpha=0.3) + theme_bw()
p = p + geom_text_repel(aes(label=ifelse(inCPTAC & !duplicated(Mutation_Gene),Mutation_Gene,NA)))
p = p + labs(x = "Protein coordinate", y="Gene")
p = p + xlim(0,1250)
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p
ggsave(file="output/kinase_domain_sites_mut_site_linear.plot.pdf", h = 5, useDingbats=FALSE)

##### find top PFAM in each type of cluster #####
for (type in unique(annotated_cluster_feature$Type)) {
  annotated_cluster_feature_t = annotated_cluster_feature[annotated_cluster_feature$Type==type,]
  top_feature = names(table(annotated_cluster_feature_t$PFAM_desc)[order(table(annotated_cluster_feature_t$PFAM_desc),decreasing = T)][2:21])
  annotated_cluster_feature_t_top = annotated_cluster_feature_t[annotated_cluster_feature_t$PFAM_desc %in% top_feature,]
  top_feature_count = data.frame(table(annotated_cluster_feature_t_top$PFAM_desc,annotated_cluster_feature_t_top$site_type))
  
  p = ggplot(top_feature_count,aes(x = Var1, y=Freq, fill=Var2))
  p = p + geom_bar(stat="identity") + theme_bw()
  p = p + labs(x = "Feature", y="Counts")
  p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
  p
  ggsave(file=paste("output/",type,"_top20_PFAM_features.pdf",sep=""), useDingbats=FALSE)
  
  top_gene_in_top_20_feature = names(table(annotated_cluster_feature_t_top$Gene_Drug)[order(table(annotated_cluster_feature_t_top$Gene_Drug),decreasing = T)][1:20])
  annotated_cluster_feature_t_top_topG = annotated_cluster_feature_t_top[annotated_cluster_feature_t_top$Gene_Drug %in% top_gene_in_top_20_feature,]
  top_feature_gene_count = data.frame(table(annotated_cluster_feature_t_top_topG$PFAM_desc,annotated_cluster_feature_t_top_topG$Gene_Drug,annotated_cluster_feature_t_top_topG$site_type))
  top_feature_gene_count$Freq_plot = top_feature_gene_count$Freq
  top_feature_gene_count$Freq_plot[top_feature_gene_count$Freq_plot>20]=20
  
  p = ggplot(top_feature_gene_count,aes(x = Var1, size=Freq_plot, y=Var2, color=Var1))
  p = p + facet_grid(.~Var3)
  p = p + geom_point(stroke = 0) + theme_bw() + scale_size_area()#range = c(0,10))
  p = p + labs(x = "Feature", y="Gene")
  p = p + geom_text(aes(label=ifelse(Freq > 10, Freq, NA)), color="black",size=3)
  p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
  p
  ggsave(file=paste("output/",type,"_top20_PFAM_features_20gene.pdf",sep=""), height=7, width=9, useDingbats=FALSE)
}

##### whether the domains are enriched in significant clusters #####
all_domains_test_enrich=vector("list")
k = 1
annotated_cluster_feature$hybrid=F
annotated_cluster_feature$hybrid[annotated_cluster_feature$Type=="Hybrid"]=T
annotated_cluster_feature_h = annotated_cluster_feature[annotated_cluster_feature$Type=="Hybrid",]
annotated_cluster_feature_h_ptm = annotated_cluster_feature_h[annotated_cluster_feature_h$Alternate=="ptm",]
#test_domain_scores = function(domain){
for (domain in (unique(annotated_cluster_feature$PFAM_desc))){
  #domain_name = annotated_cluster_feature$name[annotated_cluster_feature$domain==domain][1]
  if (sum(annotated_cluster_feature$PFAM_desc==domain, na.rm=T)<5){next}
  # # compare to both mutations/sites in other type of clusters
  # test.table = table(annotated_cluster_feature$PFAM_desc==domain,annotated_cluster_feature$hybrid)
  # NotDomainNotHybrid = test.table[1,1] 
  # DomainNotHybrid = test.table[2,1]
  # NotDomainHybrid = test.table[1,2]
  # DomainHybrid = test.table[2,2]
  # compare to clustered phosphosites to mapped sites
  NotDomainNotHybrid = sum(sitePFAM1$PFAM_desc != domain,na.rm=T) + sum(is.na(sitePFAM1$PFAM_desc))
  DomainNotHybrid = sum(sitePFAM1$PFAM_desc == domain,na.rm=T) - sum(annotated_cluster_feature_h_ptm$PFAM_desc==domain,na.rm=T)
  NotDomainHybrid = sum(annotated_cluster_feature_h_ptm$PFAM_desc!=domain,na.rm=T) + sum(is.na(annotated_cluster_feature_h_ptm$PFAM_desc))
  DomainHybrid = sum(annotated_cluster_feature_h_ptm$PFAM_desc==domain,na.rm=T)
  test.table = matrix(nrow=2,c(NotDomainNotHybrid,DomainNotHybrid,NotDomainHybrid,DomainHybrid))
  # one-tailed fisher's exact test
  P = NA; OR = NA
  f.test = fisher.test(test.table,alternative ="greater")
  OR = f.test$estimate
  P = f.test$p.value
  
  row = c(domain, NotDomainNotHybrid, DomainNotHybrid, NotDomainHybrid, DomainHybrid, OR, P)
  all_domains_test_enrich[[k]] = row
  k = k + 1
}
all_domains_test_enrich_m = do.call(rbind,all_domains_test_enrich)
all_domains_test_enrich_m = data.frame(all_domains_test_enrich_m)
colnames(all_domains_test_enrich_m) = c("domain", "numSite_NotDomainNotHybrid", "numSite_DomainNotHybrid", "numSite_NotDomainHybrid", "numSite_DomainHybrid", "OR", "P")
for (i in 2:7){
  all_domains_test_enrich_m[,i] = as.numeric(as.character(all_domains_test_enrich_m[,i]))
}
all_domains_test_enrich_m$FDR=p.adjust(all_domains_test_enrich_m$P, method="BH")
all_domains_test_enrich_m = all_domains_test_enrich_m[order(as.numeric(as.character(all_domains_test_enrich_m$P)), decreasing=FALSE),]
tn = "output/hybrid_cluster_phosites_pfam_domain_enrichment_fisher.tsv"
write.table(all_domains_test_enrich_m, file=tn, quote=F, sep = '\t', row.names=F)

##### plotting the top features against the top genes #####


top_feature = all_domains_test_enrich_m$domain[!is.na(all_domains_test_enrich_m$domain) & all_domains_test_enrich_m$domain != "N/A"][1:20]
all_domains_m_hybrid_uniprot_top = annotated_cluster_feature[annotated_cluster_feature$PFAM_desc %in% top_feature & annotated_cluster_feature$Type=="Hybrid",]
top_gene_in_top_feature = names(table(all_domains_m_hybrid_uniprot_top$Gene_Drug)[order(table(all_domains_m_hybrid_uniprot_top$Gene_Drug),decreasing = T)][1:20])
all_domains_m_hybrid_uniprot_top_topG = all_domains_m_hybrid_uniprot_top[all_domains_m_hybrid_uniprot_top$Gene_Drug %in% top_gene_in_top_feature,] #check vector etc
top_feature_gene_count_uniprot = data.frame(table(all_domains_m_hybrid_uniprot_top_topG$PFAM_desc,all_domains_m_hybrid_uniprot_top_topG$Gene_Drug,all_domains_m_hybrid_uniprot_top_topG$site_type))

top_feature_gene_count_uniprot$Var1 = factor(top_feature_gene_count_uniprot$Var1,levels = names(table(all_domains_m_hybrid_uniprot_top$PFAM_desc)[order(table(all_domains_m_hybrid_uniprot_top$PFAM_desc))]))
p = ggplot(top_feature_gene_count_uniprot,aes(x = Var1, y=Freq, fill = Var2))
p = p + facet_grid(.~Var3, scale="free", space="free")
p = p + geom_bar(stat = "identity")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme_nogrid()+ theme(legend.position = "bottom") 
p = p + coord_flip()  + labs(x="PFAM domain",y="Number of sites")
p
ggsave(file="output/top20sig_fisher_pfam_domain_features_20_gene_gene_panel_barplot.pdf", h=4,w=7,useDingbats=FALSE)

top_feature_gene_count_uniprot$Freq_plot = top_feature_gene_count_uniprot$Freq
top_feature_gene_count_uniprot$Freq_plot[top_feature_gene_count_uniprot$Freq_plot>20]=20

#### show how many clusters they are from
gene_count = data.frame(table(all_domains_m_hybrid_uniprot_top_topG$Gene_Drug[!duplicated(paste(all_domains_m_hybrid_uniprot_top_topG$Gene_Drug,all_domains_m_hybrid_uniprot_top_topG$Cluster))]))
p = ggplot(gene_count,aes(x = Var1, y=Freq))
p = p + geom_bar(stat = "identity")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "none") 
p = p + coord_flip() + theme_nogrid() + labs(x="Gene",y="Number of clusters in domains")
p
ggsave(file="output/top20sig_fisher_pfam_domain_features_20_gene_gene_panel.pdf", h=4,w=3,useDingbats=FALSE)

p = ggplot(top_feature_gene_count_uniprot,aes(x = Var1, size=Freq_plot, y=Var2, color=Var1))
p = p + facet_grid(.~Var3)
p = p + geom_point(stroke = 0, alpha=0.7) + theme_bw() + scale_size_area()
p = p + labs(x = "Feature", y="Gene")
p = p + geom_text(aes(label=ifelse(Freq > 0, Freq, NA)), color="black",size=3)
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "none")
p
ggsave(file="output/top20sig_fisher_pfam_domain_features_20_gene.pdf", h=7,w=8,useDingbats=FALSE)

getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
p = ggplot(top_feature_gene_count_uniprot,aes(x = Var1, fill=Freq_plot, y=Var2, color=Var1))
p = p + facet_grid(.~Var3)
p = p + geom_tile(linetype="blank")
p = p + geom_text(aes(label = ifelse(Freq!=0,Freq,NA)), color="black", size=2)
p = p + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
p = p + labs(x = "Feature", y="Gene")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "none")
p
ggsave(file="output/top20sig_fisher_pfam_domain_features_20_gene_heatmap.pdf", h=7,w=8,useDingbats=FALSE)


### non histone domains ###
all_domains_test_enrich_m = all_domains_test_enrich_m[!(all_domains_test_enrich_m$domain %in% hist_sele_domain),]

top_feature = all_domains_test_enrich_m$domain[!is.na(all_domains_test_enrich_m$domain) & all_domains_test_enrich_m$domain != "N/A"][1:20]
all_domains_m_hybrid_uniprot_top = annotated_cluster_feature[annotated_cluster_feature$PFAM_desc %in% top_feature & annotated_cluster_feature$Type=="Hybrid",]

top_gene_in_top_feature = names(table(all_domains_m_hybrid_uniprot_top$Gene_Drug)[order(table(all_domains_m_hybrid_uniprot_top$Gene_Drug),decreasing = T)][1:20])
all_domains_m_hybrid_uniprot_top_topG = all_domains_m_hybrid_uniprot_top[all_domains_m_hybrid_uniprot_top$Gene_Drug %in% top_gene_in_top_feature,] #check vector etc
top_feature_gene_count_uniprot = data.frame(table(all_domains_m_hybrid_uniprot_top_topG$PFAM_desc,all_domains_m_hybrid_uniprot_top_topG$Gene_Drug,all_domains_m_hybrid_uniprot_top_topG$site_type))

top_feature_gene_count_uniprot$Var1 = factor(top_feature_gene_count_uniprot$Var1,levels = names(table(all_domains_m_hybrid_uniprot_top$PFAM_desc)[order(table(all_domains_m_hybrid_uniprot_top$PFAM_desc))]))
p = ggplot(top_feature_gene_count_uniprot,aes(x = Var1, y=Freq, fill = Var2))
p = p + facet_grid(.~Var3, scale="free", space="free")
p = p + geom_bar(stat = "identity")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme_nogrid()+ theme(legend.position = "bottom") 
p = p + coord_flip()  + labs(x="PFAM domain",y="Number of sites")
p
ggsave(file="output/top20sig_fisher_pfam_domain_features_20_gene_gene_panel_nonhist_barplot.pdf", h=4,w=7,useDingbats=FALSE)

top_feature_gene_count_uniprot$Freq_plot = top_feature_gene_count_uniprot$Freq
top_feature_gene_count_uniprot$Freq_plot[top_feature_gene_count_uniprot$Freq_plot>20]=20

#### show how many clusters they are from
gene_count = data.frame(table(all_domains_m_hybrid_uniprot_top_topG$Gene_Drug[!duplicated(paste(all_domains_m_hybrid_uniprot_top_topG$Gene_Drug,all_domains_m_hybrid_uniprot_top_topG$Cluster))]))
p = ggplot(gene_count,aes(x = Var1, y=Freq))
p = p + geom_bar(stat = "identity")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "none") 
p = p + coord_flip() + theme_nogrid() + labs(x="Gene",y="Number of clusters in domains")
p
ggsave(file="output/top20sig_fisher_pfam_domain_features_20_gene_gene_panel_nonhist.pdf", h=4,w=3,useDingbats=FALSE)

p = ggplot(top_feature_gene_count_uniprot,aes(x = Var1, size=Freq_plot, y=Var2, color=Var1))
p = p + facet_grid(.~Var3)
p = p + geom_point(stroke = 0, alpha=0.7) + theme_bw() + scale_size_area()
p = p + labs(x = "Feature", y="Gene")
p = p + geom_text(aes(label=ifelse(Freq > 0, Freq, NA)), color="black",size=3)
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "none")
p
ggsave(file="output/top20sig_fisher_pfam_domain_features_20_gene_nonhist.pdf", h=7,w=8,useDingbats=FALSE)

getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
p = ggplot(top_feature_gene_count_uniprot,aes(x = Var1, fill=Freq_plot, y=Var2, color=Var1))
p = p + facet_grid(.~Var3)
p = p + geom_tile(linetype="blank")
p = p + geom_text(aes(label = ifelse(Freq!=0,Freq,NA)), color="black", size=2)
p = p + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
p = p + labs(x = "Feature", y="Gene")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position = "none")
p
ggsave(file="output/top20sig_fisher_pfam_domain_features_20_gene_heatmap_nonhist.pdf", h=7,w=8,useDingbats=FALSE)


