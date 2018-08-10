##### analyze_omnipath.R #####
# Kuan-lin Huang @ WashU 2017 Dec
# updated 2018 May
# analyze coclusters in the context of functional validated mutation
# analyze sites next to functional mutations

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")

# omni path
omnipath_f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/OmniPath/ptms.txt"
omnipath = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, file = omnipath_f)
  
omnipath_tabf = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/OmniPath/op_substrate_siteTab_wpid.txt"
omnipath_tab = read.table(header=T, quote = "", sep=",", stringsAsFactors = F, file = omnipath_tabf)

library(biomaRt)
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
mapTab = getBM(attributes = c("hgnc_symbol", "uniprotswissprot"), filters = "uniprotswissprot", values = omnipath$substrate, mart = ensembl, uniqueRows=FALSE)
dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
if (length(dupRows)==0){mapTabNoDup = mapTab} else {mapTabNoDup = mapTab[-dupRows, ]}
colnames(omnipath)[2]="uniprotswissprot"
omnipath_wgene = merge(omnipath, mapTabNoDup, by="uniprotswissprot")
colnames(omnipath_wgene)[which(colnames(omnipath_wgene )== "hgnc_symbol")] = "SUB_GENE"

mapTab2 = getBM(attributes = c("hgnc_symbol", "uniprotswissprot"), filters = "uniprotswissprot", values = omnipath$enzyme, mart = ensembl, uniqueRows=FALSE)
dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab2[,2])))
if (length(dupRows)==0){mapTabNoDup = mapTab2NoDup} else {mapTab2NoDup = mapTab2[-dupRows, ]}
colnames(mapTab2NoDup) = paste("enzyme",colnames(mapTab2NoDup),sep="_")
colnames(omnipath_wgene)[2]="enzyme_uniprotswissprot"
omnipath_wgene2 = merge(omnipath_wgene, mapTab2NoDup, by="enzyme_uniprotswissprot")

omniPathCombined = merge(omnipath_wgene2,omnipath_tab, by=c("SUB_GENE","residue_type","residue_offset"),all.x=T)
omniPathCombined$HGVSg = paste("chr",omniPathCombined$seq_region_name, ":g.",omniPathCombined$start,"_", omniPathCombined$end,sep="")
omniPathCombined = omniPathCombined[!duplicated(paste(omniPathCombined$enzyme_hgnc_symbol,omniPathCombined$HGVSg)),] # can be on two matched_ens_pid
  
# validation result v2
validation_all_f = "validation_data/AllValidationLevels"
validation_all = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, file = validation_all_f)
colnames(validation_all)[4] = "Transcript"
colnames(validation_all)[2] = "Mutation_Gene"
cat("Number of mutations with validation_all",nrow(validation_all),"\n")
table(validation_all$Level<= 6,validation_all$Effect)
## note: we are getting some MET mutations if we use evidence other than MDanderson
validation_all_act = validation_all[validation_all$Effect == "Activating" & validation_all$Level <= 6,]

# functional mutation in clusters
annotated_cluster_h$HGVSg = gsub("\\/.*","",annotated_cluster_h$GenomicPosition)
annotated_cluster_h_score = merge(annotated_cluster_h,validation_all_act,by=c("Transcript","Mutation_Gene"),all.x=T)
annotated_cluster_h_score_reg = merge(annotated_cluster_h_score, omniPathCombined,by=c("HGVSg"),all.x=T)
table(annotated_cluster_h_score_reg$enzyme_hgnc_symbol)

activating_clusters = annotated_cluster_h_score$Cluster[!is.na(annotated_cluster_h_score$Effect) & annotated_cluster_h_score$Effect=="Activating"]
annotated_cluster_h_score_reg$activating_cluster = annotated_cluster_h_score_reg$Cluster %in% activating_clusters
reg_clusters = annotated_cluster_h_score_reg$Cluster[!is.na(annotated_cluster_h_score_reg$enzyme_hgnc_symbol)]
annotated_cluster_h_score_reg$regulated_cluster = annotated_cluster_h_score_reg$Cluster %in% reg_clusters
annotated_cluster_h_score_reg_plot = annotated_cluster_h_score_reg[annotated_cluster_h_score_reg$regulated_cluster,]

cat("Number of activating clusters",length(unique(activating_clusters)),"\n")
cat("Number of regulated clusters",length(unique(reg_clusters)),"\n")
cat("Number of regulated and activating clusters",length(intersect(unique(activating_clusters),unique(reg_clusters))),"\n")
cat("Number of regulated coclustered phosphosites",sum(!is.na(annotated_cluster_h_score_reg$enzyme_hgnc_symbol) & annotated_cluster_h_score_reg$Alternate == "ptm" & !duplicated(annotated_cluster_h_score_reg$Mutation_Gene)),"\n")
cat("Number of unique genes of regulated coclustered phosphosites",length(unique(annotated_cluster_h_score_reg_plot$Gene_Drug[!is.na(annotated_cluster_h_score_reg_plot$enzyme_hgnc_symbol)])),"\n")
cat("Number of regulated coclustered phosphosites by substrate gene\n")
table(annotated_cluster_h_score_reg_plot$Gene_Drug[!duplicated(annotated_cluster_h_score_reg_plot$Mutation_Gene) & !is.na(annotated_cluster_h_score_reg_plot$enzyme_hgnc_symbol) & annotated_cluster_h_score_reg_plot$Alternate == "ptm"])[order(table(annotated_cluster_h_score_reg_plot$Gene_Drug[!is.na(annotated_cluster_h_score_reg_plot$enzyme_hgnc_symbol) & annotated_cluster_h_score_reg_plot$Alternate == "ptm"]))]
cat("Number of regulated coclustered phosphosites by enzyme gene\n")
table(annotated_cluster_h_score_reg_plot$enzyme_hgnc_symbol[!is.na(annotated_cluster_h_score_reg_plot$enzyme_hgnc_symbol) & annotated_cluster_h_score_reg_plot$Alternate == "ptm"])[order(table(annotated_cluster_h_score_reg_plot$enzyme_hgnc_symbol[!is.na(annotated_cluster_h_score_reg_plot$enzyme_hgnc_symbol) & annotated_cluster_h_score_reg_plot$Alternate == "ptm"]))]

cat("Fisher's test for enrichment of regulated clusters in activated clusters\n")
fisher_elements = c(sum(!unique(annotated_cluster_h$Cluster) %in% c(unique(reg_clusters),unique(activating_clusters))),sum((unique(annotated_cluster_h$Cluster) %in% unique(reg_clusters)) & !(unique(annotated_cluster_h$Cluster) %in% unique(activating_clusters))),
                    sum(!(unique(annotated_cluster_h$Cluster) %in% unique(reg_clusters)) & (unique(annotated_cluster_h$Cluster) %in% unique(activating_clusters))),sum((unique(annotated_cluster_h$Cluster) %in% unique(reg_clusters)) & (unique(annotated_cluster_h$Cluster) %in% unique(activating_clusters))))
test.table = matrix(as.numeric(fisher_elements), nrow=2)
fisher.test(test.table, alternative = "greater")
  

exclude_gene = c("TRIO","TP53BP1","FN1","MYH9","APC","PREX1","ERRFI1","S100A4")
p = ggplot(annotated_cluster_h_score_reg_plot[!(annotated_cluster_h_score_reg_plot$Gene_Drug %in% exclude_gene),],aes(x = as.numeric(Position), y = paste(Gene_Drug,Cluster), color=Alternate=="ptm"))
p = p + facet_grid(Gene_Drug~.,drop = T, scale = "free", space = "free")
p = p + geom_point(alpha=0.3,stroke=0,size=2)# + coord_flip()
p = p + geom_text_repel(aes(label=ifelse(Alternate=="ptm" & !is.na(enzyme_hgnc_symbol) & !duplicated(paste(enzyme_hgnc_symbol,Mutation_Gene)),paste(enzyme_hgnc_symbol,gsub("p.","",Mutation_Gene),sep=" > "),NA)))
p = p + theme_nogrid() + geom_line(color="black",alpha=0.5)
p = p + labs(y="Cluster", x = "Protein coordinate") #+ xlim(0,2000)
p
fn = paste("output/regulated_cluster_by_coordinate.pdf",sep="_")
ggsave(fn, w=10,h=8,useDingbat=F)

write.table(file = "output/annotated_clusters_w_functional_mutations_regulation.tsv", annotated_cluster_h_score_reg[annotated_cluster_h_score_reg$regulated_cluster,],quote=F, sep = '\t',row.names = F)

