##### analyze_validation.R #####
# Kuan-lin Huang @ WashU 2017 Dec
# updated 2018 May
# analyze coclusters in the context of functional validated mutation
# analyze sites next to functional mutations

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")

# # validation result
# validation_f = "Validation_data/validatedMutations_050918"
# validation = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, file = validation_f)
# colnames(validation)[4] = "Transcript"
# colnames(validation)[2] = "Mutation_Gene"
# cat("Number of mutations with validation",nrow(validation),"\n")
# table(validation$consensus_call)

# PhosphositePlus for cancer sites
hugo = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", fill=T,file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/reference_files/hugo_IDs.txt")
colnames(hugo)[2] = "hgnc_symbol"
colnames(hugo)[18] = "uniprotswissprot"
hugo_sele = hugo[,c("hgnc_symbol","uniprotswissprot")]

library(biomaRt)
pps_cancer = read.table(row.names=NULL,header=TRUE, sep="\t", quote = "", file= "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Disease-associated_sites_term_cancer.oma.leuk_wGPos.txt")
colnames(pps_cancer) <- c(colnames(pps_cancer)[-1],"x")
# biomart conversion: modularize this to a function sometime
# ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org") # this is not working
# mapTab = getBM(attributes = c("ensembl_transcript_id", "uniprotswissprot"), filters = "uniprotswissprot", values = pps_cancer$ACC_ID, mart = ensembl, uniqueRows=FALSE)
# ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hmale_gene_ensembl",host="www.ensembl.org") # this is not working
# mapTab = getBM(attributes = c("ensembl_transcript_id", "uniprotswissprot"), filters = "uniprotswissprot", values = pps_cancer$ACC_ID, mart = ensembl, uniqueRows=FALSE)
# dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
# if (length(dupRows)==0){mapTabNoDup = mapTab} else {mapTabNoDup = mapTab[-dupRows, ]}
colnames(pps_cancer)[5]="uniprotswissprot"
# pps_cancer_wgene = merge(pps_cancer, mapTabNoDup, by="uniprotswissprot")

pps_cancer_wgene = merge(pps_cancer, hugo_sele, by="uniprotswissprot")
pps_cancer_wgene$geneSite = paste(pps_cancer_wgene$hgnc_symbol,tolower(gsub("-p","",pps_cancer_wgene$MOD_RSD)),sep=".")
pps_cancer_wgene = pps_cancer_wgene[-grep("p",pps_cancer_wgene$geneSite),]
pps_cancer_wgene = pps_cancer_wgene[-grep("k",pps_cancer_wgene$geneSite),]
pps_cancer_wgene$Mutation_Gene = paste("p.",gsub("-p","",pps_cancer_wgene$MOD_RSD),sep="")

# phosphosites in hybrid clusters
cat("Phosphosites in known cancer associated sites in PhosphositePlus\n")
table(paste(annotated_cluster_h$Gene_Drug,annotated_cluster_h$Mutation_Gene) %in% paste(pps_cancer_wgene$GENE,pps_cancer_wgene$Mutation_Gene))
#annotated_cluster_h[paste(annotated_cluster_h$Gene_Drug,annotated_cluster_h$Mutation_Gene) %in% paste(pps_cancer_wgene$GENE,pps_cancer_wgene$Mutation_Gene),]
pps_highlighted = pps_cancer_wgene[paste(pps_cancer_wgene$GENE,pps_cancer_wgene$Mutation_Gene) %in% paste(annotated_cluster_h$Gene_Drug,annotated_cluster_h$Mutation_Gene),]
write.table(pps_highlighted, col.names=NA, quote=F, sep = '\t', file="output/known_cancer_coclustering_sites.tsv")

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
annotated_cluster_score = merge(annotated_cluster,validation_all_act,by=c("Transcript","Mutation_Gene"),all.x=T)
annotated_cluster_score = annotated_cluster_score[annotated_cluster_score$Type == "Hybrid",]

activating_clusters = annotated_cluster_score$Cluster[!is.na(annotated_cluster_score$Effect) & annotated_cluster_score$Effect=="Activating"]
activating_annotated_cluster_score = annotated_cluster_score[annotated_cluster_score$Cluster %in% activating_clusters,]

cat("Mutation distribution of activating-mutation clusters:","\n")
table(activating_annotated_cluster_score$Effect)

cat(length(unique(activating_annotated_cluster_score$Gene_Drug))," unique genes.","\n")
cat("Gene distribution of activating-mutation clusters:","\n")
table(activating_annotated_cluster_score$Gene_Drug[!duplicated(paste(activating_annotated_cluster_score$Gene_Drug_Drug,activating_annotated_cluster_score$Cluster))])

cat(length(unique(activating_annotated_cluster_score$Cluster)),"Unique activating-mutation clusters:","\n")
table(activating_annotated_cluster_score$Cluster)

write.table(file = "output/annotated_clusters_w_functional_mutations.tsv", activating_annotated_cluster_score,quote=F, sep = '\t',row.names = F)

activating_annotated_cluster_score$site_type = "other mutation"
activating_annotated_cluster_score$site_type[activating_annotated_cluster_score$Alternate=="ptm"] = "phosphosite"
activating_annotated_cluster_score$site_type[!is.na(activating_annotated_cluster_score$Effect) &
                                               activating_annotated_cluster_score$Effect=="Activating"] = "activating mutation"

numOfSites = data.frame(table(activating_annotated_cluster_score$Cluster,activating_annotated_cluster_score$Gene_Drug))
numOfSites = numOfSites[numOfSites$Freq != 0,]
colnames(numOfSites) = c("Cluster","Gene","SiteCount")
numOfSites = numOfSites[order(numOfSites$Cluster),]
numOfSites$ActivatingMutation = NA
numOfSites$Phosphosite = NA
for (i in 1:nrow(numOfSites)){
  cluster = numOfSites$Cluster[i]
  gene = numOfSites$Gene[i]
  clust_all = activating_annotated_cluster_score[activating_annotated_cluster_score$Cluster==cluster & activating_annotated_cluster_score$Gene_Drug == gene,]
  act = paste(sapply(clust_all$Mutation_Gene[clust_all$site_type=="activating mutation"], paste, collapse=":"), collapse="\t")
  pho = paste(sapply(clust_all$Mutation_Gene[clust_all$Alternate=="ptm"], paste, collapse=":"), collapse="\t")
  numOfSites$ActivatingMutation[i] = act
  numOfSites$Phosphosite[i] = pho
}
write.table(file = "output/functional_mutation_cluster_wKeySites.tsv", numOfSites,quote=F, sep = '\t',row.names = T)

cat("Site distribution within hybrid clusters of activating mutations:","\n")
table(activating_annotated_cluster_score$site_type)
activating_clusters_two_mut = activating_clusters[duplicated(activating_clusters)]
activating_annotated_cluster_score_duo = activating_annotated_cluster_score[activating_annotated_cluster_score$Cluster %in% activating_clusters_two_mut,]
# summarize: number of phosphosites and driver mutation and non-driver mutation in each clusters
counts = data.frame(table(activating_annotated_cluster_score$Gene_Drug,activating_annotated_cluster_score$Cluster,activating_annotated_cluster_score$site_type))
total_counts = counts
colnames(total_counts) = c("gene","cluster","site_type","count")
total_counts_clean = total_counts[total_counts$count>0,]
p = ggplot(total_counts_clean,aes(x=paste(gene), y=count,fill=site_type))
p = p + facet_grid(.~cluster,scale="free",space="free",drop=T) 
p = p + geom_bar(stat = "identity") + theme_nogrid() #+scale_size_area() + scale_size_continuous(range = c(0,max(total_counts$count)))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom") #+ ylim(0,5)# note that this will cut off some values
p = p + labs(y= "Count", x= "")
p
outFile = paste("output/coclustered_residues_to_functional_mutation_counts.pdf")
ggsave(file=outFile, w=9, h = 5, useDingbats=FALSE)

## phosphosites next to activating mutations
# plot and save outputs
p = ggplot(activating_annotated_cluster_score,aes(x=as.factor(Gene_Drug), y=Geodesic_From_Centroid,color=site_type))
p = p + facet_grid(.~Cluster,scale="free",space="free")
#p = p + facet_wrap(~Cluster,drop=T)
p = p + geom_point(alpha=0.3, height = 0, width = 0.2) + theme_nogrid()
p = p + geom_text_repel(size=2,aes(label = ifelse(!duplicated(Mutation_Gene) & (site_type=="phosphosite" | Effect=="activating"), Mutation_Gene, NA)),alpha=0.7) + theme_bw() #+ theme_nogrid()
#p = p + geom_text_repel(aes(label = ifelse(!duplicated(Mutation_Gene) & (site_type=="phosphosite"), Mutation_Gene, NA)),alpha=0.7) + theme_bw() #+ theme_nogrid()
#p = p + geom_text(aes(label = ifelse(!duplicated(Mutation_Gene) & (site_type=="phosphosite"), Mutation_Gene, NA)),alpha=0.9, vjust=-0.5) + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom") #+ ylim(0,5)# note that this will cut off some values
p = p + labs(y= "Geodesic From Centroid", x= "Gene")
# p = p + geom_hline(yintercept = 1, alpha=0.3)
p
outFile = paste("output/coclustered_residues_to_functional_mutation.pdf")
ggsave(file=outFile, w=10, h = 6, useDingbats=FALSE)
