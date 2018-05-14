##### analyze_validation.R #####
# Kuan-lin Huang @ WashU 2017 Dec
# updated 2018 May
# analyze coclusters in the context of functional validated mutation
# analyze sites next to functional mutations

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")

# validation result
validation_f = "Validation_data/validatedMutations_050918"
validation = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, file = validation_f)
colnames(validation)[4] = "Transcript"
colnames(validation)[2] = "Mutation_Gene"
cat("Number of mutations with validation",nrow(validation),"\n")
table(validation$consensus_call)

# functional mutation in clusters
annotated_cluster_score = merge(annotated_cluster,validation,by=c("Transcript","Mutation_Gene"),all.x=T)
annotated_cluster_score = annotated_cluster_score[annotated_cluster_score$Type == "Hybrid",]

activating_clusters = annotated_cluster_score$Cluster[!is.na(annotated_cluster_score$consensus_call) & annotated_cluster_score$consensus_call=="activating"]
activating_annotated_cluster_score = annotated_cluster_score[annotated_cluster_score$Cluster %in% activating_clusters,]

cat("Mutation distribution of activating-mutation clusters:","\n")
table(activating_annotated_cluster_score$consensus_call)

cat("Gene distribution of activating-mutation clusters:","\n")
table(activating_annotated_cluster_score$Gene_Drug[!duplicated(paste(activating_annotated_cluster_score$Gene_Drug_Drug,activating_annotated_cluster_score$Cluster))])

write.table(file = "output/annotated_clusters_w_functional_mutations.tsv", activating_annotated_cluster_score,quote=F, sep = '\t',row.names = F)

activating_annotated_cluster_score$site_type = "other mutation"
activating_annotated_cluster_score$site_type[activating_annotated_cluster_score$Alternate=="ptm"] = "phosphosite"
activating_annotated_cluster_score$site_type[!is.na(activating_annotated_cluster_score$consensus_call) &
                                               activating_annotated_cluster_score$consensus_call=="activating"] = "activating mutation"
activating_clusters_two_mut = activating_clusters[duplicated(activating_clusters)]
activating_annotated_cluster_score_duo = activating_annotated_cluster_score[activating_annotated_cluster_score$Cluster %in% activating_clusters_two_mut,]
# summarize: number of phosphosites and driver mutation and non-driver mutation in each clusters
counts = data.frame(table(activating_annotated_cluster_score$Gene_Drug,activating_annotated_cluster_score$Cluster,activating_annotated_cluster_score$site_type))
total_counts = counts
colnames(total_counts) = c("gene","cluster","site_type","count")
total_counts_clean = total_counts[total_counts$count>0,]
p = ggplot(total_counts_clean,aes(x=paste(gene), y=count,fill=site_type))
p = p + facet_grid(.~cluster,scale="free",space="free",drop=T) 
p = p + geom_bar(stat = "identity") + theme_bw() #+scale_size_area() + scale_size_continuous(range = c(0,max(total_counts$count)))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom") #+ ylim(0,5)# note that this will cut off some values
p = p + labs(y= "Count", x= "")
p
outFile = paste("output/coclustered_residues_to_functional_mutation_counts.pdf")
ggsave(file=outFile, w=6, h = 4, useDingbats=FALSE)

## phosphosites next to activating mutations
# plot and save outputs
p = ggplot(activating_annotated_cluster_score,aes(x=as.factor(Gene_Drug), y=Geodesic_From_Centroid,color=site_type))
p = p + facet_grid(.~Cluster,scale="free",space="free")
#p = p + facet_wrap(~Cluster,drop=T)
p = p + geom_point(alpha=0.3, height = 0, width = 0.2) + theme_bw() #+ theme_nogrid()
p = p + geom_text_repel(aes(label = ifelse(!duplicated(Mutation_Gene) & (site_type=="phosphosite" | consensus_call=="activating"), Mutation_Gene, NA)),alpha=0.7) + theme_bw() #+ theme_nogrid()
#p = p + geom_text_repel(aes(label = ifelse(!duplicated(Mutation_Gene) & (site_type=="phosphosite"), Mutation_Gene, NA)),alpha=0.7) + theme_bw() #+ theme_nogrid()
#p = p + geom_text(aes(label = ifelse(!duplicated(Mutation_Gene) & (site_type=="phosphosite"), Mutation_Gene, NA)),alpha=0.9, vjust=-0.5) + theme_bw() #+ theme_nogrid()
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom") #+ ylim(0,5)# note that this will cut off some values
p = p + labs(y= "Geodesic From Centroid", x= "Gene")
# p = p + geom_hline(yintercept = 1, alpha=0.3)
p
outFile = paste("output/coclustered_residues_to_functional_mutation.pdf")
ggsave(file=outFile, w=8, h = 6, useDingbats=FALSE)
