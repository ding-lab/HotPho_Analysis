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
validation_fasmic = validation[!is.na(validation$baf3_call) | !is.na(validation$mcflo_call),]
table(validation_fasmic$baf3_call,validation_fasmic$mcflo_call)
cat("number of BAF3 calls",sum(!is.na(validation$baf3_call)),"\n")
cat("number of MCF calls",sum(!is.na(validation$mcflo_call)),"\n")

# functional mutation in clusters
annotated_cluster_score = merge(annotated_cluster,validation_fasmic,by=c("Transcript","Mutation_Gene"),all.y=T)
geneList = unique(annotated_cluster_score$Gene_Drug)

annotated_cluster_score_g = annotated_cluster_score[annotated_cluster_score$Gene %in% geneList,]
dim(annotated_cluster_score_g)
annotated_cluster_score_g$coclustering = "None"
annotated_cluster_score_g$coclustering[!is.na(annotated_cluster_score_g$Gene_Drug)] = "Co-clustered with Phosphosites"
table(annotated_cluster_score_g$coclustering)

cat("Distribution of mutational calls in each of the dataset\n")
table(annotated_cluster_score_g$baf3_call)
table(annotated_cluster_score_g$mcflo_call)

# output percentage table
annotated_cluster_score_g_count_baf3 = data.frame(table(annotated_cluster_score_g$Gene,annotated_cluster_score_g$coclustering,annotated_cluster_score_g$baf3_call == "activating"))
colnames(annotated_cluster_score_g_count_baf3) = c("gene","coclustering","activating","count")
annotated_cluster_score_g_count_baf3_act = annotated_cluster_score_g_count_baf3[annotated_cluster_score_g_count_baf3$activating==TRUE,c("gene","coclustering","count")]
annotated_cluster_score_g_count_baf3_nonact = annotated_cluster_score_g_count_baf3[annotated_cluster_score_g_count_baf3$activating!=TRUE,c("gene","coclustering","count")]
annotated_cluster_score_g_count_baf3_combined = merge(annotated_cluster_score_g_count_baf3_act,annotated_cluster_score_g_count_baf3_nonact, by = c("gene","coclustering"))
colnames(annotated_cluster_score_g_count_baf3_combined)[3:4] = c("activating_count","nonactivating_count")
annotated_cluster_score_g_count_baf3_combined$activating_percentage = annotated_cluster_score_g_count_baf3_combined$activating_count/(annotated_cluster_score_g_count_baf3_combined$activating_count +annotated_cluster_score_g_count_baf3_combined$nonactivating_count)
annotated_cluster_score_g_count_baf3_combined$cell_line = "BAF3"
tn = "output/baf3_cocluster_vs_activating_percentage.txt"
write.table(annotated_cluster_score_g_count_baf3_combined[order(annotated_cluster_score_g_count_baf3_combined$coclustering),], quote=F, sep="\t", file = tn, row.names = F)

annotated_cluster_score_g_count_mcflo = data.frame(table(annotated_cluster_score_g$Gene,annotated_cluster_score_g$coclustering,annotated_cluster_score_g$mcflo_call == "activating"))
colnames(annotated_cluster_score_g_count_mcflo) = c("gene","coclustering","activating","count")
annotated_cluster_score_g_count_mcflo_act = annotated_cluster_score_g_count_mcflo[annotated_cluster_score_g_count_mcflo$activating==TRUE,c("gene","coclustering","count")]
annotated_cluster_score_g_count_mcflo_nonact = annotated_cluster_score_g_count_mcflo[annotated_cluster_score_g_count_mcflo$activating!=TRUE,c("gene","coclustering","count")]
annotated_cluster_score_g_count_mcflo_combined = merge(annotated_cluster_score_g_count_mcflo_act,annotated_cluster_score_g_count_mcflo_nonact, by = c("gene","coclustering"))
colnames(annotated_cluster_score_g_count_mcflo_combined)[3:4] = c("activating_count","nonactivating_count")
annotated_cluster_score_g_count_mcflo_combined$activating_percentage = annotated_cluster_score_g_count_mcflo_combined$activating_count/(annotated_cluster_score_g_count_mcflo_combined$activating_count +annotated_cluster_score_g_count_mcflo_combined$nonactivating_count)
annotated_cluster_score_g_count_mcflo_combined$cell_line = "MCF10A"
tn = "output/mcflo_cocluster_vs_activating_percentage.txt"
write.table(annotated_cluster_score_g_count_mcflo_combined[order(annotated_cluster_score_g_count_mcflo_combined$coclustering),], quote=F, sep="\t", file = tn, row.names = F)

##### statistical test #####
out_table=character(0);
cell = "Ba/F3"
overall_co = colSums(annotated_cluster_score_g_count_baf3_combined[annotated_cluster_score_g_count_baf3_combined$coclustering!="None",c("activating_count","nonactivating_count")])
overall_no = colSums(annotated_cluster_score_g_count_baf3_combined[annotated_cluster_score_g_count_baf3_combined$coclustering=="None",c("activating_count","nonactivating_count")])
cat("Test for enrichment of activating mutations among co-clustering mutations in ",cell,"\n")
rbind(overall_co,overall_no)
fisher.test(rbind(overall_co,overall_no),alternative = "greater")

for(gene in unique(annotated_cluster_score_g_count_baf3_combined$gene)){
  test_gene = annotated_cluster_score_g_count_baf3_combined[annotated_cluster_score_g_count_baf3_combined$gene==gene,c("activating_count","nonactivating_count")]
  f.test = fisher.test(test_gene,alternative = "greater")
  OR = f.test$estimate
  p = f.test$p.value
  
  out_row = c(cell, gene, OR, p)
  out_table = rbind(out_table,out_row)
}

cell = "MCF10A"
overall_co = colSums(annotated_cluster_score_g_count_mcflo_combined[annotated_cluster_score_g_count_mcflo_combined$coclustering!="None",c("activating_count","nonactivating_count")])
overall_no = colSums(annotated_cluster_score_g_count_mcflo_combined[annotated_cluster_score_g_count_mcflo_combined$coclustering=="None",c("activating_count","nonactivating_count")])
cat("Test for enrichment of activating mutations among co-clustering mutations in ",cell,"\n")
rbind(overall_co,overall_no)
fisher.test(rbind(overall_co,overall_no),alternative = "greater")

for(gene in unique(annotated_cluster_score_g_count_mcflo_combined$gene)){
  test_gene = annotated_cluster_score_g_count_mcflo_combined[annotated_cluster_score_g_count_mcflo_combined$gene==gene,c("activating_count","nonactivating_count")]
  f.test = fisher.test(test_gene,alternative = "greater")
  OR = f.test$estimate
  p = f.test$p.value
  
  out_row = c(cell, gene, OR, p)
  out_table = rbind(out_table,out_row)
}
tn = "output/cocluster_vs_activating_percentage_tests.txt"
write.table(out_table, quote=F, sep="\t", file = tn, row.names = F,col.names = F)

##### plotting #####
annotated_cluster_score_g_plot = annotated_cluster_score_g
annotated_cluster_score_g_plot$baf3_call[!is.na(annotated_cluster_score_g_plot$baf3_call) & (annotated_cluster_score_g_plot$baf3_call %in% c("inhibitory","non-inhibitory","noninfo"))] = "other"
annotated_cluster_score_g_plot$mcflo_call[!is.na(annotated_cluster_score_g_plot$mcflo_call) & (annotated_cluster_score_g_plot$mcflo_call %in% c("inhibitory","non-inhibitory","noninfo"))] = "other"

p = ggplot(annotated_cluster_score_g_plot[!is.na(annotated_cluster_score_g_plot$baf3_call),],aes(x=Gene, fill = baf3_call))#,fill=site_type))
p = p + facet_grid(coclustering~.,scale="free",drop=T) 
p = p + geom_bar() + theme_nogrid() #+scale_size_area() + scale_size_continuous(range = c(0,max(total_counts$count)))
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom") #+ ylim(0,5)# note that this will cut off some values
p = p + labs(y= "Count", x= "Gene")
p
outFile = paste("output/fasmic_baf3_call_coclustered_vs_not.pdf")
ggsave(file=outFile, w=5,h=5,useDingbats=FALSE)

p = ggplot(annotated_cluster_score_g_plot[!is.na(annotated_cluster_score_g_plot$mcflo_call),],aes(x=Gene, fill = mcflo_call))#,fill=site_type))
p = p + facet_grid(coclustering~.,scale="free",drop=T) 
p = p + geom_bar() + theme_nogrid() 
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90,vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p = p + theme(legend.position="bottom") 
p = p + labs(y= "Count", x= "Gene")
p
outFile = paste("output/fasmic_mcflo_call_coclustered_vs_not.pdf")
ggsave(file=outFile, w=5,h=5, useDingbats=FALSE)


# cat("Mutation distribution of activating-mutation clusters:","\n")
# table(activating_annotated_cluster_score$Effect)

