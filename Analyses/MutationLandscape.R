##### pathVarRPPAAssoc.R #####
# Kuan-lin Huang @ WashU 2017 Dec
# statistically evaluate clustered mutations and their RPPA effect

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")

##### MAIN #####

### read input files ###
# ## phopshosites
# site_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201805/3D_Proximity_cleaned.musites.gz"
# site = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(site_f))
# site_uniq = site[!duplicated(paste(site$Gene1,site$Mutation1,site$Transcript2,site$TranscriptPosition2)),]
#
# phosphosites = site[,c("Transcript","Site2")]
# phosphosites = phosphosites[!duplicated(paste(phosphosites$Transcript,phosphosites$Site2)),]
# phosphosites$Position = as.numeric(gsub("p.[A-Z](.*)","\\1",phosphosites$Site2))
annotated_clusterPTM_h = annotated_cluster_h[annotated_cluster_h$Alternate=="ptm",]
annotated_clusterPTM_h$Position = as.numeric(annotated_clusterPTM_h$Position)

## somatic mutations
somatic_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/somatic/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz"
somatic = read.table(header=T, quote = "", sep="\t", file = gzfile(somatic_f), stringsAsFactors=FALSE)
colnames(somatic) = c("HUGO_Symbol","Somatic_Variant_Classification","TumorSample","Somatic_HGVSp") #Hugo_Symbol     Variant_Classification  Tumor_Sample_Barcode    HGVSp_Short
somatic$bcr_patient_barcode = gsub("(^TCGA-[A-Z0-9][A-Z0-9]-[A-Z0-9][A-Z0-9][A-Z0-9][A-Z0-9])-.*","\\1",somatic$TumorSample)
somatic = somatic[somatic$Somatic_Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                                                                "Missense_Mutation","Nonsense_Mutation","Splice_Site"),]
somatic = somatic[somatic$HUGO_Symbol %in% annotated_cluster$Gene_Drug,]

somatic$Position = as.numeric(gsub("p.[A-Z](.*)[A-Z]","\\1",somatic$Somatic_HGVSp))
somatic$Type = "None"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Somatic_HGVSp) %in% paste(annotated_cluster_h$Gene_Drug,annotated_cluster_h$Mutation_Gene)] = "Clustered"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(annotated_clusterPTM_h$Gene_Drug,annotated_clusterPTM_h$Position-1)] = "Proximal"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(annotated_clusterPTM_h$Gene_Drug,annotated_clusterPTM_h$Position-2)] = "Proximal"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(annotated_clusterPTM_h$Gene_Drug,annotated_clusterPTM_h$Position+1)] = "Proximal"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(annotated_clusterPTM_h$Gene_Drug,annotated_clusterPTM_h$Position+2)] = "Proximal"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(annotated_clusterPTM_h$Gene_Drug,annotated_clusterPTM_h$Position)] = "Direct"

clin_f = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/clinical/PanCan_ClinicalData_V4_wAIM.txt"
clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)
clin_som = clin[clin$bcr_patient_barcode %in% somatic$bcr_patient_barcode,]
clin_tab = clin_som[,c('type',"bcr_patient_barcode")]
colnames(clin_tab)[1] = "cancer"
sample_size = data.frame(table(clin_som$type))
colnames(sample_size) = c("Cancer","Cohort_size")

somatic_clin = merge(somatic,clin_tab,by="bcr_patient_barcode")
somatic_clin_uniq = somatic_clin[!duplicated(paste(somatic_clin$HUGO_Symbol,somatic_clin$bcr_patient_barcode,somatic_clin$Type)),]
somatic_clin_count = data.frame(table(somatic_clin_uniq$HUGO_Symbol,somatic_clin_uniq$Type,somatic_clin_uniq$cancer))
colnames(somatic_clin_count) = c("Gene","Type","Cancer","Count")

somatic_clin_count_freq = merge(somatic_clin_count,sample_size,by="Cancer")
somatic_clin_count_freq$Freq = round(somatic_clin_count_freq$Count/somatic_clin_count_freq$Cohort_size,digits = 3)
sele_g = names(rank_vectors(somatic_clin_uniq$HUGO_Symbol[somatic_clin_uniq$Type!="None"]))


somatic_clin_count_freq_g = somatic_clin_count_freq[somatic_clin_count_freq$Gene %in% sele_g & somatic_clin_count_freq$Type != "None",]

getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
somatic_clin_count_freq_g$Count_plot = somatic_clin_count_freq_g$Count
somatic_clin_count_freq_g$Count_plot[somatic_clin_count_freq_g$Count_plot>30]= 30
somatic_clin_count_freq_g$Type = factor(somatic_clin_count_freq_g$Type, levels = c("Direct","Proximal","Clustered"))
# combined_sum_f_added_g$Cancer = factor(combined_sum_f_added_g$Cancer, levels = unique(cancer_count$Cancer[cancer_count$Cancer!="All"]))

p = ggplot(data=somatic_clin_count_freq_g)
p = p + facet_grid(Type~.,drop=T,scales = "free_y", space = "free_y")
p = p + geom_tile(aes(x=Cancer, y=Gene, fill=Count_plot), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
p = p + geom_text(aes(x=Cancer, y=Gene, label = ifelse(Count>0,Count,NA)), color="black", size=3)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'output/mutation_count_direct2clustered_heatmap.pdf'
ggsave(file=fn, height=10, width=9, useDingbats=FALSE)


somatic_clin_count_freq_g$Freq[somatic_clin_count_freq_g$Freq>0.1] = round(somatic_clin_count_freq_g$Freq[somatic_clin_count_freq_g$Freq>0.1],digits = 2)

somatic_clin_count_freq_g$Freq_plot = somatic_clin_count_freq_g$Freq
somatic_clin_count_freq_g$Freq_plot[somatic_clin_count_freq_g$Freq_plot>1]= 1
somatic_clin_count_freq_g$Type = factor(somatic_clin_count_freq_g$Type, levels = c("Direct","Proximal","Clustered"))
# combined_sum_f_added_g$Cancer = factor(combined_sum_f_added_g$Cancer, levels = unique(cancer_count$Cancer[cancer_count$Cancer!="All"]))

p = ggplot(data=somatic_clin_count_freq_g)
p = p + facet_grid(Type~.,drop=T,scales = "free_y", space = "free_y")
p = p + geom_tile(aes(x=Cancer, y=Gene, fill=Freq_plot), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
p = p + geom_text(aes(x=Cancer, y=Gene, label = ifelse(Freq>0.01,Freq*100,NA)), color="black", size=3)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'output/mutation_freq_direct2clustered_heatmap.pdf'
ggsave(file=fn, height=10, width=9, useDingbats=FALSE)
