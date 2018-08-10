##### clustered_CPTAC_phosphosites.R #####
# Kuan-lin Huang @ WashU 201712

####### positional effect: mutated aa vs. phospho aa ###
##### dependencies #####
### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")

baseD = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/"

##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_SOMATIC_formatted_amino_acid.txt",sep=""))
colnames(BRCA_mut)[1]= "Gene"
BRCA_mut_m = melt(BRCA_mut,id.vars = "Gene")

### Phosphosites ###
BRCA_Pho = read.table(stringsAsFactors = F,header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/input/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv")
colnames(BRCA_Pho) = gsub(".Log.Ratio","",colnames(BRCA_Pho))
clustered_sites = read.table(stringsAsFactors = F,header=TRUE, sep="\t", file="output/clusteredCPTACsites.txt")
BRCA_Pho$clustered = F
BRCA_Pho$clustered[BRCA_Pho$Phosphosite %in% clustered_sites$Phosphosite] = T
BRCA_Pho = BRCA_Pho[rowSums(!is.na(BRCA_Pho))>9,]
BRCA_Pho$SD = rowSds(data.matrix(BRCA_Pho[,c(2:109)]),na.rm = T)
BRCA_Pho_g = BRCA_Pho[BRCA_Pho$Gene %in% clustered_sites$hgnc_symbol,]

p = ggplot(BRCA_Pho_g,aes(y = SD, x = clustered, fill=clustered))
p = p + geom_violin(alpha=0.5)  #+ ylim(0,20)#+ scale_y_log10()
p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.8)
p = p + theme_bw()
p = p + labs(y = "SD for the phosphosites")
p = p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=14),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = paste("output/sd_for_clustered_phosphosites.pdf",sep="_")
ggsave(fn,useDingbat=F)

BRCA_Pho_sites = BRCA_Pho[BRCA_Pho$Phosphosite %in% clustered_sites$Phosphosite,]
BRCA_Pho_sites = BRCA_Pho_sites[rowSums(!is.na(BRCA_Pho_sites))>9,]
BRCA_Pho_sites_m = melt(BRCA_Pho_sites, id.vars = c("Phosphosite","Peptide","Gene","Organism"))

BRCA_Pho_sites_mut = merge(BRCA_Pho_sites_m,BRCA_mut_m,by=c("Gene","variable"))
colnames(BRCA_Pho_sites_mut) = c("Gene","Sample","Phosphosite","Peptide","Organism","Expression","Mutation") 
cat("incidence of mutations and phosphosites in the same sample","\n")
table(BRCA_Pho_sites_mut$Gene[BRCA_Pho_sites_mut$Mutation!="wt"])

##### plot collective violin plot #####
### make violin plots for selected mutated gene-marker combo###
p = ggplot(data=BRCA_Pho_sites_mut,aes(x=Phosphosite, y=Expression, color=Mutation, fill= Gene))
p = p + facet_grid(.~Gene, scale="free", space="free")
#p = p + geom_jitter(alpha=0.3,height=0) #+ guides(fill=FALSE) 
p = p + geom_violin(alpha=0.3)
p = p + geom_text(aes(label = ifelse(Mutation != "wt",Mutation, NA)), size=3,alpha=0.8)
p = p + theme_bw()
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14,angle=90,vjust=0.5), 
              axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
p = p + theme(legend.position = "none")
p
fn = "output/clustered_phosphosites_expression.pdf"
ggsave(file=fn, w=16,useDingbats=FALSE)

### clustering with subtype
brca_clin = read.table(row.names=1, header=TRUE, sep="\t", quote = "", 
                       file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/BRCA_clinical_summary.txt")
basal = as.vector(colnames(brca_clin)[brca_clin["pam50",]=="Basal"])
LumA  = as.vector(colnames(brca_clin)[brca_clin["pam50",]=="LumA"])
LumB  = as.vector(colnames(brca_clin)[brca_clin["pam50",]=="LumB"])
Her2 = as.vector(colnames(brca_clin)[brca_clin["pam50",]=="Her2"])


BRCA_Pho_sites_clust = BRCA_Pho_sites[,-c(which(colnames(BRCA_Pho_sites) %in% c("Phosphosite","Peptide","Gene","Organism")))]
row.names(BRCA_Pho_sites_clust) = paste(BRCA_Pho_sites$Gene,BRCA_Pho_sites$Phosphosite)

BRCA_Pho_sites_clust_lowNA = BRCA_Pho_sites_clust[rowSums(!is.na(BRCA_Pho_sites_clust))>50,colnames(BRCA_Pho_sites_clust) %in%colnames(brca_clin)]
BRCA_Pho_sites_clust_lowNA = as.matrix(BRCA_Pho_sites_clust_lowNA)
samples=colnames(BRCA_Pho_sites_clust_lowNA)
samples[samples %in% basal]="red"
samples[samples %in% LumA]="darkblue"
samples[samples %in% LumB]="blue"
samples[samples %in% Her2]="pink"

pdf("output/clustered_phosphosites_heatmap.pdf",width=15)
par(oma=c(3,5,3,5))
hm = heatmap.2(BRCA_Pho_sites_clust_lowNA, trace="none", na.color="white", notecol="black",
                        cexRow=0.8,cexCol=0.8, scale="none",dendrogram='column', ColSideColors = samples,
               col=getPalette, margins=c(5,5))
                       # labRow=NA,labCol=NA,

par(lend = 1)  
# legend("bottomleft",    # location of the legend on the heatmap plot
#        # legend = c("Basal BRCA", "BRCA", "OV", "CRC"), # category labels
#        # col = c("green","forestgreen", "orange","purple"),  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )

dev.off()


### do they affect transcription? ###
histones = read.table("input/histones.tsv",header=F,sep="\t")
histone_genes = histones[,1]
### expression ###
BRCA_rna = read.table(header=TRUE, sep="\t", file=paste(baseD,"pan3can_shared_data/BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))
colnames(BRCA_rna)[1]= "Gene"
row.names(BRCA_rna) = BRCA_rna$Gene
BRCA_rna = BRCA_rna[,-which(colnames(BRCA_rna)=="Gene")]
BRCA_rna_PCA = prcomp(t(BRCA_rna))
BRCA_rna_PCA_data = data.frame(BRCA_rna_PCA$x)


BRCA_Pho_sites_interest = BRCA_Pho_sites[,-c(which(colnames(BRCA_Pho_sites) %in% c("Peptide","Organism","clustered","SD")))]
#histone_site = BRCA_Pho_sites_interest[BRCA_Pho_sites_interest$Gene %in% histone_genes,]
BRCA_Pho_sites_interest_HI = t(BRCA_Pho_sites_interest[BRCA_Pho_sites_interest$Phosphosite=="NP_778224.1:s48",])
BRCA_Pho_sites_interest_HI = t(BRCA_Pho_sites_interest[BRCA_Pho_sites_interest$Phosphosite=="NP_778224.1:y89",])
colnames(BRCA_Pho_sites_interest_HI) = "phosphosite"

BRCA_rna_pho = merge(BRCA_rna_PCA_data,BRCA_Pho_sites_interest_HI,by="row.names")
row.names(BRCA_rna_pho) = BRCA_rna_pho[,1]
BRCA_rna_pho = BRCA_rna_pho[,-which(colnames(BRCA_rna_pho)=="Row.names")]
pam50 = t(brca_clin[1,])
BRCA_rna_pho_pam50 = merge(BRCA_rna_pho,pam50,by="row.names")

p = ggplot(data=BRCA_rna_pho_pam50,aes(x=PC1, y=PC2, color=as.numeric(phosphosite)))
p = p + facet_grid(pam50~.)
p = p + geom_point() #+ guides(fill=FALSE)
p = p + theme_bw() + scale_color_gradientn(colors=RdBu1024,n=1000)
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14,angle=90,vjust=0.5),
              axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
p = p + labs(title = BRCA_Pho_sites_interest_HI[1,1])
p
fn = paste("output/expression_PCs_vs_histone_phosphosites_HIST4H4p.y89.pdf",sep=".")
ggsave(file=fn, w=6, h=5, useDingbats=FALSE)

# # no obvious findings
# 
# # germline
# fn = "/Users/khuang/Box\ Sync/PhD/germline/PanCanAtlasGermline/TCGA_data/germline/PCA_pathVar_integrated_filtered_adjusted_wSomaticCounts.tsv"
# pathVar = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=fn)
# pathVarP = pathVar[pathVar$Overall_Classification %in% c("Pathogenic","Likely Pathogenic"),]
# pathVarPbrca = pathVarP$bcr_patient_barcode[pathVarP$HUGO_Symbol %in% c("BRCA1","BRCA2") & pathVarP$LOH_FDR < 0.05]
# BRCA_rna_pho_pam50$brca = "WT"
# BRCA_rna_pho_pam50$brca[paste("TCGA.",gsub(".01A","",BRCA_rna_pho_pam50$Row.names),sep="") %in% gsub("-",".",pathVarPbrca)] = "BRCA"
# 
# p = ggplot(data=BRCA_rna_pho_pam50,aes(x = brca, y=as.numeric(phosphosite)))
# p = p + facet_grid(.~pam50)
# p = p + geom_point() #+ guides(fill=FALSE)
# p = p + theme_bw() + scale_color_gradientn(colors=RdBu1024,n=1000)
# p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14,angle=90,vjust=0.5),
#               axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
# #p = p + theme(legend.position = "none")
# p
# # no obvious findings

# BRCA_Pho_sites_interest_HI = t(BRCA_Pho_sites_interest[BRCA_Pho_sites_interest$Gene=="HIST4H4",])
# colnames(BRCA_Pho_sites_interest_HI) = c("S48","Y89")
# data = data.frame(BRCA_Pho_sites_interest_HI[-c(1,110),])
# data$S48 = as.numeric(data$S48)
# data$Y89 = as.numeric(data$Y89)
# p = ggplot(data,aes(x = S48, y=Y89))
# #p = p + facet_grid(.~pam50)
# p = p + geom_point() #+ guides(fill=FALSE)
# #p = p + theme_bw() + scale_color_gradientn(colors=RdBu1024,n=1000)
# p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14,angle=90,vjust=0.5),
#               axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
# #p = p + theme(legend.position = "none")
# p
# # the two phosphosites on HIST4H4 are not well-correlated!
