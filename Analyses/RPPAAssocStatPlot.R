##### pathVarRPPAAssoc.R #####
# Kuan-lin Huang @ WashU 2017 Dec
# statistically evaluate clustered mutations and their RPPA effect

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")


#FUNCTION myglm#
myglm=function(z,trait,variant,covar=NA,ytype) {
  if (is.na(covar) | is.null(covar) | nchar(covar)==0  ) { 
    model=formula(paste(trait,"~",variant)) 
  } else {
    model=formula(paste(trait,"~",variant,"+",covar))
  }
  if (ytype=="B") fit=glm(formula=model,data=z,family=binomial(link = "logit"))
  if (ytype=="Q") fit=glm(formula=model,data=z,family=gaussian(link = "identity"))
  fit
}
#END myglm

#FUNCTION plot_violin#
### plot quantitative y in relationship to x
### faceted by covariate covi
plot_violin = function(y,yi,xi,covi){ 
  # change height based on the number of markers
  dat = y[,c(yi, xi, covi)]
  dat[,1] = as.numeric(dat[,1])
  #dat.m = melt(mut_exp, id.var="Status")
  # plot violin plots faceted by marker genes
  n_facet = length(unique(dat[,covi]))
  p = ggplot(data=dat)
  p = p + facet_grid(as.formula(paste(". ~", covi)))
  p = p + geom_violin(aes_string(x=xi, y=yi, fill=xi),alpha=0.5) + guides(fill=FALSE, color =FALSE) 
  p = p + geom_jitter(aes_string(x=xi, y=yi, color=xi), alpha = 0.2) #+ geom_point(aes(x=Status, y=value)) 
  p = p + guides(fill=FALSE, color =FALSE) 
  p = p + labs(x = paste(xi,"pathVarPline Variant"), y = yi) + theme_bw()
  p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8, angle = 90))
  p
  fn = paste(pd, yi, "by", xi, "cross", covi, "violin.pdf", sep="_")
  ggsave(file=fn, width = n_facet, useDingbats=FALSE)
}

run_glm = function(data=NULL, covi="") {
  
  ### data input #####
  data_g = data#[, c("age_at_initial_pathologic_diagnosis","type",gene)]
  
  ######### analysis ##########
  row_stat=NULL
  
  # analysis_type   clinical_data_trait_name    variant/gene_name   covariates  memo
  ytype="Q";yi="expression";xi="Clustered";covi=covi;#covi=md[i,4];memo=md[i,5]
  # DEBUG
  #cat(paste("    Processing: yi =", yi, " xi =", xi, " covi =", covi, "\n") )
  
  if (ytype=="Q") {
    test = "F"
  } else if (ytype=="B") {
    test = "Chisq"
  } else {
    stop("Unknown model ytype ", ytype) 
  }
  
  glm = try(myglm(data_g,yi,xi,covi,ytype))  # MAW new
  if(class(glm)[1] == "try-error") {
    cat(paste("    Error caught, continuing.  yi =", yi, " xi =", xi, " covi =", covi, "\n") )
    next
  }
  try(anova(glm,test=test))->fit
  
#   if (xi %in% names(coefficients(glm))){
#     coeff = coefficients(glm)[[xi]]
  if (length(names(coefficients(glm)))>1){
    coeff = coefficients(glm)[[2]]
  } else {coeff = NA}
  
  if (class(fit)[1]!="try-error")
  {
    fit=as.matrix(fit)
    if (xi %in% rownames(fit)) (row_stat = cbind(yi,ytype,xi,as.data.frame(t(fit[xi,])),coeff,covi))
  }
  return(row_stat)
  
}

### MAIN ###

## read input files
## phopshosites
site_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201803/3D_Proximity.musites.gz"
site = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(site_f))
site_uniq = site[!duplicated(paste(site$Gene1,site$Mutation1,site$Transcript2,site$TranscriptPosition2)),]

phosphosites = site[,c("Gene2","Site2")]
phosphosites = phosphosites[!duplicated(paste(phosphosites$Gene2,phosphosites$Site2)),]
phosphosites$Position = as.numeric(gsub("p.[A-Z](.*)","\\1",phosphosites$Site2))

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
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Somatic_HGVSp) %in% paste(annotated_cluster$Gene_Drug,annotated_cluster$Mutation_Gene)] = "Clustered"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(phosphosites$Gene2,phosphosites$Position-1)] = "Proximal"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(phosphosites$Gene2,phosphosites$Position-2)] = "Proximal"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(phosphosites$Gene2,phosphosites$Position+1)] = "Proximal"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(phosphosites$Gene2,phosphosites$Position+2)] = "Proximal"
somatic$Type[paste(somatic$HUGO_Symbol,somatic$Position) %in% paste(phosphosites$Gene2,phosphosites$Position)] = "Direct"

## RPPA (processed data from germline project)
RPPA = read.table("/Users/khuang/Box Sync/PhD/germline/PanCanAtlasGermline/analysis/RPPA_effect/out/pancan_RPPA_quantile_all.tsv",header=T, stringsAsFactors = F, quote = "", sep = "\t")

RPPA_g = RPPA[RPPA$genes %in% annotated_cluster$Gene_Drug,]
RPPA_g_m = RPPA_g[,c("marker","genes","bcr_patient_barcode","expression","cancer")]
RPPA_g_m$marker = gsub(".*\\|","",RPPA_g_m$marker)

colnames(RPPA_g_m)=c("marker","Gene","bcr_patient_barcode","expression","cancer")
RPPA_g_m$expression = as.numeric(RPPA_g_m$expression)


##### individual cancer type analysis #####
cancers = unique(RPPA_g_m$cancer)
cancers = cancers[-which(cancers %in% c("COADREAD","GBMLGG","STES","KIPAN"))]
markers = unique(RPPA_g_m$marker)
# limit runs to cancers with at least 3 likely patho/pathogenic variants 
tt = NULL
for (marker in markers){
  RPPA_g_m_g = RPPA_g_m[RPPA_g_m$marker==marker,]
  gene = RPPA_g_m_g$Gene[1]
  gene_sample_g = somatic[somatic$HUGO_Symbol==gene & somatic$Type == "Clustered",]
  var_exp_g = merge(RPPA_g_m_g,gene_sample_g,by="bcr_patient_barcode",all.x=T)
  var_exp_g$Clustered = 0
  var_exp_g$Clustered[!is.na(var_exp_g$Type) & var_exp_g$Type=="Clustered"] = 1
  for (cancer in cancers){
    var_exp_g_c = var_exp_g[var_exp_g$cancer %in% cancer,]
    gene_path_count = sum(var_exp_g_c$Clustered)
    if (gene_path_count > 2){
      # run GLM
      w = wilcox.test(var_exp_g_c$expression[var_exp_g_c$Clustered==0],var_exp_g_c$expression[var_exp_g_c$Clustered!=0])
      wP = w$p.value
      wWstat = w$statistic
      cancer_gene_stat = run_glm(var_exp_g_c)
      # compile results
      full_cancer_gene_stat = cbind(cancer,marker,gene,gene_path_count,wP,wWstat,cancer_gene_stat)
      tt = rbind(tt, full_cancer_gene_stat)
    }
  }
}

#"yi","ytype","xi","Df","Deviance","Resid. Df","Resid. Dev","F","Pr(>F)","covi","memo"
colnames(tt) = c("cancer","marker","gene","gene_path_count","wilcoxP","W_stat","y","y_type","Gene","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F_statistic","p-value","coefficient","covariants");
tt$FDR = p.adjust(tt[,"p-value"], method="fdr") # MAW new, calculates FDR based on the method from,
# Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57, 289–300.
tt$wilcoxFDR = p.adjust(tt[,"wilcoxP"], method="fdr")

tt=tt[order(tt$FDR, decreasing=FALSE),]
tn = "output/clusterMutRPPAAssoc.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

### plotting ###
tt$gene = as.character(tt$gene)
tt$marker = as.character(tt$marker)
tt$FDR_plot = tt$FDR
tt$FDR_plot[tt$FDR_plot<10^(-8)]= 0.95*10^(-8)
p = ggplot(data=tt)
p = p + geom_point(aes(y=-log10(FDR_plot),x= coefficient,color = cancer),alpha=0.5)
#p = p + geom_text_repel(aes(y=-log10(FDR),x= cohort_AF,label=ifelse(FDR<0.05, Gene,NA),color = Cancer))#,alpha=1.3)
#p = p + geom_point(aes(y=cohort_AF,x=Cancer,size=-log10(FDR),color = Cancer))
p = p + geom_text_repel(aes(y=-log10(FDR_plot),x= coefficient,color = cancer,label=ifelse(FDR<0.05, marker,NA)))
p = p + getPCACancerColor() + xlim(-1.9,1.9)
p = p + labs(x="Coefficient",y= "-log10(FDR)")
p = p + geom_vline(xintercept = 0, alpha=0.5)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p = p + theme(legend.position = "bottom")
p
fn = 'output/RPPAAssocVolcanoGLM.pdf'
ggsave(fn,w = 4, h = 6, useDingbat=F)

tt$wilcoxFDR_plot = tt$wilcoxFDR
tt$wilcoxFDR_plot[tt$wilcoxFDR_plot<10^(-8)]= 0.95*10^(-8)
p = ggplot(data=tt)
p = p + geom_point(aes(y=-log10(wilcoxFDR_plot),x= coefficient,color = cancer),alpha=0.5)
#p = p + geom_text_repel(aes(y=-log10(FDR),x= cohort_AF,label=ifelse(FDR<0.05, Gene,NA),color = Cancer))#,alpha=1.3)
#p = p + geom_point(aes(y=cohort_AF,x=Cancer,size=-log10(FDR),color = Cancer))
p = p + geom_text_repel(aes(y=-log10(wilcoxFDR_plot),x= coefficient,color = cancer,label=ifelse(FDR<0.05, marker,NA)))
p = p + getPCACancerColor()
p = p + labs(x="Coefficient",y= "-log10(FDR)")
p = p + geom_vline(xintercept = 0, alpha=0.5) + xlim(-1.6,1.6)
p = p  + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'output/RPPAAssocVolcanoWCOX.pdf'
ggsave(fn,w = 6, h = 5, useDingbat=F)

colnames(RPPA_g)[9] = "HUGO_Symbol"
RPPA_g_mut = merge(RPPA_g,somatic,by=c("HUGO_Symbol","bcr_patient_barcode"))
featMarkers = tt$marker[tt$FDR< 0.05]
featMarkers = featMarkers[which(featMarkers!="B-Raf" & featMarkers!="VHL")]
RPPA_g_mut_fg = RPPA_g_mut[gsub(".*\\|","",RPPA_g_mut$marker) %in% featMarkers,]
RPPA_g_mut_fg$Type = factor(RPPA_g_mut_fg$Type,levels = c("None","Clustered","Proximal","Direct"))
RPPA_g_mut_fg$marker2 = gsub(".*\\|","",RPPA_g_mut_fg$marker)
#RPPA_g_mut_fg$quantile = as.numeric(RPPA_g_mut_fg$quantile)

for ( marker in tt$marker[tt$FDR<0.05]){
  cancers = tt$cancer[tt$marker==marker & tt$FDR<0.05]
  RPPA_g_mut_fg_c = RPPA_g_mut_fg[RPPA_g_mut_fg$cancer %in% cancers & RPPA_g_mut_fg$marker2==marker,]
  p = ggplot(RPPA_g_mut_fg_c,aes(x=Type,y=quantile, fill=Type, color = Type))
  p = p + facet_grid(.~cancer, scale = "free", space = "free", drop=T)
  p = p + geom_dotplot(dotsize=1.2,binwidth=.01, binaxis= "y",colour=NA,stackdir ="centerwhole",alpha=0.5)
  p = p + geom_violin(alpha=0.6)
  #p = p + geom_jitter(alpha = 0.2,height = 0)
  p = p + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.8)
  p = p + geom_text_repel(aes(color=Type, label=ifelse(Type=="Clustered" & quantile > 0.95 & HUGO_Symbol != "TP53",as.character(Somatic_HGVSp),NA)), alpha = 0.7)
  p = p + theme_bw() 
  p = p + getTypeFillScale() + getTypeColorScale()
  p = p + ylab("RPPA Expression Percentile within the Cancer Type") + xlab("Interaction of Mutation with Phosphosites")
  p = p + scale_y_continuous(breaks = seq(0,1, by= 0.25))
  #p = p + getVarColorScale()
  p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))
  p
  fn = paste("output/clusteredMutRPPAExpression_byGene",marker,"pdf",sep=".")
  ggsave(file=fn, useDingbats=FALSE)
}

p = ggplot(RPPA_g_mut_fg,aes(x=Type,y=quantile, fill=Type))
p = p + facet_grid(.~marker, scale = "free", space = "free", drop=T)
#p = p + geom_dotplot(dotsize=1.2,binwidth=.01, binaxis= "y",colour=NA,stackdir ="centerwhole")
p = p + geom_violin(alpha=0.2)
#p = p + geom_point(aes(color=ifelse(Type!="None" & quantile > 0.97 & HUGO_Symbol != "TP53",Type,NA)), alpha = 0.7)
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.8)
p = p + geom_text_repel(aes(color=Type, label=ifelse(Type=="Clustered" & quantile > 0.97 & HUGO_Symbol != "TP53",as.character(Somatic_HGVSp),NA)), alpha = 0.7)
p = p + theme_bw() 
p = p + getTypeFillScale() + getTypeColorScale()
p = p + ylab("RPPA Expression Percentile within the Cancer Type") + xlab("Interaction of Mutation with Phosphosites")
p = p + scale_y_continuous(breaks = seq(0,1, by= 0.25))
#p = p + getVarColorScale()
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))
p
fn = "output/clusteredMutRPPAExpression_byGene.pdf"
ggsave(file=fn, w=15, h=6,useDingbats=FALSE)

RPPA_g_mut_fg_smut = RPPA_g_mut_fg[RPPA_g_mut_fg$Type!="None" & RPPA_g_mut_fg$quantile> 0.95,]
RPPA_g_mut_fg_smut = RPPA_g_mut_fg_smut[order(RPPA_g_mut_fg_smut$quantile, decreasing = T),]
tn = "output/clusterMutHighRPPA.txt"
write.table(RPPA_g_mut_fg_smut, quote=F, sep="\t", file = tn, row.names = F)
