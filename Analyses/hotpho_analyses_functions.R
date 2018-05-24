### HotPho analyses dependency and functions ###

### common dependencies for plotting ###

# dependencies

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library("gplots")
library(gridExtra)
library(matrixStats)
library(grid)
library(ggrepel)

# output folder
system("mkdir output")

### aesthetics  ###
col_paletteB = colorRampPalette(brewer.pal(9,"Blues"))
col_paletteR = colorRampPalette(brewer.pal(9,"Reds"))
RdBu = brewer.pal(9, "RdBu") 
getPalette = colorRampPalette(rev(RdBu))
RdBu1024 = colorRampPalette(rev(RdBu))(1024)

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank())

theme0 = function(...) theme( legend.position = "none",
                              panel.background = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              #panel.margin = unit(0,"null"),
                              axis.ticks = element_blank(),
                              axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.ticks.length = unit(0,"null"),
                              axis.ticks.margin = unit(0,"null"),
                              panel.border=element_rect(color=NA),...)

### functions ###
plot_barplot = function(matrix, x_string, fill_string=NULL, fileName="data.pdf", trans=NULL){
  fn = paste("output", fileName,sep ="/")
  
  if (is.null(fill_string)){
    p = ggplot(matrix,aes_string(x = x_string))
  } else{
    p = ggplot(matrix,aes_string(x = x_string, fill = fill_string))
  }
  
  p = p + geom_bar() + theme_bw()
  p = p + labs(x = x_string, y="Counts")
  p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
  if(!is.null(trans)){p = p + coord_trans(y=trans)}
  p
  ggsave(file=fn, w=5,h=5,useDingbats=FALSE)
}

rank_vectors = function(vector,n=20){
  table(vector)[order(table(vector),decreasing = T)][1:n]
}

plot_heatmap = function(data){
  getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
  p = ggplot(data=data)
  p = p + facet_grid(.~classification,drop=T,scales = "free", space = "free")
  p = p + geom_tile(aes(x=cancer, y=gene, fill=count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
  p = p + geom_text(aes(x=cancer, y=gene, label = count, stringsAsFactors=FALSE), color="black", size=3)
  p = p  + theme_bw() + theme_nogrid() +
    theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=14),axis.ticks = element_blank())#element_text(colour="black", size=14))
  return(p)
}

annotate_cluster = function(cluster){
  cluster$Total_count = NA
  cluster$Mut_count = NA
  cluster$Site_count = NA
  cluster$Type = NA
  
  for (cluster_num in unique(cluster$Cluster)){
    # count mutation and sites in each cluster
    # # previous version
    # cluster_sites = cluster$Mutation_Gene[cluster$Cluster == cluster_num]
    # alt_res = gsub("p.[A-Z][0-9]+","",cluster_sites)
    # num_total = length(alt_res)
    # num_site = sum(alt_res=="")
    # num_mut = num_total - num_site
    # new version simply using Alt
    cluster_alt = cluster$Alternate[cluster$Cluster == cluster_num]
    num_total = length(cluster_alt)
    num_site = sum(cluster_alt=="ptm")
    num_mut = num_total - num_site
    
    cluster$Total_count[cluster$Cluster == cluster_num] = num_total
    cluster$Mut_count[cluster$Cluster == cluster_num] = num_mut
    cluster$Site_count[cluster$Cluster == cluster_num] = num_site
    
    # classify cluster type
    if (num_mut > 0){
      if (num_site > 0){
        cluster$Type[cluster$Cluster == cluster_num] = "Hybrid"
      } else if (num_site == 0){
        cluster$Type[cluster$Cluster == cluster_num] = "Mut_Only"
      }
    } else if (num_site > 0){
      cluster$Type[cluster$Cluster == cluster_num] = "Site_Only"
    }
  }
  return(cluster)
}

getPCACancerColor = function() {
  # according to google spreadsheet: https://docs.google.com/spreadsheets/d/1Nb9mMkonAhZR1_2OI9nv4ylCei0LZjAf-2vYTRQcXKw/edit#gid=1704872109
  colors = c(
    "#C1A72F",
    "#FAD2D9",
    "#ED2891",
    "#F6B667",
    "#104A7F",
    "#9EDDF9",
    "#3953A4",
    "#007EB5",
    "#B2509E",
    "#97D1A9",
    "#ED1C24",
    "#F8AFB3",
    "#EA7075",
    "#754C29",
    "#D49DC7",
    "#CACCDB",
    "#D3C3E0",
    "#A084BD",
    "#542C88",
    "#D97D25",
    "#6E7BA2",
    "#E8C51D",
    "#7E1918",
    "#DAF1FC",
    "#00A99D",
    "#BBD642",
    "#00AEEF",
    "#BE1E2D",
    "#F9ED32",
    "#CEAC8F",
    "#FBE3C7",
    "#F89420",
    "#009444")    
  color.names = c("ACC",
                  "BLCA",
                  "BRCA",
                  "CESC",
                  "CHOL",
                  "COAD",
                  "DLBC",
                  "ESCA",
                  "GBM",
                  "HNSC",
                  "KICH",
                  "KIRC",
                  "KIRP",
                  "LAML",
                  "LGG",
                  "LIHC",
                  "LUAD",
                  "LUSC",
                  "MESO",
                  "OV",
                  "PAAD",
                  "PCPG",
                  "PRAD",
                  "READ",
                  "SARC",
                  "SKCM",
                  "STAD",
                  "TGCT",
                  "THCA",
                  "THYM",
                  "UCEC",
                  "UCS",
                  "UVM")
  names(colors) = color.names
  color.scale = scale_color_manual(name="Cancer", values=colors)
  return(color.scale)
}

getPCACancerFill = function() {
  # according to google spreadsheet: https://docs.google.com/spreadsheets/d/1Nb9mMkonAhZR1_2OI9nv4ylCei0LZjAf-2vYTRQcXKw/edit#gid=1704872109
  colors = c(
    "#C1A72F",
    "#FAD2D9",
    "#ED2891",
    "#F6B667",
    "#104A7F",
    "#9EDDF9",
    "#3953A4",
    "#007EB5",
    "#B2509E",
    "#97D1A9",
    "#ED1C24",
    "#F8AFB3",
    "#EA7075",
    "#754C29",
    "#D49DC7",
    "#CACCDB",
    "#D3C3E0",
    "#A084BD",
    "#542C88",
    "#D97D25",
    "#6E7BA2",
    "#E8C51D",
    "#7E1918",
    "#DAF1FC",
    "#00A99D",
    "#BBD642",
    "#00AEEF",
    "#BE1E2D",
    "#F9ED32",
    "#CEAC8F",
    "#FBE3C7",
    "#F89420",
    "#009444")    
  color.names = c("ACC",
                  "BLCA",
                  "BRCA",
                  "CESC",
                  "CHOL",
                  "COAD",
                  "DLBC",
                  "ESCA",
                  "GBM",
                  "HNSC",
                  "KICH",
                  "KIRC",
                  "KIRP",
                  "LAML",
                  "LGG",
                  "LIHC",
                  "LUAD",
                  "LUSC",
                  "MESO",
                  "OV",
                  "PAAD",
                  "PCPG",
                  "PRAD",
                  "READ",
                  "SARC",
                  "SKCM",
                  "STAD",
                  "TGCT",
                  "THCA",
                  "THYM",
                  "UCEC",
                  "UCS",
                  "UVM")
  names(colors) = color.names
  color.scale = scale_fill_manual(name="Cancer", values=colors)
  return(color.scale)
}

getTypeColorScale = function() {
  colors = c("#f0f0f0", "#4575b4", "#1a9850","#d73027") #positive is dark grey       
  color.names = c("None","Direct","Proximal","Clustered")
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="Interaction with Phosphosites", values=colors)
  return(clinical.color.scale)
}

getTypeFillScale = function() {
  colors = c("#f0f0f0", "#4575b4", "#1a9850","#d73027") #positive is dark grey       
  color.names = c("None","Direct","Proximal","Clustered")
  names(colors) = color.names
  clinical.color.scale = scale_fill_manual(name="Interaction with Phosphosites", values=colors)
  return(clinical.color.scale)
}

### reference files ###
# kinome file
manning_kinome_wgene = read.table(header =T, quote = "",sep = '\t', row.names = NULL,file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Gene_family/knome_render_peerj-01-126-s001_wgene.txt")
manning_kinome_wgene_map = manning_kinome_wgene[,9:12]
cancer_gene = as.vector(t(read.table(header=F, file="input/CancerGeneListV5-2014-04-18.add-Rahman.add-Fanconi-Gene.txt")))

# cptac file
cptac_site_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/input/CPTAC.BRCA.OV.phosphoproteome.phosphosite.itraq.wgene.tsv"
cptac_site = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = cptac_site_f)

# hotspot3d cluster file
cluster_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201805/MC3.maf.3D_Proximity.pairwise.3D_Proximity_cleaned.sites.3D_Proximity_cleaned.musites.site.l0.ad10.r10.clusters" # load the ref corrected cluster file
cluster = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, file = cluster_f, colClasses=c(rep('character',3),rep("numeric",4),rep("character",7)))
pass_cluster_f = "output/Data_201805_cc.p0.05.cluster.tsv"
annotated_cluster = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = pass_cluster_f, colClasses=c(rep('character',3),rep("numeric",4),rep("character",7),rep("numeric",4),"character"))
annotated_cluster_h = annotated_cluster[annotated_cluster$Type=="Hybrid",]
