##### preprocess_cluster.R #####
# Kuan-lin Huang @ WashU 2017 May; updated 2018 April

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data"
setwd(bdir)
source("Analyses/hotpho_analyses_functions.R")

### common dependencies for plotting ###

# hotspot3d cluster file
cluster_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201807/MC3.maf.3D_Proximity.pairwise.3D_Proximity_cleaned.sites.3D_Proximity_cleaned.musites.recurrence.l0.ad10.r10.clusters"

##### CLUSTERs #####

cluster = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = cluster_f, colClasses=c(rep('character',3),rep("numeric",4),rep("character",7)))
colnames(cluster)=gsub("\\.","_",colnames(cluster)) # period cuz troubles for R data frame manipulation 

cluster_summary = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = paste(cluster_f,"summary",sep="."), colClasses=c(rep('character',2),rep("numeric",13),rep("character",7)))

### clean up clusters ###
# cluster$Alt_Mutation_Gene = gsub(".*:","",cluster$Alternative_Transcripts) # always get the last one as that seems to be the correct one; ex. ENST00000320868:p.T42|ENST00000320868:p.Y42
# cluster$Original_Mutation_Gene = cluster$Mutation_Gene
# cluster$Mutation_Gene[cluster$Alt_Mutation_Gene != ""] = cluster$Alt_Mutation_Gene[cluster$Alt_Mutation_Gene != ""]
# #write.table(cluster, quote=F, sep="\t", file = "HotSpot3D/Data_201807/convertToRef.filtered.0705.pass.fast.MC3.combined.mumu.musite.sitesite.max20.ld0.ad10.r10.net.recur.unspec.strInd.subInd.noSingletons.clusters", row.names = F)

# some REF residues for sites seem to be off during generation of pairwise or cluster file
# thankfully hotspot3d only considers position
# post-hoc fixing them here, the PTMs
PTMcosmo_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/input/PTMPortal_ptm_sites_dump.phospho_only.with_enst.tsv"
PTMcosmo = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = PTMcosmo_f)
PTMcosmo$Mutation_Gene = paste("p.",PTMcosmo$residue,PTMcosmo$p_coord,sep="")
PTMcosmo_map = PTMcosmo[,c("Ensembl.Transcript","p_coord","Mutation_Gene")]
cptac_site$Position = gsub("p.[A-Z]","",cptac_site$amino_acid_residue)
cptac_site_map = cptac_site[,c("ensembl_transcript_id","Position","amino_acid_residue")]
colnames(PTMcosmo_map) = c("Transcript","Position","originalLabel")
colnames(cptac_site_map) = c("Transcript","Position","originalLabel")
phosphosites_map = rbind(cptac_site_map,PTMcosmo_map)
phosphosites_map = phosphosites_map[!duplicated(paste(phosphosites_map$Transcript,phosphosites_map$Position)),]

cluster$Mutation_Gene = gsub("p. ","p.",cluster$Mutation_Gene) # some residues had a gap...
cluster$Position = gsub("p.[A-Z]*([0-9]*)[A-Z]*","\\1",cluster$Mutation_Gene)
cluster_merge = merge(cluster,phosphosites_map,by=c("Transcript","Position"),all.x=T)
cat("Number of residues with inconsistent residues from the original phosphosite files (True being inconsistent):\n")
table(cluster_merge$Mutation_Gene[cluster_merge$Alternate=="ptm"] != cluster_merge$originalLabel[cluster_merge$Alternate=="ptm"])
# cluster_merge$Mutation_Gene[cluster_merge$Alternate=="ptm"] = cluster_merge$originalLabel[cluster_merge$Alternate=="ptm"]
# table(cluster_merge$Mutation_Gene[cluster_merge$Alternate=="ptm"] != cluster_merge$originalLabel[cluster_merge$Alternate=="ptm"])
# cluster_merge = cluster_merge[,-c(which(colnames(cluster_merge) =="originalLabel"))]
# write.table(cluster_merge, quote=F, sep="\t", file = "HotSpot3D/Data_201807/PTM_MC3_noFs.maf.3D_Proximity.pairwise.3D_Proximity.sites.3D_Proximity.musites.site.l0.ad10.r10.cleaned.clusters", row.names = F)

cat("Unique clusters (unfiltered):",length(unique(cluster$Cluster)),"\n")
cat("Unique genes:",length(unique(cluster$Gene_Drug)),"\n")

annotated_cluster = annotate_cluster(cluster)
cat("Hybrid clusters (unfiltered):",length(unique(annotated_cluster$Cluster[annotated_cluster$Type=="Hybrid"])),"\n")
cat("Mutation-only clusters (unfiltered):",length(unique(annotated_cluster$Cluster[annotated_cluster$Type=="Mut_Only"])),"\n")
cat("Site-only clusters (unfiltered):",length(unique(annotated_cluster$Cluster[annotated_cluster$Type=="Site_Only"])),"\n")

# annotated_clusterPTM$ref_org = gsub("p.([A-Z])[0-9]+","\\1",annotated_clusterPTM$Original_Mutation_Gene)
# table(annotated_clusterPTM$ref_org)

annotated_cluster_centroids = annotated_cluster[annotated_cluster$Geodesic_From_Centroid==0,]
annotated_cluster_centroids_unique = annotated_cluster_centroids[!duplicated(annotated_cluster_centroids$Cluster),]

# Use cluster closeness to find the top 5% clusters
# hybrid_clus = unique(annotated_cluster$Cluster[annotated_cluster$Type == "Hybrid"])
# mut_clus = unique(annotated_cluster$Cluster[annotated_cluster$Type == "Mut_Only"])
# site_clus = unique(annotated_cluster$Cluster[annotated_cluster$Type == "Site_Only"])
# 
top_clust = c()
cluster_types = c("Hybrid","Mut_Only","Site_Only")
thres = 0.95
h_thres = 0
cluster_type_thres = list()
cluster_summary$Type = NA

for (type in cluster_types){
  cluster_summary$Type[cluster_summary$Cluster_ID %in% unique(annotated_cluster$Cluster[annotated_cluster$Type == type])] = type
  clust_type = cluster_summary[cluster_summary$Cluster_ID %in% 
                                 unique(annotated_cluster$Cluster[annotated_cluster$Type == type]),]
  type_thres = quantile(clust_type$Cluster_Closeness,probs=thres)
  cluster_type_thres[type]=type_thres
  if (type =="Hybrid"){h_thres= type_thres}
  top_clust = c(top_clust,clust_type$Cluster_ID[clust_type$Cluster_Closeness > type_thres])
  cat(type, ":\t",type_thres, "\n")
}
annotated_cluster_centroids_unique_pass = annotated_cluster_centroids_unique[annotated_cluster_centroids_unique$Cluster %in% top_clust,]

# # get the 5% most significant clusters for each category for subsequent analysis
# cluster_types = c("Hybrid","Mut_Only","Site_Only")
# thres = 0.95
# h_thres = 0
# cluster_type_thres = list()
# for (type in cluster_types){
#   type_thres = quantile(annotated_cluster_centroids_unique$Closeness_Centrality[annotated_cluster_centroids_unique$Type == type],probs=thres)
#   cluster_type_thres[type]=type_thres
#   if (type =="Hybrid"){h_thres= type_thres}
# }
# annotated_cluster_centroids_unique_pass = annotated_cluster_centroids_unique[
#   annotated_cluster_centroids_unique$Type == "Hybrid" & annotated_cluster_centroids_unique$Closeness_Centrality > cluster_type_thres[["Hybrid"]] |
#     annotated_cluster_centroids_unique$Type == "Mut_Only" & annotated_cluster_centroids_unique$Closeness_Centrality > cluster_type_thres[["Mut_Only"]] |
#     annotated_cluster_centroids_unique$Type == "Site_Only" & annotated_cluster_centroids_unique$Closeness_Centrality > cluster_type_thres[["Site_Only"]], 
#   ]


# take a look at the density
p = ggplot(cluster_summary,aes(x = log10(Cluster_Closeness), fill=Type))
p = p + geom_density(alpha=0.2,size=0.5)
p = p + theme_bw() #+ xlim(0,5)
p = p + geom_vline(xintercept = log10(h_thres),alpha=0.5)
p
fn = paste("output/Data_201807_cc_dist_by_cluster_type.pdf",sep="_")
ggsave(fn, useDingbat=F)

cat("\n")
cat("Hybrid clusters (filtered):",length(unique(annotated_cluster_centroids_unique_pass$Cluster[annotated_cluster_centroids_unique_pass$Type=="Hybrid"])),"\n")
rank_vectors(annotated_cluster_centroids_unique_pass$Gene_Drug[annotated_cluster_centroids_unique_pass$Type=="Hybrid"])
cat("Mutation-only clusters (filtered):",length(unique(annotated_cluster_centroids_unique_pass$Cluster[annotated_cluster_centroids_unique_pass$Type=="Mut_Only"])),"\n")
rank_vectors(annotated_cluster_centroids_unique_pass$Gene_Drug[annotated_cluster_centroids_unique_pass$Type=="Mut_Only"])
cat("Site-only clusters (filtered):",length(unique(annotated_cluster_centroids_unique_pass$Cluster[annotated_cluster_centroids_unique_pass$Type=="Site_Only"])),"\n")
rank_vectors(annotated_cluster_centroids_unique_pass$Gene_Drug[annotated_cluster_centroids_unique_pass$Type=="Site_Only"])

annotated_cluster_pass = annotated_cluster[annotated_cluster$Cluster %in% annotated_cluster_centroids_unique_pass$Cluster, ]
write.table(annotated_cluster_pass, quote=F, sep="\t", file = "output/Data_201807_cc.p0.05.cluster.tsv", row.names = F)

# sync up transcripts within the same cluster; when the PTM sites are on a different transcript
transvarIn_f = "HotSpot3D/Data_201807/PTM_Site_transvar.txt.gz"
transvarIn = read.table(header=F, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(transvarIn_f))
transvarIn_anno = transvarIn[,c(3,6,12)]
colnames(transvarIn_anno) = c("Transcript","Position","GenomicPosition")
# transvarIn_anno$Start = gsub("chr.*g.([0-9]*)_.*","\\1",transvarIn_anno$GenomicPosition)
# transvarIn_anno$Stop = gsub("chr.*g.([0-9]*)_([0-9]*)/.*","\\2",transvarIn_anno$GenomicPosition)

annotated_cluster_pass = merge(annotated_cluster_pass, transvarIn_anno, by=c("Transcript","Position"), all.x=T)
annotated_cluster_pass$Start[annotated_cluster_pass$Alternate == "ptm"] = gsub("chr.*g.([0-9]*)_.*","\\1",annotated_cluster_pass$GenomicPosition[annotated_cluster_pass$Alternate == "ptm"])
annotated_cluster_pass$Stop[annotated_cluster_pass$Alternate == "ptm"] = gsub("chr.*g.([0-9]*)_([0-9]*)/.*","\\2",annotated_cluster_pass$GenomicPosition[annotated_cluster_pass$Alternate == "ptm"])

transvar_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/output/annotated_cluster_h_PTM_transvarOut.txt"
transvar = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = transvar_f)
transvar$Transcript = gsub(" .*","",transvar$transcript)
transvar$Mutation_Gene = gsub(".*/(p.[A-Z][0-9]+)","\\1",transvar$coordinates.gDNA.cDNA.protein.)
transvar_anno = transvar[grep("p.",transvar$Mutation_Gene),c("input","Transcript","Mutation_Gene")]
transvar_anno$Start = gsub(".*:g.([0-9]+)_.*","\\1",transvar_anno$input)

cluster = "10394.1"
for (cluster in annotated_cluster_pass$Cluster){
  annotated_cluster_pass_c = annotated_cluster_pass[annotated_cluster_pass$Cluster == cluster,]
  if (annotated_cluster_pass_c$Type != "Site_Only"){
    if (length(unique(annotated_cluster_pass_c$Transcript))>1){
      mutTranscript = annotated_cluster_pass_c$Transcript[annotated_cluster_pass_c$Alternate!="ptm"][1]
      PTMs = annotated_cluster_pass_c$Mutation_Gene[annotated_cluster_pass_c$Alternate=="ptm"]
      for (PTM in PTMs){
        if (annotated_cluster_pass_c$Transcript[annotated_cluster_pass_c$Mutation_Gene == PTM] != mutTranscript){
          
          updatedSite = gsub("(p.[A-Z][0-9]+).*","\\1",transvar_anno$Mutation_Gene[transvar_anno$Start == annotated_cluster_pass_c$Start[annotated_cluster_pass_c$Mutation_Gene == PTM] &
                                                                                     transvar_anno$Transcript == mutTranscript])
          if (length(updatedSite>0)){
            annotated_cluster_pass$Transcript[annotated_cluster_pass$Cluster == cluster & annotated_cluster_pass$Mutation_Gene == PTM] = mutTranscript
            annotated_cluster_pass$Mutation_Gene[annotated_cluster_pass$Cluster == cluster & annotated_cluster_pass$Mutation_Gene == PTM] = updatedSite
          }
        }
      }
    }
  }
}

annotated_cluster_pass$Position = gsub("p.[A-Z]([0-9]+).*","\\1",annotated_cluster_pass$Mutation_Gene)
write.table(annotated_cluster_pass, quote=F, sep="\t", file = "output/Data_201807_cc.p0.05.cluster_transcriptSynced.tsv", row.names = F)
