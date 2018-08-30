##### preprocess_cluster.R #####
# Kuan-lin Huang @ WashU 2017 May; updated 2018 April

### dependencies ###
bdir = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201805"
setwd(bdir)

# cptac file
cptac_site_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/input/CPTAC.BRCA.OV.phosphoproteome.phosphosite.itraq.wgene.tsv"
cptac_site = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = cptac_site_f)


# some REF residues for sites seem to be off during generation of pairwise or cluster file
# we found this is largely due to mis-reporting on mapping in PDB
# thankfully mostly non-important proteins
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

##### pairwise #####

site_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201807/3D_Proximity.musites.gz"
site = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(site_f))

cat("number of musite pairs dis-qualified from comparing to original phosphosites\n")
table(paste(site$Transcript2,site$Site2) %in% paste(phosphosites_map$Transcript,phosphosites_map$originalLabel))
site_cleaned = site[paste(site$Transcript2,site$Site2) %in% paste(phosphosites_map$Transcript,phosphosites_map$originalLabel),]
write.table(site_cleaned, quote=F, sep="\t", file = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201807/3D_Proximity_cleaned.musites", row.names = F)

sitesite_f = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201807/3D_Proximity.sites.gz"
sitesite = read.table(header=T, quote = "", sep="\t", stringsAsFactors = F, fill =T, file = gzfile(sitesite_f))

cat("number of site-site pairs dis-qualified from comparing to original phosphosites\n")
table(paste(sitesite$Transcript2,sitesite$Site2) %in% paste(phosphosites_map$Transcript,phosphosites_map$originalLabel) & paste(sitesite$Transcript1,sitesite$Site1) %in% paste(phosphosites_map$Transcript,phosphosites_map$originalLabel))
sitesite_cleaned = sitesite[paste(sitesite$Transcript2,sitesite$Site2) %in% paste(phosphosites_map$Transcript,phosphosites_map$originalLabel) & paste(sitesite$Transcript1,sitesite$Site1) %in% paste(phosphosites_map$Transcript,phosphosites_map$originalLabel),]
write.table(sitesite_cleaned, quote=F, sep="\t", file = "/Users/khuang/Box\ Sync/Ding_Lab/Projects_Current/hotpho_data/HotSpot3D/Data_201807/3D_Proximity_cleaned.sites", row.names = F)
