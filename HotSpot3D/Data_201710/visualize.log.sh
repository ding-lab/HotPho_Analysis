# TODO: pending update to this version
# get presence
perl ~/bin/hotspot3d_KH/scripts/clusterPDBPresence_musite.pl 20171013.search.dc20.musites onlyCanonical.20171013.cluster.ad10.r10.net.rec.strInd.subInd.uns.clusters 20171013_musite

# plot some of the previous top candidates
musite=20171013.search.dc20.musites
cluster=onlyCanonical.20171013.cluster.ad10.r10.net.rec.strInd.subInd.uns.clusters

# limit cluster file to clusters of interest first

# visualize
## cancer phosphosites
grep1 734.10 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=1H8F --output-file=pml_scripts/GSK3B_734.10_1H8F.pml --script-only 
grep1 8977.2 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=1QZJ --output-file=pml_scripts/KIT_8977.2_1QZJ.pml --script-only 
grep1 5090.3 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=1M14 --output-file=pml_scripts/EGFR_5090.3_1M14.pml --script-only 
grep1 1967.0 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=4E4X --output-file=pml_scripts/BRAF_1967.0_4E4X.pml --script-only 

## mutations with highly expressed protein/phosphoprotein
grep1 5086.0 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=1YY9 --output-file=pml_scripts/EGFR_5086.0_1YY9.pml --script-only 
grep1 5090.1 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=3KEX --output-file=pml_scripts/ERBB3_5090.1_3KEX.pml --script-only 
grep1 8977.0 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=1QZJ --output-file=pml_scripts/KIT_8977.0_1QZJ.pml --script-only

## phosphosites with functional mutations
grep1 734.0 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=4EJN --output-file=pml_scripts/AKT1_734.0_4EJN.pml --script-only 
grep1 2572.0 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=1LD2 --output-file=pml_scripts/CDK4_2572.0_1LD2.pml --script-only 
grep1 5090.0 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=1M14 --output-file=pml_scripts/EGFR_5090.0_1M14.pml --script-only
grep1 5483.10 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=1OVC --output-file=pml_scripts/ERBB2_5483.1_1OVC.pml --script-only 
grep1 12478.2 ${cluster} > tmp.clusters
hotspot3d visual --musites-file=${musite} --clusters-file=tmp.clusters --pdb=3HHM --output-file=pml_scripts/PIK3CA_12478.2_3HHM.pml --script-only