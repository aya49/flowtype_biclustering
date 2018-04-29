# Input: protein gene interaction network --> Output: distance between genes (only those in largest connected component)
# aya43@sfu.ca 20151228

## root directory
root = "~/projects/network/innatedb"
result_dir = "result"; dir.create(result_dir, showWarnings=F)
setwd(root)

rootD = "~/projects/IMPC"
panelL = c("P1")
centreL = c("Sanger_SPLEEN")#,"Sanger_MLN","CIPHE","TCP","H")
ci = 1; panel = panelL[ci]; centre = centreL[ci]

resultD_dir = paste0(rootD,"/result/", panel, "/", centre)

## input directories
# Data format is IRefTab: http://irefindex.org/wiki/index.php?title=README_MITAB2.6_for_iRefIndex#What_each_line_represents
data_dir = paste0("all.mitab") # All Imported Experimentally Validated Interactions (don't just get rows for mouse only!)
# data_dir = paste0("innatedb_all.mitab") #all

meta_dir = paste0(resultD_dir,"/meta")
meta_file_dir = paste0(meta_dir,"/file")


## output directories
edge_dir = paste(result_dir, "/edge", sep=""); dir.create(edge_dir, showWarnings=F)
DSD_dir = paste(result_dir, "/dist_DSD", sep=""); dir.create(DSD_dir, showWarnings=F)


## libraries
library("stringr")
library("colorspace")
library("lubridate") #if there are date variables
library("data.table") #to read tables
library("robustHD")
library("rols") #ontology look up service
library("GO.db") #GO only
library("iRefR") #MITAB format parser
library("rentrez") #NCBI database look up service
library("RISmed")
library("taxize")
library("igraph")
library("arules") #discretize
library("foreach")
library("doMC")
# devtools::install_github("altayg/ganet")
# library(ganet)
# devtools::install_github("R3CPET","sirusb")
# library(R3CPET)

source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")




#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)




## Options
options(na.strings = c("","-", " ", "?", "NULL", NA))
options(stringsAsFactors = F)

control = "[+]_[+]|[+]_Y"
kt = km = kp = 1/3 #paramters for calculating edge weights










## get genes
meta_file = get(load(paste0(meta_file_dir,".Rdata")))
gene0 = unique(meta_file$gene)
gene = gene0[!grepl(control, gene0)]
# gene = str_split(gene0,"_")
# gene1 = sapply(gene, function(x) x[1])
# gene2 = sapply(gene, function(x) x[2])
# identical(gene1,gene2)
# sort(unique(gene1[gene2=="+"]))
# for (g in unique(gene1[gene2=="+"])) cat("\n(",sum(gene1==g & gene2=="+"),"/",sum(gene1==g),") ", g, sep="")
# sort(unique(gene1[gene2=="Y"]))
# for (g in unique(gene1[gene2=="Y"])) cat("\n(",sum(gene1==g & gene2=="Y"),"/",sum(gene1==g),") ", g, sep="")
# gene = unique(gene1[gene2!="+" & gene2!="Y"]) #delete 118 genes because mostly 1/2 or 1/1 -/+
# # gene = gene1
# gene = str_split(gene,"[(]") #take out "(b)" at the end of gene name, number of genes dont decrease
# gene = sapply(gene, function(x) x[1])
# gene = toupper(gene)
# # write.table(gene,file=gene_dir,sep="\t",row.names=F,col.names=F,quote=F)







## Ontology ---------------------------------

## HTTP 404 MI ontology not found :(
mi <- Ontology("mi")
miterm <- term(mi,"MI:0794")

qry = olsSearch(OlsSearch(q = gene[1]))

entrez_dbs() #databases in NCBI






## Import network & Specify column names ----------------------------------

b0 = data.frame(fread(input=data_dir))
colnames(b0)



## delete columns
#where all values are the same
b0 = delcols(b0)


aliasA_col = "alias_A"
aliasB_col = "alias_B"

interactA_col = "X.unique_identifier_A"
interactB_col = "unique_identifier_B"

taxidA_col = "ncbi_taxid_A"
taxidB_col = "ncbi_taxid_B"

conf_col = "confidence_score"
intertype_col = "interaction_type"
intermethod_col = "interaction_detection_method"

taxidM = "Mouse"
taxidH = "Human"
taxid = ""
# taxid = taxidM

# "X.unique_identifier_A"               "unique_identifier_B"   e.g. innatedb:IDBG-90895
# "alt_identifier_A"                    "alt_identifier_B"      e.g. innatedb:IDBG-90895
# "alias_A"                             "alias_B"               e.g. uniprotkb:CASP9_HUMAN|refseq:NP_127463|uniprotkb:P55211|refseq:NP_001220|refseq:NP_001264983|hgnc:CASP9(display_short), refseq:NP_982281|uniprotkb:PPIE_HUMAN|uniprotkb:Q9UNP9|refseq:NP_006103|refseq:NP_001181936|hgnc:PPIE(display_short), uniprotkb:AUHM_HUMAN|refseq:NP_001689|uniprotkb:Q13825|hgnc:AUH(display_short)
# "interaction_detection_method"                                e.g. psi-mi:"MI:0114"(x-ray crystallography), psi-mi:"MI:0018"(two hybrid)
# "author"                              "pmid"                  e.g. Fair et al.(2001), Renatus et al. (2001); pubmed:11313484
# "ncbi_taxid_A"                        "ncbi_taxid_B"          e.g. taxid:9606(Human)
# !!missing interactorA/B
# "interaction_type"                                            e.g. psi-mi:"MI:0915"(physical association)
# "source_database"                                             e.g. MI:0463(biogrid)
# "idinteraction_in_source_db"                                  e.g. BIOGRID:10022
# "confidence_score"                                            e.g. lpr:6|hpr:6|np:1|
# 
# EXTRA COLUMNS:
# "expansion_method"                                            e.g. "-"     "spoke"; Model used to convert n-ary data into binary data for purpose of export in MITAB file
# "biological_role_A"                   "biological_role_B"     e.g. psi-mi:"MI:0499"(unspecified role)
# "exp_role_A"                          "exp_role_B"            e.g. psi-mi:"MI:0499"(unspecified role); psi-mi:"MI:0498"(prey), "MI:0496"(bait); "psi-mi:\"MI:0581\"(suppressor gene)"  
# !!"interactor_type_A"                 "interactor_type_B"     e.g. psi-mi:"MI:0326"(protein); "psi-mi:\"MI:0319\"(deoxyribonucleic acid)" "psi-mi:\"MI:0320\"(ribonucleic acid)"
# "annotations_interaction"                                     e.g. not used by IRefIndex; comment:"BIND interaction division: BIND-3DBP; "
# "ncbi_taxid_host_organism"                                    e.g. may differ from taxid; taxid:-1
# "parameters_interaction"                                      e.g. not used by IRefIndex
# "creation_date"                       "update_date"           e.g. 2014/11/29
# "participant_identification_method_A" "participant_identification_method_B"            e.g. psi-mi:"MI:0363"(inferred by author), psi-mi:"MI:0396"(predetermined participant)


## delete rows: see data regarding taxonomy -----------------------------------
sum(b0[,taxidA_col]!=b0[,taxidB_col])
sum(b0[,taxidA_col]!=b0[,taxidB_col] &
      (!grepl(taxid,b0[,taxidA_col]) & 
         !grepl(taxid,b0[,taxidB_col]))) # where two interactors have different taxonomy and they are not mouse or human

bmouse = grepl(taxid,(b0[,taxidA_col]),ignore.case=T) | grepl(taxid,(b0[,taxidB_col]),ignore.case=T)
b00 = b0
b1 = b0 = b00[bmouse,]

cat("\n# deleted rows: ", sum(!bmouse), "/", nrow(b0), sep="")




## create edge list ---------------------------------------------

# make ordering the same, identify non/duplicates
# edge list with interactor names that is lookup-able on GO
b0a = gsub("[:]|[(]display_short[)]","",str_extract(b0[,aliasA_col],":(?!.*:).*?[(]display_short[)]")) #get the thing between ":" and "(display_short)"
b0a = toupper(b0a)
names(b0a) = NULL
sum(!is.na(match(gene,b0a))) #number of gene matches between network and data

b0b = gsub("[:]|[(]display_short[)]","",str_extract(b0[,aliasB_col],":(?!.*:).*?[(]display_short[)]")) #get the thing between ":" and "(display_short)"
b0b = toupper(b0b)
names(b0b) = NULL
sum(!is.na(match(gene,b0b))) #number of gene matches between network and data

b0edgel = cbind(b0a,b0b)
sum(duplicated(b0edgel) | duplicated(b0edgel, fromLast = TRUE)) #see if there are any duplicate edges

b1edgel = cbind(b0a,b0b)
b1abi = b1edgel[,1]>b1edgel[,2]
b1edgel[b1abi,] = cbind(b1edgel[b1abi,2],b1edgel[b1abi,1])







## merge rows: https://academic.oup.com/database/article/2433131
# number of publications
# interaction type:
# - scv1=0.10 || cv1=MI:0208 | genetic interaction
# - scv2=0.33 || cv2=MI:0403 | colocalization
# - scv3=0.33 || cv3=MI:0914 | association
# - scv4=0.66 || cv4=MI:0915 | physical association
# - scv5=1.00 || cv5=MI:0407 | direct interaction
# - scv6=0.05 || cv6=unknown | unknown
# - Gscv1=scv1 | Gscv2=scv2 | Gscv3=scv3, scv4, scv5
# interaction detection method:
# - scv1 = 1.00 || cv1 = MI:0013 | biophysical 
# - scv2 = 0.66 || cv2 = MI:0090 | protein complementation assay 
# - scv3 = 0.10 || cv3 = MI:0254 | genetic interference 
# - scv4 = 0.10 || cv4 = MI:0255 | post transcriptional interference 
# - scv5 = 1.00 || cv5 = MI:0401 | biochemical 
# - scv6 = 0.33 || cv6 = MI:0428 | imaging technique 
# - scv7 = 0.05 || cv7 = unknown | unknown 
# - Gscv1 = scv1 | Gscv2 = scv2 | Gscv3 = scv3 | Gscv4 = scv4|Gscv5 = scv5 | Gscv6 = scv6




## Publication scores
b1np = as.numeric(gsub("np[:]","",str_extract(b1[,conf_col],"np[:][0-9]+")))




## Ineraction type scores
b1it = b1[,intertype_col]
sort(unique(b1it))
b1it[b1it%in%c("psi-mi:\"MI:0208\"(genetic interaction)", #? not sure
               "psi-mi:\"MI:0794\"(synthetic genetic interaction defined by inequality)",
               "psi-mi:\"MI:0796\"(suppressive genetic interaction defined by inequality)",
               "psi-mi:\"MI:0799\"(additive genetic interaction defined by inequality)")] = "0.1"
b1it[b1it%in%c("psi-mi:\"MI:0403\"(colocalization)")] = "0.33"
b1it[b1it%in%c("psi-mi:\"MI:0407\"(direct interaction)",
               "psi-mi:\"MI:1126\"(self interaction)",
               "psi-mi:\"MI:0196\"(covalent interaction)", #MI0195 covalent binding
               "psi-mi:\"MI:0195\"(covalent binding)",
               "psi-mi:\"MI:0408\"(disulfide bond)",
               "psi-mi:\"MI:0414\"(enzymatic reaction)",
               "psi-mi:\"MI:0566\"(sumoylation reaction)",
               "psi-mi:\"MI:1250\"(isomerase reaction)",
               "psi-mi:\"MI:0192\"(acetylation reaction)",
               "psi-mi:\"MI:0220\"(ubiquitination reaction)",
               "psi-mi:\"MI:0204\"(deubiquitination reaction)",
               "psi-mi:\"MI:0220\"(ubiquitination reaction)",
               "psi-mi:\"MI:0213\"(methylation reaction)",
               "psi-mi:\"MI:0194\"(cleavage reaction)",
               "psi-mi:\"MI:0570\"(protein cleavage)",
               "psi-mi:\"MI:0217\"(phosphorylation reaction)",
               "psi-mi:\"MI:0203\"(dephosphorylation reaction)")] = "1"
b1it[b1it%in%c("psi-mi:\"MI:0914\"(association)")] = "0.33"
b1it[b1it%in%c("psi-mi:\"MI:0915\"(physical association)")] = "0.66"
b1it = as.numeric(b1it)
b1it[is.na(b1it)] = 0 #"", "psi-mi:"MI:0179"(other modification)"

# b1itG = b1it #supercategory
# b1itG[b1$interaction_type%in%c("psi-mi:\"MI:0914\"(association)",
# "psi-mi:\"MI:0407\"(direct interaction)",
# "psi-mi:\"MI:0915\"(physical association)")] = "0.33"


## Interaction detection method scores
b1im = b1[,intermethod_col]
sort(unique(b1im))
b1im[b1im%in%c("psi-mi:\"MI:0401\"(biochemical)",
               "psi-mi:\"MI:2197\"(probe interaction assay)",
               "psi-mi:\"MI:2198\"(labelling assay)",
               "psi-mi:\"MI:1313\"(proximity labelling technology)",
               "psi-mi:\"MI:1314\"(proximity-dependent biotin identification)",
               "psi-mi:\"MI:1036\"(nucleotide exchange assay)",
               "psi-mi:\"MI:0949\"(gdp/gtp exchange assay)",
               "psi-mi:\"MI:0027\"(cosedimentation)",
               "psi-mi:\"MI:0028\"(cosedimentation in solution)",
               "psi-mi:\"MI:0029\"(cosedimentation through density gradient)",
               "psi-mi:\"MI:0417\"(footprinting)",
               "psi-mi:\"MI:0605\"(enzymatic footprinting)",
               "psi-mi:\"MI:1183\"(nuclease footprinting)",
               "psi-mi:\"MI:0606\"(DNase I footprinting)",
               "psi-mi:\"MI:0030\"(cross-linking study)",
               "psi-mi:\"MI:0031\"(protein cross-linking with a bifunctional reagent)",
               "psi-mi:\"MI:0430\"(nucleic acid uv cross-linking assay)",
               "psi-mi:\"MI:0982\"(electrophoretic mobility-based method)",
               "psi-mi:\"MI:0412\"(electrophoretic mobility supershift assay)",
               "psi-mi:\"MI:0807\"(comigration in gel electrophoresis)",
               "psi-mi:\"MI:0404\"(comigration in non denaturing gel electrophoresis)",
               "psi-mi:\"MI:0276\"(blue native page)",
               "psi-mi:\"MI:0413\"(electrophoretic mobility shift assay)",
               "psi-mi:\"MI:0808\"(comigration in sds page)",
               "psi-mi:\"MI:0091\"(chromatography technology)",
               "psi-mi:\"MI:0226\"(ion exchange chromatography)",
               "psi-mi:\"MI:0047\"(far western blotting)",
               "psi-mi:\"MI:0415\"(enzymatic study)",
               "psi-mi:\"MI:0997\"(ubiquitinase assay)",
               "psi-mi:\"MI:0434\"(phosphatase assay)",
               "psi-mi:\"MI:0841\"(phosphotransferase assay)",
               "psi-mi:\"MI:0424\"(protein kinase assay)",
               "psi-mi:\"MI:1005\"(adp ribosylase assay)",
               "psi-mi:\"MI:0889\"(acetylase assay)",
               "psi-mi:\"MI:0887\"(histone acetylase assay)",
               "psi-mi:\"MI:0406\"(deacetylase assay)",
               "psi-mi:\"MI:0870\"(demethylase assay)",
               "psi-mi:\"MI:0515\"(methyltransferase assay)",
               "psi-mi:\"MI:0516\"(methyltransferase radiometric assay)",
               "psi-mi:\"MI:0990\"(cleavage assay)",
               "psi-mi:\"MI:0071\"(molecular sieving)",
               "psi-mi:\"MI:0400\"(affinity technology)",
               "psi-mi:\"MI:0405\"(competition binding)",
               "psi-mi:\"MI:0034\"(display technology)",
               "psi-mi:\"MI:0084\"(phage display)",
               "psi-mi:\"MI:0066\"(lambda phage display)",
               "psi-mi:\"MI:0048\"(filamentous phage display)",
               "psi-mi:\"MI:0440\"(saturation binding)",
               "psi-mi:\"MI:0892\"(solid phase assay)",
               "psi-mi:\"MI:0047\"(far western blotting)",
               "psi-mi:\"MI:0411\"(enzyme linked immunosorbent assay)",
               "psi-mi:\"MI:0049\"(filter binding)",
               "psi-mi:\"MI:0008\"(array technology)",
               "psi-mi:\"MI:0081\"(peptide array)",
               "psi-mi:\"MI:0678\"(antibody array)",
               "psi-mi:\"MI:0695\"(sandwich immunoassay)",
               "psi-mi:\"MI:0089\"(protein array)",
               "psi-mi:\"MI:0004\"(affinity chromatography technology)",
               "psi-mi:\"MI:0729\"(luminescence based mammalian interactome mapping)",
               "psi-mi:\"MI:0676\"(tandem affinity purification)",
               "psi-mi:\"MI:0019\"(coimmunoprecipitation)",
               "psi-mi:\"MI:0007\"(anti tag coimmunoprecipitation)",
               "psi-mi:\"MI:0096\"(pull down)",
               "psi-mi:\"MI:0059\"(gst pull down)",
               "psi-mi:\"MI:1313\"(bioid)")] = "1"
b1im[b1im%in%c("psi-mi:\"MI:0055\"(fluorescent resonance energy transfer)", #biophysical
               "psi-mi:\"MI:0114\"(x-ray crystallography)")] = "1"
b1im[b1im%in%c("psi-mi:\"MI:0090\"(protein complementation assay)",
               "psi-mi:\"MI:0018\"(two hybrid)",
               "psi-mi:\"MI:0398\"(two hybrid pooling approach)",
               "psi-mi:\"MI:0728\"(gal4 vp16 complementation)",
               "psi-mi:\"MI:0399\"(two hybrid fragment pooling approach)")] = "0.66"
b1im[b1im%in%c("psi-mi:\"MI:0428\"(imaging technique)",
               "psi-mi:\"MI:0663\"(confocal microscopy)",
               "psi-mi:\"MI:0040\"(electron microscopy)",
               "psi-mi:\"MI:0020\"(transmission electron microscopy)",
               "psi-mi:\"MI:0410\"(electron tomography)")] = "0.33"
b1im[b1im%in%c("psi-mi:\"MI:0254\"(genetic interference)")] = "0.1"
b1im[b1im%in%c("psi-mi:\"MI:0686\"(unspecified method)")] = "0.05"
b1im[b1im%in%c("psi-mi:\"MI:0203\"(dephosphorylation reaction)",
               "psi-mi:\"MI:0414\"(enzymatic reaction)",
               "psi-mi:\"MI:0194\"(cleavage reaction)",
               "psi-mi:\"MI:0204\"(deubiquitination reaction)",
               "psi-mi:\"MI:0220\"(ubiquitination reaction)")] = "0.66"
b1im[b1im%in%c("psi-mi:\"MI:0001\"(interaction detection method)",
               "psi-mi:\"MI:0045\"(experimental interaction detection)",
               "psi-mi:\"MI:0010\"(beta galactosidase complementation)",
               "psi-mi:\"MI:0255\"(post transcriptional interference)",
               "psi-mi:\"MI:0231\"(mammalian protein protein interaction trap)",
               "psi-mi:\"MI:0011\"(beta lactamase complementation)",
               "psi-mi:\"MI:0232\"(transcriptional complementation assay)",
               "psi-mi:\"MI:0112\"(ubiquitin reconstruction)",
               "psi-mi:\"MI:0113\"(western blot)",
               "psi-mi:\"MI:0630\"(3d-structure)",
               "psi-mi:\"MI:0686\"(unspecified method)",
               "psi-mi:\"MI:0432\"(one hybrid)",
               "psi-mi:\"MI:0256\"(rna interference)",
               "psi-mi:\"MI:0254\"(genetic interference)",
               "psi-mi:\"MI:0813\"(proximity ligation assay)",
               "psi-mi:\"MI:0363\"(inferred by author)",
               "psi-mi:\"MI:0090\"(protein complementation assay)")] = "0"
b1im[b1im%in%c("psi-mi:\"MI:0013\"(biophysical)",
               "psi-mi:\"MI:0114\"(x-ray crystallography)",
               "psi-mi:\"MI:0043\"(electron resonance)",
               "psi-mi:\"MI:0042\"(electron paramagnetic resonance)",
               "psi-mi:\"MI:1247\"(microscale thermophoresis)",
               "psi-mi:\"MI:0016\"(circular dichroism)",
               "psi-mi:\"MI:0051\"(fluorescence technology)",
               "psi-mi:\"MI:0510\"(homogeneous time resolved fluorescence)",
               "psi-mi:\"MI:0420\"(kinase homogeneous time resolved fluorescence)",
               "psi-mi:\"MI:0053\"(fluorescence polarization spectroscopy)",
               "psi-mi:\"MI:0416\"(fluorescence microscopy)",
               "psi-mi:\"MI:0809\"(bimolecular fluorescence complementation)",
               "psi-mi:\"MI:0054\"(fluorescence-activated cell sorting)",
               "psi-mi:\"MI:0115\"(yeast display)",
               "psi-mi:\"MI:0017\"(classical fluorescence spectroscopy)",
               "psi-mi:\"MI:0012\"(bioluminescence resonance energy transfer)",
               "psi-mi:\"MI:0077\"(nuclear magnetic resonance)",
               "psi-mi:\"MI:0067\"(light scattering)",
               "psi-mi:\"MI:0826\"(x ray scattering)",
               "psi-mi:\"MI:0104\"(static light scattering)",
               "psi-mi:\"MI:0038\"(dynamic light scattering)",
               "psi-mi:\"MI:0968\"(biosensor)",
               "psi-mi:\"MI:0969\"(bio-layer interferometry)",
               "psi-mi:\"MI:0107\"(surface plasmon resonance)",
               "psi-mi:\"MI:0921\"(surface plasmon resonance array)",
               "psi-mi:\"MI:0888\"(small angle neutron scattering)",
               "psi-mi:\"MI:0071\"(molecular sieving)",
               "psi-mi:\"MI:0943\"(detection by mass spectrometry)",
               "psi-mi:\"MI:0069\"(mass spectrometry studies of complexes)")] = "1"
# ? "psi-mi:\"MI:0423\"(in-gel kinase assay)", "psi-mi:\"MI:0097\"(reverse ras recruitment system)"
b1im = as.numeric(b1im)
b1im[is.na(b1im)] = 0

b1imG = b1im #supercategory





start = Sys.time()

# b1edgeldupi = duplicated(b1edgel) | duplicated(b1edgel, fromLast = TRUE) #rows to delete (duplicated)
b1edgeldupi = duplicated(paste(b1edgel[,1],b1edgel[,2])) #rows to delete (duplicated)
b1edgeldup = which(b1edgeldupi) #rows to delete (duplicated)
b1edgeldupnot = which(!b1edgeldupi) #rows not to delete (not duplicated)
# b1edgedupnot = b1edgel[!b1edgeldupi] 

TimeOutput(start)


start = Sys.time()

#for every duplicate, find its orignal row; i.e. dplicate row 11, is duplicate of row 10
b1edgerepi = sapply(b1edgeldup, function(x)
  which(b1edgel[,1]==b1edgel[x,1] & b1edgel[,2]==b1edgel[x,2] & !b1edgeldupi) )
b1edgerepui = unique(b1edgerepi) #unqiue rows with duplicates; calculate scores!

TimeOutput(start) #few min



start = Sys.time()

## calculate scores!
interactionTypeScore = methodScore = pubScore = rep(NA,nrow(b1edgel))

# calculate cores for unique rows with duplications
# for (b1edgerepuii in b1edgerepui) {
loop.ind = loopInd(b1edgerepui,no_cores)
result = foreach (b1edgerepuiiL = loop.ind) %dopar% {
  pubScore = methodScore = interactionTypeScore = NULL
  for (b1edgerepuiii in 1:length(b1edgerepuiiL)) {
    b1edgerepuii = b1edgerepuiiL[b1edgerepuiii]
    edgei = append(b1edgerepuii, b1edgeldup[b1edgerepi==b1edgerepuii]) #rows in edgel that duplicates b1edgerepuii
    at = sum(b1it[edgei])
    bt = at*2 #+ sum(as.numeric(b1itG[edgei]))
    interactionTypeScore[b1edgerepuiii] = log(at+1,bt+1)
    
    am = sum(b1im[edgei])
    bm = am*2 #+ sum(b1imG[edgei])
    methodScore[b1edgerepuiii] = log(am+1,bm+1)
    
    cm = sum(b1np[edgei])
    bp = 7 # number of pubs with max score?
    pubScore[b1edgerepuiii] = log(cm+1,bp+1)
  }
  return(list(pubScore =pubScore, methodScore =methodScore, interactionTypeScore=interactionTypeScore))
}

pubScore1 = unlist(sapply(result, function(x) x$pubScore))
methodScore1 = unlist(sapply(result, function(x) x$methodScore))
interactionTypeScore1 = unlist(sapply(result, function(x) x$interactionTypeScore))

pubScore1[is.na(pubScore1)] = 0
methodScore1[is.na(methodScore1)] = 0
interactionTypeScore1[is.na(interactionTypeScore1)] = 0

pubScore[b1edgerepui] = pubScore1
methodScore[b1edgerepui] = methodScore1
interactionTypeScore[b1edgerepui] = interactionTypeScore1



# calculate scores for unique rows without duplications
edgei0 = setdiff(which(!b1edgeldupi), b1edgerepui) # rows not duplicated
edgei1 = rep(F,length(b1it)); edgei1[edgei0] = T
at0 = b1it[edgei1]
bt0 = at0*2 #+ as.numeric(b1itG[edgei1])
interactionTypeScore2 = log(at0+1,bt0+1)

ab1 = b1im[edgei1]
bb1 = ab1*2 #+ b1imG[edgei1]
methodScore2 = log(ab1+1,bb1+1)

ap0 = b1np[edgei1]
bp0 = 7 # number of pubs with max score?
pubScore2 = log(ap0+1,bp0+1)

pubScore2[is.na(pubScore2)] = 0
methodScore2[is.na(methodScore2)] = 0
interactionTypeScore2[is.na(interactionTypeScore2)] = 0

interactionTypeScore[edgei1] = interactionTypeScore2
methodScore[edgei1] = methodScore2
pubScore[edgei1] = pubScore2

# calculate final scores for all unique rows
finalScore = ((kt*interactionTypeScore) + (km*methodScore) + (kp*pubScore)) / (kt+km+kp)

# check if we've calculated the scores for all unique rows; 
# final scores contains scores for unique rows and NA for duplicate rows
b1edgelunique = is.na(finalScore)
identical(which(b1edgelunique), b1edgeldup)

TimeOutput(start) #23min







## make edge list (multiple edges for each edge with positive weight, DSD doesn't do weighted edges)

finalScoreP = finalScore[!b1edgeldupi]
require(arules)
finalScorePD = as.numeric(discretize(finalScoreP, method="interval", breaks=10)) #discretize edge weights

b1edgelP = b1edgel[!b1edgeldupi,]
rownames(b1edgelP) = 1:nrow(b1edgelP)

b1edgelP_expanded <- b1edgelP[rep(row.names(b1edgelP), finalScorePD),]
write.table(b1edgelP_expanded, file=paste0(edge_dir,"/",data_dir), row.names=F, col.names=F)

system(paste0("python ../DSD-src-0.50/DSDmain.py ", edge_dir,"/",data_dir, " -n 10 -o ", DSD_dir,"/",data_dir))























## Network ------------------------------------------------------------
## plot the network

b0edges = cbind(b1edgelP,finalScoreP)
b0nodes = matrix(unique(as.vector(b0edges)),ncol=1)
g = graph_from_data_frame(b0edges, directed=F)#, vertices=b0nodes)

#which edges goes from and to the same node
b0edgesSame = apply(b0edges, 1, function(x) x[1]==x[2]) 
gplot = graph_from_data_frame(b0edges[!b0edgesSame,], directed=F)

#which nodes aren't connected to any other nodes except itself
b0nodesSame = !b0nodes%in%unique(as.vector(b0edges[!b0edgesSame,]))

#how many of interested genes are amongst those nodes that do have connections to other nodes
# sum(genee%in%b0nodes[!b0nodesSame,])
sum(gene%in%b0nodes[!b0nodesSame,])
# V(gplot)[names(V(gplot))[names(V(gplot))%in%genee]]$color = "orange"
V(gplot)[names(V(gplot))[names(V(gplot))%in%gene]]$color = "red"
x11()
plot(gplot, vertex.size=.1, vertex.label=NA, edge.arrow.size=0)



#breadth first search; find connected components
# bfs()

g <- graph.ring(10)
are.connected(g, 1, 2)
are.connected(g, 2, 4)


