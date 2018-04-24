# using innatedb data
# aya43@sfu.ca 20171017
# Data format is IRefTab: http://irefindex.org/wiki/index.php?title=README_MITAB2.6_for_iRefIndex#What_each_line_represents

root = "~/projects/network/"
result_dir = "result"
setwd(root)



library("stringr")
library("colorspace")
library("lubridate") #if there are date variables
library("data.table")
library("robustHD")
library("foreach")
library("doMC")
library("rols") #ontology look up service
library("GO.db")
library("iRefR") #MITAB format parser
library("rentrez") #NCBI database look up service
library("RISmed")
library("taxize")
library("igraph")
# devtools::install_github("altayg/ganet")
# library(ganet)
# devtools::install_github("R3CPET","sirusb")
# library(R3CPET)

options(na.strings = c("","-", " ", "?", "NULL", NA))
options(stringsAsFactors = F)

source("~/projects/IMPC/code/_funcAlice.R")
source("~/projects/IMPC/code/_funcdist.R")






## Options
control = "[+]_[+]|[+]_Y"
kt = km = kp = 1/3 #paramters for calculating edge weights

## get genes
sampleMeta = get(load("~/projects/IMPC/result/P1/Sanger_SPLEEN/sampleMeta.Rdata")) 
gene0 = unique(sampleMeta$gene)
gene0 = gene0[!grepl(" KO$",gene0)]
gene0 = gene0[!grepl("failed", gene0)]
gene0 = gene0[!grepl(control, gene0)]
gene = str_split(gene0,"_")
gene1 = sapply(gene, function(x) x[1])
gene2 = sapply(gene, function(x) x[2])
identical(gene1,gene2)
sort(unique(gene1[gene2=="+"]))
for (g in unique(gene1[gene2=="+"])) cat("\n(",sum(gene1==g & gene2=="+"),"/",sum(gene1==g),") ", g, sep="")
sort(unique(gene1[gene2=="Y"]))
for (g in unique(gene1[gene2=="Y"])) cat("\n(",sum(gene1==g & gene2=="Y"),"/",sum(gene1==g),") ", g, sep="")
gene = unique(gene1[gene2!="+" & gene2!="Y"]) #delete 118 genes because mostly 1/2 or 1/1 -/+
# gene = gene1
gene = str_split(gene,"[(]") #take out "(b)" at the end of gene name, number of genes dont decrease
gene = sapply(gene, function(x) x[1])
gene = toupper(gene)
# write.table(gene,file=gene_dir,sep="\t",row.names=F,col.names=F,quote=F)




#Setup Cores
no_cores = detectCores()-1
registerDoMC(no_cores)












## Biogrid --------------------------------------------

m0 = data.frame(fread(input="biogrid/BIOGRID-ALL-3.4.154.mitab.txt")) #all
m0 = data.frame(fread(input="biogrid/BIOGRID-ORGANISM-3.4.154.mitab/BIOGRID-ORGANISM-Mus_musculus-3.4.154.mitab.txt")) #mouse
m0 = data.frame(fread(input="biogrid/BIOGRID-MV-Physical-3.4.154.mitab.txt")) # https://wiki.thebiogrid.org/doku.php/biogrid_mv
colnames(m0)

# unique interactors
interactors = unique(append(m0$X.ID.Interactor.A,m0$ID.Interactor.B))
length(unique(append(m0$X.ID.Interactor.A[m0$Taxid.Interactor.A=="taxid:9606" | m0$Taxid.Interactor.A=="taxid:10090"],m0$ID.Interactor.B[m0$ID.Interactor.B=="taxid:9606" | m0$ID.Interactor.B=="taxid:10090"]))) #human | mouse interactors



# "X.ID.Interactor.A"            "ID.Interactor.B"              e.g. entrez gene/locuslink:6416
# "Alt.IDs.Interactor.A"         "Alt.IDs.Interactor.B"         e.g. database:accession;      biogrid:112315|entrez gene/locuslink:MAP2K4
# "Aliases.Interactor.A"         "Aliases.Interactor.B"         e.g. database name:alias;     "entrez gene/locuslink:JNKK(gene name synonym)|entrez gene/locuslink:JNKK1(gene name synonym)|entrez gene/locuslink:MAPKK4(gene name synonym)|entrez gene/locuslink:MEK4(gene name synonym)|entrez gene/locuslink:MKK4(gene name synonym)|entrez gene/locuslink:PRKMK4(gene name synonym)|entrez gene/locuslink:SAPKK-1(gene name synonym)|entrez gene/locuslink:SAPKK1(gene name synonym)|entrez gene/locuslink:SEK1(gene name synonym)|entrez gene/locuslink:SERK1(gene name synonym)|entrez gene/locuslink:SKK1(gene name synonym)"
# "Interaction.Detection.Method"                                e.g. database name:alias;     psi-mi:"MI:0018"(two hybrid); look up controlled vocab term id if not in bracket: https://www.ebi.ac.uk/ols/ontologies/mi (NA may be listed here if aliases are not available); experimental methods (like immunoprecipitations) provide evidence that a list of 3 or more proteins are associated but no evidence for a direct interaction between any given pair of proteins in that list
# "Publication.1st.Author"       "Publication.Identifiers"      e.g. DellAngelica EC (2000), pubmed:10875894
# "Taxid.Interactor.A"           "Taxid.Interactor.B"           e.g. databaseName:identifier; taxid:83333 (taxid since always NCBI's taxonomy database) (NA if it's a complex)
# "Interaction.Types"                                           e.g. database:identifier(interaction type); psi-mi:"MI:0407"(direct interaction). NA may be listed here if the interaction type is not available
# "Source.Database"                                             e.g. database:identifier(source name);      database:identifier(source name)
# "Interaction.Identifiers"                                     e.g. name:identifier;        biogrid:103
# "Confidence.Values"                                           e.g. lpr:1|hpr:12|np:1 (many -; or no score value)
# - lpr score (lowest PMID re-use) is the lowest number of distinct interactions any one PMID (supporting the interaction in this row) is used to support. e.g. 1 = at least one of the supporting PMIDs has never been used to support any other interaction. This likely indicates that only one interaction was described by that reference and that the present interaction is not derived from high throughput methods.
# - hpr score (highest PMID re-use) is the highest ... e.g. high value (e.g. >50) indicates that one PMID describes >50 other interactions and it is more likely that high-throughput methods were used.
# - np score (number PMIDs) is the total number of unique PMIDs used to support the interaction

#Interaction types:
# - Binary interaction data: edgetype=X (Note columns: interactionType, Method)
# - Complexes: n-ary interaction data: edgetype=C. expansion=bipartite. interactor A is placeholder while interactor B list one of the members of the list, therefore, the entire n-ary interaction record is described using one row for each interactor. Each of these rows will have the same interactor A. This method of representation is a bi-partite model since there are two kinds of nodes corresponding to complexes and proteins.
# - Intramolecular interactions and multimers: edgetype=Y / X / C; numParticipants=1 (intra-molecular) / 2 (binary) / 3+ (bipartite; multimer (3 or more) of some protein is being represented; complexes (a.k.a. n-ary data); C - A, edge type 'C', numParticipants=3); interactors = same protein  or  interaction. 
#
#Redundant interactions: group rows with evidence for same/related set of proteins. See columns 33-35 (Checksum_A, Checksum_B and Checksum_Interaction) and 43-51 (integer identifier and canonical data columns): http://irefindex.org/wiki/index.php?title=Canonicalization

#Confidence Value: 
# - https://wiki.thebiogrid.org/doku.php/psi_mitab_file Confidence score. Denoted as scoreType:value. There are many different types of confidence score, but so far no controlled vocabulary. Thus the only current recommendation is to use score types consistently within one source. Multiple scores separated by "|".
# - High throughput screens https://wiki.thebiogrid.org/doku.php/highthroughput: As well as small-scale curated experiments, BioGRID contains data from papers reporting high throughput screens. Since the experimental setup of high throughput screens for interactions can be adjusted to give variable numbers of "hits", it is common for publications to report a confidence score for each interaction. Before loading any high throughput data into BioGRID the data is assessed and interactions considered low confidence by the authors are not loaded. Generally, the confidence score used by the author is loaded, and on BioGRID gene pages containing interactions from that publication the scores can be displayed by choosing the publication from the "Displaying:" drop-down list. In addition, the number of interactions in BioGRID for each publication is displayed next to that publication.



## Delete columns
m0 = delcols(m0)

# Interactor names
y = str_extract(m0$Alt.IDs.Interactor.A,"locuslink:(.*)")
yi = grepl("[|]",y)
y[yi] = gsub("[|]","",str_extract(y[yi],"^.*?[|]"))
m0a = gsub("locuslink[:]","",y) #get the thing between the first "locuslink:" and "| / $"
m0a = toupper(m0a)
names(m0a) = NULL
sum(!is.na(match(gene,m0a)))

y = str_extract(m0$Alt.IDs.Interactor.B,"locuslink:(.*)")
yi = grepl("[|]",y)
y[yi] = gsub("[|]","",str_extract(y[yi],"^.*?[|]"))
m0b = gsub("locuslink[:]","",y) #get the thing between the first "locuslink:" and "| / $"
m0b = toupper(m0b)
names(m0b) = NULL
sum(!is.na(match(gene,m0b)))





## Biogrid: Calculate scores for each unique row/interaction/edge ---------------------------------



# Ineraction type scores
sort(unique(m0$Interaction.Types))
m0it = m0$Interaction.Types
m0it[m0$Interaction.Types%in%c("psi-mi:\"MI:0794\"(synthetic genetic interaction defined by inequality)",
                               "psi-mi:\"MI:0796\"(suppressive genetic interaction defined by inequality)",
                               "psi-mi:\"MI:0799\"(additive genetic interaction defined by inequality)")] = "0.1"
m0it[m0$Interaction.Types%in%c("psi-mi:\"MI:0403\"(colocalization)")] = "0.33"
m0it[m0$Interaction.Types%in%c("psi-mi:\"MI:0407\"(direct interaction)")] = "1"
m0it[m0$Interaction.Types%in%c("psi-mi:\"MI:0914\"(association)")] = "0.33"
m0it[m0$Interaction.Types%in%c("psi-mi:\"MI:0915\"(physical association)")] = "0.66"
m0it = as.numeric(m0it)

m0itG = m0it #supercategory
m0itG[m0$Interaction.Types%in%c("psi-mi:\"MI:0914\"(association)","psi-mi:\"MI:0407\"(direct interaction)","psi-mi:\"MI:0915\"(physical association)")] = "0.33"


# Interaction detection method scores
sort(unique(m0$Interaction.Detection.Method))
m0im = m0$Interaction.Detection.Method
m0im[m0$Interaction.Detection.Method%in%c("psi-mi:\"MI:0401\"(biochemical)",
                                          "psi-mi:\"MI:0047\"(far western blotting)",
                                          "psi-mi:\"MI:0415\"(enzymatic study)",
                                          "psi-mi:\"MI:0004\"(affinity chromatography technology)",
                                          "psi-mi:\"MI:0096\"(pull down)",
                                          "psi-mi:\"MI:1313\"(bioid)")] = "1"
m0im[m0$Interaction.Detection.Method%in%c("psi-mi:\"MI:0055\"(fluorescent resonance energy transfer)", #biophysical
                                          "psi-mi:\"MI:0114\"(x-ray crystallography)")] = "1"
m0im[m0$Interaction.Detection.Method%in%c("psi-mi:\"MI:0090\"(protein complementation assay)",
                                          "psi-mi:\"MI:0018\"(two hybrid)")] = "0.66"
m0im[m0$Interaction.Detection.Method%in%c("psi-mi:\"MI:0428\"(imaging technique)")] = "0.33"
m0im[m0$Interaction.Detection.Method%in%c("psi-mi:\"MI:0254\"(genetic interference)")] = "0.1"
m0im[m0$Interaction.Detection.Method%in%c("psi-mi:\"MI:0686\"(unspecified method)")] = "0.05"
m0im = as.numeric(m0im)

m0imG = m0im #supercategory





start = Sys.time()

m0edgel = cbind(m0a,m0b)
m0abi = m0edgel[,1]>m0edgel[,2]
m0edgel[m0abi,] = cbind(m0edgel[m0abi,2],m0edgel[m0abi,1])

# m0edgeldupi = duplicated(m0edgel) | duplicated(m0edgel, fromLast = TRUE) #rows to delete (duplicated)
m0edgeldupi = duplicated(paste(m0edgel[,1],m0edgel[,2])) #rows to delete (duplicated)
m0edgeldup = which(m0edgeldupi) #rows to delete (duplicated)
m0edgeldupnot = which(!m0edgeldupi) #rows not to delete (not duplicated)
m0edgeu = m0edgel[!m0edgeldupi] 

TimeOutput(start)

start = Sys.time()

m0edgerepi = sapply(m0edgeldup, function(x) which(m0edgel[,1]==m0edgel[x,1] & m0edgel[,2]==m0edgel[x,2] & !m0edgeldupi) ) #for each m0edgeldup, row in edgel (& m0edgeldupnot) it duplicates
m0edgerepui = unique(m0edgerepi)

interactionTypeScore = methodScore = pubScore = rep(NA,nrow(m0edgel))

TimeOutput(start) #23min

start = Sys.time()

# calculate initial scores for rows with duplications
# for (m0edgerepuii in m0edgerepui) {
loop.ind = loopInd(m0edgerepui,no_cores)
result = foreach (m0edgerepuiiL = loop.ind, .combine="c") %dopar% {
  pubScore0 = NULL
  for (m0edgerepuiii in 1:length(m0edgerepuiiL)) {
    m0edgerepuii = m0edgerepuiiL[m0edgerepuiii]
    edgei = m0edgeldup[m0edgerepi==m0edgerepuii] #rows in edgel that duplicates m0edgerepuii
    at = sum(m0it[edgei])
    bt = at + sum(as.numeric(m0itG[edgei]))
    interactionTypeScore[m0edgerepuii] = log(at+1,bt+1)
    
    am = sum(m0im[edgei])
    bm = am + sum(m0imG[edgei])
    methodScore[m0edgerepuii] = log(am+1,bm+1)
    
    bp = 7 # number of pubs with max score?
    pubScore0[m0edgerepuiii] = log(length(edgei)+1,bp+1)
  }
  return(pubScore0)
}
pubScore[m0edgerepui] = result


# calculate initial scores for rows without duplications
edgei0 = intersect(which(!m0edgeldupi), m0edgerepui) # rows not duplicated
edgei1 = rep(F,length(m0it)); edgei1[edgei0] = T
at0 = m0it[edgei1]
bt0 = at0 + as.numeric(m0itG[edgei1])
interactionTypeScore[edgei1] = log(at0+1,bt0+1)

am0 = m0im[edgei1]
bm0 = am0 + m0imG[edgei1]
methodScore[edgei1] = log(am0+1,bm0+1)

bp0 = 7 # number of pubs with max score?
pubScore[edgei1] = log(1+1,bp0+1)


# calculate final scores for all unique rows
finalScore = ((kt*interactionTypeScore) + (km*methodScore) + (kp*pubScore)) / (kt+km+kp)

# check if we've calculated the scores for all unique rows; final scores contains scores for unique rows and NA for duplicate rows
m0edgelunique = is.na(finalScore)
identical(which(m0edgelunique), m0edgeldup)







## Biogrid: Convert to graph (unique edges) ----------------------------------------
# m0edges = m0edgel
m0edges = cbind(m0edgel[!m0edgelunique,],finalScore[m0edgelunique])
m0nodes = matrix(unique(as.vector(m0edges)),ncol=1)
g = graph_from_data_frame(m0edges, directed=F)#, vertices=m0nodes)

#which edges goes from and to the same node
m0edgesSame = apply(m0edges, 1, function(x) x[1]==x[2]) 
gplot = graph_from_data_frame(m0edges[!m0edgesSame,], directed=F)

#which nodes aren't connected to any other nodes except itself
m0nodesSame = !m0nodes%in%unique(as.vector(m0edges[!m0edgesSame,]))

#how many of interested genes are amongst those nodes that do have connections to other nodes
sum(genee%in%m0nodes[!m0nodesSame,])
sum(gene%in%m0nodes[!m0nodesSame,])
V(gplot)[names(V(gplot))[names(V(gplot))%in%genee]]$color = "orange"
V(gplot)[names(V(gplot))[names(V(gplot))%in%gene]]$color = "red"
plot(gplot, vertex.size=.1, vertex.label=NA, edge.arrow.size=0)








## Biogrid: Human/Mice interactors ------------------------------------


taxid = "taxid:9606" #human
taxid = "taxid:10090" #mice

genei = ( sort(intersect( unique(m0a[m0a%in%gene & m0$Taxid.Interactor.A!=taxid]), unique(m0a[m0a%in%gene & m0$Taxid.Interactor.A==taxid]) )) )  #genes belonging to both mice and other taxonomy species
geneinonmice = table(m0a[m0a%in%genei & m0$Taxid.Interactor.A!=taxid])
geneimice = table(m0a[m0a%in%genei & m0$Taxid.Interactor.A==taxid])
paste0(genei,": ", geneimice, "/", geneinonmice)[order(geneimice/geneinonmice)] #ratio of the gene in mice against non-mice

setdiff(unique(m0a[m0a%in%gene & m0$Taxid.Interactor.A!=taxid]), unique(m0a[m0a%in%gene & m0$Taxid.Interactor.A==taxid])) #genes belonging to non-mice only
# without removing 1/2 genes: "STXBP4" "RNF157" "CLPP"   "TUBA3A" "BRD" 
# with: "SUPT5"  "BARHL1" "STXBP4" "PIGF"   "RNF157" "CLPP"   "TUBA3A" "BRD"  




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















## Innatedb ---------------------------------------------------

# All Imported Experimentally Validated Interactions (get rows for mouse only)
b0 = data.frame(fread(input="innatedb/all.mitab"))
colnames(b0)

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



## Delete columns
b0 = delcols(b0)

# What database names can be used to look up in ontology
b0a = gsub("[:]|[(]display_short[)]","",str_extract(b0$alias_A,":(?!.*:).*?[(]display_short[)]")) #get the thing between ":" and "(display_short)"
b0a = toupper(b0a)
names(b0a) = NULL
sum(!is.na(match(gene,b0a)))

b0b = gsub("[:]|[(]display_short[)]","",str_extract(b0$alias_B,":(?!.*:).*?[(]display_short[)]")) #get the thing between ":" and "(display_short)"
b0b = toupper(b0b)
names(b0b) = NULL
sum(!is.na(match(gene,b0b)))





b0edgel = cbind(b0a,b0b)
sum(duplicated(b0edgel) | duplicated(b0edgel, fromLast = TRUE)) #see if there are any duplicate edges






#see data regarding taxonomy
sum(b0$ncbi_taxid_A!=b0$ncbi_taxid_B)
sum(b0$ncbi_taxid_A!=b0$ncbi_taxid_B & ((!grepl("Human",b0$ncbi_taxid_A) & !grepl("Mouse",b0$ncbi_taxid_A)) | (!grepl("Mouse",b0$ncbi_taxid_A) & !grepl("Human",b0$ncbi_taxid_A)))) # where two interactors have different taxonomy and they are not mouse or human





bmouse = grepl("Mouse",(b0[,"ncbi_taxid_A"]),ignore.case=T) | grepl("Mouse",(b0[,"ncbi_taxid_B"]),ignore.case=T)
b = b0[bmouse,]
cat("\n# deleted rows: ", sum(!bmouse), "/", nrow(b0), sep="")











## InnateDB Curated Interactions (get rows for mouse only)
a0 = data.frame(fread(input="innatedb/innatedb_all.mitab"))

# Delete columns
a0 = delcols(a0)

# What database names can be used to look up in ontology
a0a = gsub("[:]|[(]display_short[)]","",str_extract(a0$alias_A,":(?!.*:).*?[(]display_short[)]")) #get the thing between ":" and "(display_short)"
a0a = toupper(a0a)
names(a0a) = NULL
sum(!is.na(match(gene,a0a)))

a0b = gsub("[:]|[(]display_short[)]","",str_extract(a0$alias_B,":(?!.*:).*?[(]display_short[)]")) #get the thing between ":" and "(display_short)"
a0b = toupper(a0b)
names(a0b) = NULL
sum(!is.na(match(gene,a0b)))





a0edgel = cbind(a0a,a0b)
sum(duplicated(a0edgel) | duplicated(a0edgel, fromLast = TRUE)) #see if there are any duplicate edges








#see data regarding taxonomy
sum(a0$ncbi_taxid_A!=a0$ncbi_taxid_B)
sum(a0$ncbi_taxid_A!=a0$ncbi_taxid_B & ((!grepl("Human",a0$ncbi_taxid_A) & !grepl("Mouse",a0$ncbi_taxid_A)) | (!grepl("Mouse",a0$ncbi_taxid_A) & !grepl("Human",a0$ncbi_taxid_A)))) # where two interactors have different taxonomy and they are not mouse or human



amouse = grepl("Mouse",(a0[,"ncbi_taxid_A"]),ignore.case=T) | grepl("Mouse",(a0[,"ncbi_taxid_B"]),ignore.case=T)
a = a0[amouse,]
cat("\n# deleted rows: ", sum(!amouse), "/", nrow(a0), sep="")













# # merge two interaction matrices
# ab0 = rbind(a,b)
# abcolunique0 = sapply(1:ncol(ab0), function(x) length(unique(ab0[,x])))
# ab = ab0[,abcolunique0>1] #get rid of columns with only one value
# cat("\ndeleted columns: ",paste(colnames(ab0)[!abcolunique0>1],collapse=", "), sep="")
# abcolunique = sapply(1:ncol(ab), function(x) length(unique(ab[,x])))
# names(abcolunique) = colnames(ab)
# ab = ab[,!colnames(ab)%in%c("alt_identifier_A","alt_identifier_B","expansion_method")]
# ab1 = ab[!duplicated(lapply(ab, digest))]














## Ontologies --------------------------------------------------------------
qry = olsSearch(OlsSearch(q = gene[1]))

entrez_dbs() #databases in NCBI















## Network ------------------------------------------------------------

#breadth first search & are connected to find connected components
bfs()

g <- graph.ring(10)
are.connected(g, 1, 2)
are.connected(g, 2, 4)


