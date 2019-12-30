# flow Feature Biclustering Pipeline
A pipeline for extracting features from and biclustering flow cytometry data


- Input: FlowType files derived from Flow Cytometry (FCM) FCS files
- Output: Biclusters (file clusters and cell population clusters)

## Input Data
The original FlowCAP-II AML Flow Cytometry FCS files can be found here: https://flowrepository.org/id/FR-FCM-ZZYA

The IMPC data set has yet to be finalized for public viewing

Each FCS file need to be first converted to a FlowType Rdata file. This can be done using the 'FlowType' function, on parameter: method='Thresholds'. Find a tutorial of the FlowType package here: https://www.bioconductor.org/packages/devel/bioc/html/flowType.html

The resulting input matrix should be a **(FCS file x flowtype cell population) matrix** i.e. the output from the FlowType procedure is a vector of cell count values for each cell population; the matrix would contain a vector for each file


## Pre-requisites

- Everything is written in R. All the libraries needed are listed at the top of every script -- please install as needed.
- Paths and Script options are listed at the top of every script -- please change as needed.

## Running the Code
The code is numbered based on dependency (e.g. 06 is depentant on 05, 04, etc.); please replace the directories as needed

### FlowCAP-II data set
To run the pipeline for the FlowCAP-II data set, run the scripts in the following order:
-	[01a_normalize.R](flowCAP-II/01a_normalize.R) **normalizes** the cell counts
- [02_pvalue_single.R](flowCAP-II/02_pvalue_single.R) creates the **p value** feature matrix
-	[02_plot.R](flowCAP-II/02_plot.R) plots statistics on cell counts in different files
-	[02_stats_count.R](flowCAP-II/02_stats_count.R) calculates statistics on cell counts in different files
-	[03_childparent_matrix.R](flowCAP-II/03_childparent_matrix.R) creates the rest of the FCM file **features**
-	[04a_biclust.R](flowCAP-II/04a_biclust.R) **biclusters** feature matrices
- [04b_biclust_score.R](flowCAP-II/04b_biclust_score.R) **scores** those clusterings
-	[04b_biclust_plot.R](flowCAP-II/04b_biclust_plot.R) creates plots for clusterings

### IMPC data set
The IMPC data set has several confounding factors, including time file was created on, therefore to account for confounding factors, we run additional scripts:
- [01c_confound_peer-.R](IMPC/01c_confound_peer-.R) adjusts normalized count feature for confounding factors.
- [02a_timegroup.R](IMPC/02a_timegroup.R) splits files by the time covariate/confounding factor
- [02b_pvalue_singlefile.R](IMPC/02b_pvalue_singlefile.R) creates the **p value** feature matrix
- [02.5_pvalue_count_distribution_plots.R](IMPC/02.5_pvalue_count_distribution_plots.R) plots distribution of p value feature value

IMPC also requires an additional scoring metric based off of the gene-protein interaction  network; in order to use this metric, run the following:
- [00_data_innatedb.R](00_data_innatedb.R) calculates distances between genes (representative of the IMPC files as each file is made from a mouse with a different gene knocked out) needed to score clusterings; we download the TAB networks from [innatedb](http://www.innatedb.ca/redirect.do?go=downloadCurated); optionally [biogrid](https://downloads.thebiogrid.org/BioGRID) (accessed 201804)
- [04b_biclust_score.R](IMPC/04b_biclust_score.R) **scores** those clusterings


