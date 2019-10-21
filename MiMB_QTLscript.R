#3. Methods
#3.1. Preparing and Loading Data
#1.Prepare the data as a .csv file with the correct headings (see Notes 1 and 2). Genotypes are standardly in the format “AA,” “AB,” and “BB.” All missing phenotypes and genotypes should be noted by “NA” to work with commands below (see Note 3).

#2.Once R (or R studio) is open, load the r/QTL software.
library(qtl)

#3.From the toolbar menu, select Session, Set Working Directory, then Choose Directory to navigate to the folder that contains your .csv data file.

#4.Read in the genotypes, genetic map (see Note 4), and phenotypes for your mapping cross from your .csv file. In the below example, you will be prompted to select your file.
data<-read.cross(format="csv", file=file.choose(), genotypes=c("AA", "AB", "BB"), na.strings=c("NA"), convertXdata=TRUE)



#3.2. Data Verification
#1.High-quality data is critical and errors in genotype or phenotype data can produce strange mapping results. As a first step, look at details of the data set and verify the type of cross, number of individuals, number of phenotypes (see Note 5), and genotypic data (see Note 6).
summary(data)

#2.Visually inspect your data including missing genotypes (see Subheading 3.2, step 4), the genetic map, and histogram distributions for each phenotype column.
plot(data)

#3.Visualize your mapping data (see Note 7).
geno.image(data)

#4.Identify markers and hybrids that have too many missing genotypes (see Note 8).
plotMissing(data)

#5.Remove individuals that are missing too many genotypes. The below will only retain in the data set individuals that have greater than 50 genotypes in the data set.
data<-subset(data, ind=(ntyped(data)>50))

#6.Remove markers that have too much missing data. The ntyped command will list the number of individuals genotyped for each marker. Any markers below a certain threshold can be removed by the drop.markers command and listing the name of the specific marker to drop, here marker1.
ntyped(data, "mar") 
data<-drop.markers(data, "marker1")

#7.Fill in your missing data (i.e., “NA” genotypes) (see Note 9).
augdata<-mqmaugment(data,minprob=0.1)



#3.3. Building a Statistical Model
#1.Conduct an initial scan, here for the phenotypic data in column 2 (see Note 10) to identify putative unlinked QTL that will be used to build a more rigorous model (see Note 11).
scan<-mqmscan(augdata,pheno.col=2,plot=T,model="dominance",verbose=FALSE)

#Running the below will generate the same plot without gridlines
plot(scan)

#2.View a summary of marker locations with the top LOD scores per chromosome (see Note 12) to identify putative QTL that you will use as cofactors to generate a more robust model (see Note 13).
summary(scan)

#3.For each putative QTL (see Notes 14 and 15), visually verify the effect of genotype on phenotype by generating an effect plot. The below command is for the locus on chromosome 1 at 5 centiMorgans (cM) (chr=1, pos=5). If at least two genotype groups do not have overlapping ranges of phenotypes, this is a putative QTL that willshould be added as a cofactor in your model. More details on effect plots are discussed in Subheading 3.5, step 3.
effect<-effectplot(augdata,pheno.col=2,mname1=find.marker(augdata, chr=1, pos=5))

#4.Establish the file of cofactors you will use.
cofactorslist<-NULL

#5.Build your list of cofactors (see Note 16), repeating the below with new chromosome and position values for all putative QTL (see Notes 17 and 18).
cofactorslist<-c(cofactorslist,find.markerindex(augdata,find.marker(augdata,chr=1,pos=5)))

#6.Generate the cofactor matrix (see Notes 19–21).
cofactors<-mqmsetcofactors(augdata,cofactors=c(cofactorslist))

#7.Run a QTL scan using your updated model with cofactors. You will use the output of this to further refine your model (see Note 22).
scan<- mqmscan(cross=augdata,pheno.col=2,cofactors=cofactors,cofactor.significance=0.002, verbose=T,plot=T,model="dominance")
summary (scan)

#8.If there are any new putative QTL, assess the effect plot at this locus (see Subheading 3.3, step 3), adding these as cofactors in the model as necessary (see Subheadings 3.3, steps 5 and 6).
#9.Repeat the process until the model stabilizes (see Notes 23–25): build a new model with additional cofactors, run the scan with the updated model, assess new putative QTL using effect plots, and build a new model with these updated cofactors.



#3.4. Statistical Significance
#1.Determine the threshold of LOD that meet the 5% and 10% statistical significance levels (i.e., p < 0.05 and p < 0.10, respectively, see Note 26) by running a permutation on your final model, here 1000 permutations (see Notes 27 and 28).
result<-mqmpermutation(cross=augdata,scanfunction=mqmscan,cofactors=cofactors,pheno.col=2,n.perm=1000,plot=F,verbose=T,model="dominance")
resultqtl<-mqmprocesspermutation(result)
summary(resultqtl)

#2.For each significant QTL, run Bayesian analysis for each peak to get 95% confidence interval and closest flanking markers (see Note 29). Replace the “X” in the below with your chromosome number.
bayesint(scan, X, qtl.index=1,prob=0.95,lodcolumn=1, expandtomarkers=T)



#3.5. Assessing QTL
#1.View the LOD plot across the genome.
plot(scan)

#2.View the LOD plot across a single chromosome (see Notes 30 and 31).
plot(scan,chr=1)

#3.Analyze the effect plots at the peak LOD of your QTL to assess mode of inheritance (see Fig. 1, Notes 32 and 33). The effect command will output the phenotypic mean and standard error values (see Note 34). The below command is for the locus on chromosome 1 at 5 cM.
effect<-effectplot(augdata,pheno.col=2,mname1= find.marker(augdata,chr=1,pos=5))
effect

#4.Calculate the additive, dominance effects, and heritability for the peak marker for each QTL locus (see Note 34), replacing AA, AB, and BB in the below formulas with the mean values outputted from the effect command.
effect
add<-((AA-BB)/2)
add
dom<-(AB -(AA+BB)/2)
dom
heritability<-(2*add^2+dom^2)/(2*add^2+dom^2+4*1)
heritability

#5.Calculate the phenotypic variance explained (PVE) by the QTL (see Note 35), replacing n in the below formula with the number of individuals included in the analysis (see Note 9) and LOD in the below with the outputted LOD score for that locus.
PVE<-(1-(10^-((2/n)*(LOD))))
PVE



#3.6. From QTL to Candidate Genes and Causative Loci
#1.QTL loci are commonly large genetic intervals (see Note 31), and considerable additional work is necessary to determine candidate genes and causative genetic changes. Often, candidate genes within a QTL interval are chosen based on bioinformatic analysis, previous literature searches, and other genetic information such as population genetics. Defining causal mutations is particularly difficult and requires a great deal of follow-up work including experiments such as gene expression analysis, embryonic manipulations, characterization of loss-of-function phenotypes, and gene editing (e.g., using the CRISPR-Cas9 system).



#4. Notes
#1.Data should be in .csv format. This format can easily be produced from a spreadsheet in OpenOffice or Microsoft Excel. The first column should be an individual ID, followed by columns as needed with data on one phenotype each. The remaining columns are genotype data, one per marker. The first row is a header, with a description of the phenotype or the marker name. The second and third rows are reserved for genetic map information, or are left blank if a genetic map is being estimated (see Note 4). This second and third row contains chromosome number as a numeral and the third row contains  centiMmorgan (cM) position, respectively, for all markers; the second and third rows must be empty for the ID and phenotype columns.
#2.A .csv formatted mapping file in the appropriate format is supplied at https://github.com/kpowder/MiMB_QTL. This file is highly modified partial data set from [41]. When including all individuals and markers as given and building a model with markers 5 and 119 as cofactors, a significant QTL is present on chr 1. This QTL has an additive effect, LOD = 3.85, peak at marker 5, and a 95% CI from marker 1 to marker 22 (0–18.2 cM).
#3.Additional genotype formats are acceptable, for example, “A,” “B,” and “C” for a microsatellite with three different lengths. If using a different set of strings for genotypes or missing data, this can be adjusted in the command in Subheading 3.1, step 4.
#4.If the genetic map has not been assembled yet, it can be estimated at this time by including estimate.map=TRUE in the read.cross command, such as below. Be warned that including a map estimation will make this command take a while to process.
data<-read.cross(format="csv", file=file.choose(), genotypes=c("AA", "AB", "BB"), na.strings=c("NA"), estimate.map=TRUE, convertXdata=TRUE)

#5.The number of phenotypes listed in by the summary(data) command will include any non-genotype columns. For instance, in the sample data file, there are three “phenotypes” of ID, the phenotype data, and sex.
#6.Full genotype data can be attained below, including genotype frequencies at each marker and a p value for tests of 1:2:1 Mendelian proportions expected in an intercross.
geno.table(data)

#7.In this visual, AA genotypes are represented by red, BB by green, AB by blue, and NA by white colors. Rows are individuals and columns are markers.
#8.In this visual, black marks indicate missing genotypic information. Rows are individuals and columns are markers.
#9.Missing genotypes would need to be excluded from the statistical model, reducing power. Using the mqmaugment function statistically predicts missing genotypes based on neighboring markers and recombination rates [11]. Genotypes are modeled and all possible genotypes are created as new individuals with weighted probabilities. When visualizing the augmented data with the below, you will note that there are no longer white (missing) genotypes and the number of individuals increases based on the weighted genotypes. Note that any time you need the number of individuals (e.g., calculating PVE in Subheading 3.5, step 5), you should use the original number of individuals included rather than the number of individuals after augmentation.
geno.image(augdata)

#10.The program counts column A as column 1 when counting columns.
#11.I generally use model = "dominant" as this will assess alleles with both dominant and additive effects. The default option model = "additive" will only assess for additive effects.
#12.This summary command will list only the top LOD scores across the genome. Additional useful commands are below.
scan #displays LOD scores at every marker 
write.csv(scan,file="scan_LODs.csv") #outputs all LOD scores as a csv file.

#13.Use of cofactors in the model eliminates phenotypic variation due to other QTL and more accurately estimates the effect of a QTL.
#14.As a rule of thumb, any locus with a LOD >2.5–3 is a putative QTL. Once the full model is built, a LOD threshold will be empirically determined using permutations (see Subheading 3.4, step 1).
#15.If there are no obvious LOD “peaks” of >2.5, you may want to run an autoscan. The below commands will randomly pick 250 loci across the genome, accounting for marker density, that can be used as an initial set of cofactors to build the model.
cofactors <- mqmautocofactors(augdata,250)

#16.In order to treat sex as a cofactor, include the sex-determining locus as a cofactor. If other fixed cofactors such as family structure, environment, or diet are necessary to include, you may need to build models and conduct backward elimination manually using the scanone function (see [9]).
#17.Running >cofactorslist once you have entered all your putative QTL will give you a list of marker numbers to ensure you’ve got them all.
cofactorslist 

#18.While loci are commonly outputted as a chromosome number and chromosomal position, the software also assigns each marker a number. I find it helpful to keep a scratch piece of paper or file with both of these pieces of information. The below will display the marker number based on chromosomal location (here, at 5 centiMorgans on chromosome 1).
find.marker(augdata,chr=1,pos=5)

#19.Running >cofactors after this will give you a matrix of all markers in your data set. Those included as a cofactor will be a “1” and those not included are “0”. This matrix can be overwhelming, but a count of “1” entries is quick check that you have the expected number of cofactors included.
cofactors

#20.If you know the marker number (see Note 18), you can quickly adjust this matrix without adjusting the list of cofactors (i.e., not going back to Subheadings 3.3, step 4 or 5). In the example below, marker 34 is switched to be included, while marker 78 is switched to be excluded.
cofactors[34]<-1

cofactors[78]<-0

#21.Use the below to visualize how your cofactors are spread across your chromosomes.
mqmplot.cofactors(augdata, cofactors, justdots=T)

#22.This analysis conducts an unsupervised backward elimination to verify cofactors, in whichmeaning a cofactor is removed and the analysis is recalculated. I generally remove the eliminated cofactors from the next iteration of the model as too many cofactors decreases power [82]. Note that cofactors eliminated from one model may not be eliminated from all models.
#23.Again, a piece of scratch paper can come in handy as you run through different iterations of the models. I generally record the locus chromosome and position information and marker number of cofactors in the model, as well as which of these were eliminated.
#24.I generally consider the model “stable” when (1) all putative QTL peaks are included as cofactors and retained in the model or (2) any putative QTL peak(s) not included in the model is eliminated when added.
#25.The total number of cofactors in the model should not exceed two times the square root of the number of individuals [94].
#26.Generally, loci that pass the 5% threshold are significant QTL, while those that are above the 10% are considered suggestive.
#27.Bonferroni correction is not appropriate to adjust for multiple testing given that markers are potentially linked. Rather, QTL mapping uses randomized permutations to empirically establish a significance threshold [95]. This permutation also accounts for variation in the data set such as number of individuals, number of markers, pattern of missing data, and the phenotypic variation. Most commonly, 1000 permutations are conducted.
#28.Remember that a LOD of >2.5–3 was used as a rule of thumb for putative QTL when building the model. If the permuted significance level is less than this value, you should look at the effect plots for the QTL that now would be considered significant and adjust your model if needed. If you adjust the model changes, you will need to re-rerun the permutation.
#29.Another option instead of calculating peak interval with Bayesian analysis is to determine where there is a drop of 1.5 in the LOD score from the peak. However, the 1.5-lod support interval can vary greatly [96]. Replace the “X” in the below with your chromosome number.
lodint(scan,X,qtl.index=1,drop=1.5,lodcolumn=1, expandtomarkers=T)

#30.Plots of LOD score can be used to help clarify if there may be two closely linked QTL. If there is a dip in the LOD score, this can be suggestive of multiple peaks in the same region rather than a single, large peak, though additional work is necessary to confirm this. Looking at effect plots across the region may also help distinguish if effects are consistent across the region or have variation that may suggest multiple peaks.
#31.The level of resolution for a QTL analysis is based on density of markers and the amount of recombination; for many F2 QTL analyses the limit of resolution is in the range of 3–20 cM [1097]. Thus, it may not be possible to distinguish between a single QTL locus and two QTL that are tightly linked. Finer mapping of the QTL interval can be conducted in a later generation after additional rounds of inbreeding (advanced intercross lines).
#32.Two-way effect plots can be used to assess nonadditive effects between two specific genetic loci (also called gene interactions or epistasis) [9899100]. The below command is for the loci on chromosome 1 at 5 cM and chromosome 3 at 25 cM.
effect2<-effectplot(data,pheno.col=2,mname1= find.marker(augdata, chr=1, pos=5), mname2= find.marker(augdata, chr=3, pos=25))
effect2

#33.It is possible to conduct genome-wide scans for epistatic interactions, in which the model is fitted to every possible pairwise combination of loci [99899]. This has a variety of challenges including increased sample sizes due to partitioning of genotype groups, corrections for multiple testing, and increased computational demands that often require parallel computing [9899].
#34.It is common practice to report these additive, dominant, and heritability values (see Subheading 3.5, step 4) for the peak marker of each QTL in tables in publications, as well as mean phenotypic values from effect plot (see Subheading 3.5, step 3). Additionally, the LOD score (see Subheading 3.3, step 2 and Note 12), confidence interval (see Subheading 3.4, step 2), and PVE (see Subheading 3.5, step 5) should be reported in publications.
#35.Small sample sizes and thus decreased power can lead to an overestimation of the phenotypic variance explained (PVE), called the Beavis effect [101102].