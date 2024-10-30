#sequences were first joined and demultiplexed by using Qiime scripts
#no quality filtering was applied at this time
#below are the commands performed in Qiime to join ends and demultiplex the samples
#join_paired_ends.py -f /Users/luizroesch/uclust8/f_seqs.fastq -r /Users/luizroesch/uclust8/r_seqs.fastq -b /Users/luizroesch/uclust8/barcodes.fastq -o joined
#split_libraries_fastq.py -i /Users/luizroesch/uclust8/Analysis/joined/fastqjoin.join.fastq -b /Users/luizroesch/uclust8/Analysis/joined/fastqjoin.join_barcodes.fastq   -o splitout -m /Users/luizroesch/uclust8/map.txt --store_demultiplexed_fastq -r 0 -q 0 -n 100 --barcode_type 11 --rev_comp_barcode
#split_sequence_file_on_sample_ids.py -i /Users/luizroesch/Desktop/Saliva/seqs/Analysis/splitout/seqs.fastq --file_type fastq -o out_fastq/

#cheeck the primers
#https://benjjneb.github.io/dada2/ITS_workflow.html#identify-primers


setwd("~/Downloads/Run19/")
library(dada2); packageVersion("dada2")
#directory containing the fastq files after unzipping
path = "combined/"
list.files(path)

#read in the names of the fastq files, 
#and perform some string manipulation to get lists of 
#the fastq files in matched order
fnFs <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), ".fastq.gz"), `[`, 1)
sample.names

#visualizing the quality profiles of the reads for the first 4 samples
plotQualityProfile(fnFs[1:5])

#Assign the filenames for the filtered fastq.gz files
#OPTIONAL: create 2 filepaths, as below, to compare results of different quality cutoffs
#filtFslength <- file.path("merged/", "filtered_length", paste0(sample.names, ".fastq"))
filtFslenient <- file.path(".", "filtered_lenient", paste0(sample.names, "_filt.fastq.gz"))
filtFsstringent <- file.path(".", "filtered_stringent", paste0(sample.names, "_filt.fastq.gz"))

#Filter and trim

#Filter just for length (toss reads > 300 bp)
#outlength <- filterAndTrim(fnFs, filtFslength, maxLen = 300, compress=FALSE, multithread=TRUE)
#For a more lenient quality cutoff, run the command below
outlenient <- filterAndTrim(fnFs, filtFslenient, trimLeft = 29, truncLen=c(225), maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
#For a more stringent quality cutoff, run the command below
outstringent <- filterAndTrim(fnFs, filtFsstringent, trimLeft = 29, truncLen=c(225), maxN=0,truncQ=11,maxEE=c(2), rm.phix=TRUE, compress=TRUE, multithread=TRUE)
#Look at the results for both quality methods and make a decision
head(outlenient)
head(outstringent)

#DECISION: for these samples we are going to use the more lenient cutoff

filtpath <- "filtered_lenient/"
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE)
sample.names <- sapply(strsplit(basename(filts), "_filt"), `[`, 1)
names(filts) <- sample.names

# Learn error rates
set.seed(2125)
err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE, MAX_CONSIST = 20)

#Look at error rate estimates and plotted values to check for good model fitting
#If the red line is not a "good fit" for the data points, increase "nbases =" in the error estimation above
plotErrors(err, nominalQ=TRUE)

# Infer sequence variants
set.seed(2125)
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}

# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "seqtab.rds")

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
dim(seqtab)

tax <- assignTaxonomy(seqtab, "~/Downloads/new_combined_data_1_4_2022/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
#inspect the taxonomic assignments
taxa.print <- tax
rownames(taxa.print) <- NULL
head(taxa.print)
table(nchar(getSequences(seqtab)))

#make species level assignments
taxa <- addSpecies(tax, "~/Downloads/new_combined_data_1_4_2022/silva_species_assignment_v138.1.fa")

#inspect the taxonomic assignments after adding species
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
table(nchar(getSequences(seqtab)))


# Write to disk
saveRDS(seqtab, "seqtab.rds") 
saveRDS(taxa, "taxa.rds")

# Make csv files
write.csv(seqtab, "seqtab.csv")
write.csv(taxa, "taxa.csv")

seqtab1=as.data.frame(seqtab)
library(tidyverse)
seqtab1=seqtab1%>%
  add_column(Count=rowSums(seqtab1),.before = 1)
hist(seqtab1$Count)
#import final rds file
#seqtab = readRDS("seqtab.rds")
#taxa = readRDS("taxa.rds")

# Construct phylogenetic tree ---------------------------------------------

#NOTE: if having trouble installing DECIPHER due to error with RcppArmadillo, run the following
#under /Applications/Utilities
#sudo curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
#sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /

#Start by making a multiple sequence alignment on the seqtab output using "DECIPHER"
#library(DECIPHER)
#seqs <- getSequences(seqtab)
#names(seqs) <- seqs # This propagates to the tip labels of the tree
#alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, processors = 3)

#Construct the phylogeny using "phangorn"
#library(phangorn)
#phang.align <- phyDat(as(alignment, "matrix"), type="DNA") #import alignment
#dm <- dist.ml(phang.align) #create distance matrix
#treeNJ <- NJ(dm) # construct tree using Neighbor Joing method # Note, tip order != sequence order
#fit = pml(treeNJ, data=phang.align) #compute tree likelihood

#fitGTR <- update(fit, k=4, inv=0.2) #update and refit substitution model
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                    rearrangement = "stochastic", control = pml.control(trace = 0)) #optimize likelihood model
#detach("package:phangorn", unload=TRUE)

# OR use Raxml as suggested in the following issue (uses the ips package)
#https://github.com/benjjneb/dada2/issues/88
#library(ips)
#alignment <- read.dna("alignment.fasta",format="fasta",as.matrix=TRUE)
#alignment.rax.gtr <- raxml(alignment,
#                           m="GTRGAMMAIX", # model
#                           f="a", # best tree and bootstrap
#                           p=1234, # random number seed
#                           x=2345, # random seed for rapid bootstrapping
#                           N=100, # number of bootstrap replicates
#                           file="alignment", # name of output files
#                           exec="raxmlHPC-PTHREADS-SSE3", # name of executable
#                           threads=20
#)

# Import into Phyloseq ----------------------------------------------------

#library(phyloseq)
#library(ggplot2)

#Create phyloseq object first out of the ASV table and the taxonomy
#then use the "merge_phyloseq" function to merge your sample metadata as QIIME map
#can also merge an ASV tree made using DECIPHER and phangorn

#samples.out <- rownames(seqtab)

#map = "metadata.txt"
#seqtab = readRDS("seqtab.rds")
#taxa = readRDS("taxa.rds")
#physeq_no_metadata <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
#          tax_table(taxa), phy_tree(fitGTR$tree)) 
#sample_metadata = import_qiime_sample_data(map)
#physeq = merge_phyloseq(physeq_no_metadata, sample_metadata)
#physeq
#save the final phyloseq object to disk
#save(physeq, file = "physeq.phyloseq")

#Convert ASV table from phyloseq object to a matrix, the turn it into a biom file below
#otu <- as(otu_table(physeq), "matrix")

#Convert ASV table "seqtab" to a biom table and save to disk
#library(biomformat)
#st.biom <- make_biom(t(otu)) #if ASV table is matrix as above, replace 'seqtab' with matrix ('otu')
#write_biom(st.biom, "physeq.biom")

#Add sample and taxonomy (observation) metadata to the biom table you just made (using biom-format package in terminal)
#follow this documentation (http://biom-format.org/documentation/adding_metadata.html) 
#to see how the metadata and tax table should be formatted
#biom add-metadata -i min_sparse_otu_table.biom -o table.w_md.biom --observation-metadata-fp obs_md.txt --sample-metadata-fp sam_md.txt
#--sc-separated taxonomy --observation-header OTUID,taxonomy

#Create multi-FASTA file with all ASVs and save to disk
#uniquesToFasta(getUniques(otu), fout="rep_seqs.fasta", ids=paste0("ASV", seq(length(getUniques(otu)))))

#Create a fasta file for each sample after processing in DADA2
#require(ShortRead)
#x = row.names(otu)
#for(i in seq(row.names(otu))) {
#ids <- paste0("s", i, "_", seq(rowSums(otu)[i]))
#seqs.re_replicated <- rep(colnames(otu), times=otu[i,])
#writeFasta(object = ShortRead(sread = DNAStringSet(seqs.re_replicated),
#                              id = BStringSet(ids)),
#                              file = paste0("~/Desktop/Preterm-antibiotics/First paper/Qiime2_Analysis/FASTAS/", i ,".fasta"), width = 20000)
#}

