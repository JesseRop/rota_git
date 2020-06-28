##############################################################################################################
# An R.script for generating SNP frequencies from sequence alignments
# It takes an alignment from CLC exported in Phylip format
# change the input file at the load data step to match your input file
##############################################################################################################
#load applicable functions and packages
library(data.table) #fread
library(tidyverse)  #%>% and unite
library(ggpubr)     #plot
library(genetics)
library(LDheatmap)

#Loading phylogenetic packages
library(devtools)
library(adegenet)  
library(NAM)    
library(dendextend)
library(circlize)
library(snow)
library(doSNOW)
library(parallel)
library(colorspace)
library(Biostrings)
library(seqinr)
library(msa)
library(ggtree)
library(treeio)
library(phangorn)
library(ape)


#clear environment
rm(list = ls())

setwd("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/rotavirus analysis/R.project - SNP Frrequency/")

#-------------------------------------------------------------------------------------------------------------

detectCores()
cl <-makeCluster(4, type = "SOCK")
registerDoSNOW(cl)
#-------------------------------------------------------------------------------------------------------------
#function for checking whether all bases are identical per locus
#load functions
source("scripts - R.functions/function_R - identical_bases_per_locus.R")   #checking for polymorphic loci
source("scripts - R.functions/function_R - snpFreq_per_locus.R")           #snp freq per column/locus

#*************************************************************************************************************
#load ligase data 
#*************************************************************************************************************
#skip line 1 as well as empty lines and select colums 1 and 2
fastaIN <- fread("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/fut2 phylip from clc/KBC_FUT2_CDS_Consensus_alignment_QCd.phy", skip =1, header = F, blank.lines.skip = T) 

#split sequence column into multiple columns
fastaIN_seqs <- fastaIN[,tstrsplit(V3, "")]  #require data.table

#select sequence names
fastaIN_names <- fastaIN[,.(V1)]

#bind the two datasets
fastaIN <- cbind(fastaIN_names, fastaIN_seqs)

#remame column 1 as it matches column 2
names(fastaIN)[1] <- c("V0")

#convert to data.frame
fastaIN <- as.data.frame(fastaIN)

#-------------------------------------------------------------------------------------------------------------
#replace all gaps with NA so that these gaps are not evaluated 
fastaIN[fastaIN == "-"] <- NA

#create new dataframe with just the sequence IDs, to be used in the subsequent loop
fastaIN_v1 <- as.data.frame(fastaIN[1])

#keep only columns that represent haplotypes (those that dont match)
for(colN in 2:dim(fastaIN)[2]) {
  if(all_identical(fastaIN[,colN]) == F) {
    fastaIN_v1 <- as.data.frame(cbind(fastaIN_v1, fastaIN[colN]))
  }
}


##Checking HWE for all markers
gtype_object_all <-  fastaIN_v1
gtype_object_all[,-1] <- as.data.frame(lapply(gtype_object_all[,-1], function(x) ifelse(x == "A", "A/A", ifelse(x == "G", "G/G", ifelse(x == "C", "C/C", ifelse(x == "T", "T/T", ifelse(x == "R", "A/G", ifelse(x == "Y", "C/T", ifelse(x == "S", "G/C", ifelse(x == "W", "A/T", NA))))))))))
gtype_object_all <- makeGenotypes(gtype_object_all)
lapply(gtype_object_all[,-1], function(x)HWE.exact(unlist(x)))

##Plotting LD for all the fut2 markers
LDplot(genetics::LD(gtype_object_all))
LDtable(genetics::LD(gtype_object_all))

#Generating a table of the allele frequencies, hwe and ld to rs601338 for all the fut2 markers
fut2_allele_freqs = bind_rows(lapply(gtype_object_all[,-1], function(x) c(summary(x)$allele.names, summary(x)$allele.freq[c(3,4)]))) %>% t() %>% `colnames<-`(c("major_allele", "minor_allele", "major_af", "minor_af")) %>% as.data.frame()%>% rownames_to_column(var = "cds_pos") %>% mutate_at(vars(cds_pos), ~str_remove(.,"V"))
fut2_hwe = bind_rows(lapply(gtype_object_all[,-1], function(x) HWE.exact(x)[3]), .id = "cds_pos")
fut2_ld_to_rs601338 = bind_rows(lapply(gtype_object_all[,-1], function(x) genetics::LD(gtype_object_all[,"V461"],x)[2:8]), .id = "cds_pos")

##########################Phylo###########################
kbc_metadata = read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/R processed tables/sequenced_ethnicities.csv", na = c("", "<NA>", "NA"))

kbc_fut2_set <-  fread("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/fut2 phylip from clc/KBC_FUT2_CDS_Consensus_alignment_QCd.phy", skip =1, header = F, blank.lines.skip = T) %>% dplyr::select(V1, V3) %>% rename(V1 = "serial", V3 = "fut2") %>% tibble::deframe()  %>% DNAStringSet()
kbc_fut2_aln  = msa(kbc_fut2_set)
# msaPrettyPrint(p_aln, y=c(164, 213), output="asis",showNames="none", showLogo="none", askForOverwrite=FALSE)
# kbc_fut2_aln_sir <- msaConvert(kbc_fut2_aln, type="seqinr::alignment")
contrast_fut2 = matrix(data = c(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1, 1,0,1,0, 0,1,0,1, 0,1,1,0, 1,0,0,1),ncol = 4, byrow = TRUE) 
dimnames(contrast_fut2) = list(c("a","c","g","t","r","y","s","w"), c("a", "c", "g", "t"))
kbc_fut2_aln_ph <- as.phyDat(kbc_fut2_aln, type="USER", contrast = contrast_fut2)
mt = modelTest(kbc_fut2_aln_ph)
kbc_fut2_dist = dist.ml(kbc_fut2_aln_ph)
kbc_fut2_tree_nj = nj(kbc_fut2_dist)
fut2_fit <- pml(kbc_fut2_tree_nj, kbc_fut2_aln_ph)
fut2_fitJC <- optim.pml(fut2_fit, model = "JC", rearrangement = "stochastic")
fut2_bs <- bootstrap.pml(fut2_fitJC, bs=100, optNni=TRUE, control = pml.control(trace=0))
kbc_fut2_ggtree = ggtree(fut2_bs)
kbc_fut2_ggtree = ggtree(as.treedata(fut2_fitJC))
kbc_fut2_ggtree = ggtree(as.treedata(plotBS(midpoint(fut2_fitJC$tree), fut2_bs, p = 50, type="p")))
write.tree(plotBS(midpoint(fut2_fitJC$tree), fut2_bs, p = 50, type="p"), file="fut2_bs.tre")
kbc_fut2_ggtree = ggtree(read.newick("fut2_bs.tre"))



kbc_fut2_tree_annot = kbc_metadata %>% dplyr::select(serial, ct_ethnic, ethnicity_grouped, final_sickle, final_thall)

kbc_fut2_ggtree_annot = kbc_fut2_ggtree %<+% kbc_fut2_tree_annot 

# kbc_fut2_ggtree_annot + geom_tiplab(aes(linetype=ptype_all_exc)) + geom_tippoint(aes(size=hospitalization_days, shape=secretor_status), alpha=0.25)
kbc_fut2_ggtree_annot + geom_tippoint(aes(color=ct_ethnic, shape = ethnicity_grouped))
kbc_fut2_ggtree_annot + geom_tippoint(aes(color=ptype_all_exc, shape = non_OvsO_jcr)) + geom_facet(panel = "dia dur", data = rota_geno_pheno_cc %>% select(serial, diarrhoea_durn, hospitalization_days), geom = ggstance::geom_barh, aes(x = diarrhoea_durn),  stat = "identity", width = .6)+ geom_facet(panel = "Trait", data = rota_geno_pheno_cc %>% select(serial, diarrhoea_durn, diarrhoea_episode), geom = ggstance::geom_barh, aes(x = diarrhoea_episode),  stat = "identity", width = .6)

snp_t = fastaIN_v1 %>% remove_rownames %>% column_to_rownames(var="V0")
gapChar <- "?"
lsnp_fut2 <- apply(snp_t, 1, function(x) {
  x != snp_t[1,] & x != gapChar & snp_t[1,] != gapChar
})
lsnp_fut2 <- as.data.frame(lsnp_fut2)
lsnp_fut2$pos <- as.numeric(rownames(lsnp_fut2))
lsnp_fut2 <- tidyr::gather(lsnp_fut2, name, value, -pos)
snp_data_fut2 <- lsnp_fut2[lsnp_fut2$value, c("name", "pos")]

kbc_fut2_ggtree_annot + geom_tippoint(aes(color=ethnicity_grouped)) + geom_facet(panel = "SNP", data = snp_data_fut2, geom = geom_point, mapping=aes(x = pos, color = ethnicity_grouped), shape = '|') 

##########################Phylo###########################


##########################Annotation###########################

fut2_cds_snps = fread("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/fut2_snps/fut2_cds_vcf", header = T) ##%>% mutate(rsid = str_extract(INFO, "rs.+(?=;.+)"))
fut2_erythrogene = read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/fut2_snps/fut2_functional_erythrogene.csv")%>% mutate(allele_cds1 = str_extract_all(`Nucleotide Change`, "[0-9]+", simplify = T)[,1], allele_cds2 = str_extract_all(`Nucleotide Change`, "[0-9]+", simplify = T)[,2], allele_cds3 = str_extract_all(`Nucleotide Change`, "[0-9]+", simplify = T)[,3] ) %>%   mutate_if(is.character, list(~na_if(., ""))) %>% pivot_longer(c(allele_cds1, allele_cds2, allele_cds3), names_to = "cds_in_allele", values_to = "pos", values_drop_na = T) ## %>% mutate(dup = duplicated(pos)) %>% pivot_wider(., id_cols = pos, names_from = dup, values_from = `Nucleotide Change`)
## removing duplicates
fut2_erythrogene = fut2_erythrogene %>% distinct(., pos, .keep_all = T)
# read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/fut2_snps/fut2_functional_erythrogene.csv")%>% mutate(allele_cds = str_extract_all(`Nucleotide Change`, "[0-9]+(?=[A-Z]>[A-Z])", simplify = T))


fut2_1kg = fread("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/fut2_snps/fut2_snps_af__1kgg.txt", header = T, skip = 256) %>% mutate_at(vars("global", "ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU",  "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"), ~as.numeric(str_split(., ":", n=3, simplify = T)[,2])/as.numeric(str_split(., ":", n=3, simplify = T)[,1])) %>% separate(INFO, into = c("AA", "AC", "AF", "AFR_AF", "AMR_AF", "AN", "DP", "EAS_AF", "EUR_AF", "EX_TARGET", "NS", "SAS_AF", "SF", "VT"), sep =  ";") %>% dplyr::select(c("ID", "REF", "ALT", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "global", "ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"))
fut2_ensemble_snps = fread("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/fut2_snps/fut2_snps_ensembl.txt", header = T) %>% distinct(`#Uploaded_variation`, .keep_all=T) %>% dplyr::select("#Uploaded_variation", "Allele", "Consequence", "IMPACT","Protein_position", "Amino_acids", "SIFT", "PolyPhen", "AF", "CLIN_SIG")

kbc_cds_pos = t(snp_t) %>% as.data.frame() %>% rownames_to_column(var = "cds_pos") %>%  mutate_at(vars(cds_pos), ~str_remove(.,"V")) %>% mutate(chr_pos = as.numeric(cds_pos)+48702956) %>% dplyr::select(cds_pos,chr_pos) %>% right_join(fut2_allele_freqs,., by = "cds_pos") %>% right_join(fut2_erythrogene,., by=c("pos" = "cds_pos")) %>% right_join(fut2_cds_snps,., by=c("POS"="chr_pos")) %>% right_join(fut2_1kg,., by = "ID") %>% right_join(fut2_ensemble_snps,., by = (c("#Uploaded_variation" = "ID"))) %>% dplyr::select(c("#Uploaded_variation","Allele.y", "REF.x", "ALT.x", "major_allele", "minor_allele", "Nucleotide Change", "Predicted Amino Acid Change", "1000G", "cds_in_allele", "pos",  "Consequence", "IMPACT", "Protein_position", "Amino_acids", "SIFT", "PolyPhen", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "global", "major_af", "minor_af", "ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI")) 
kbc_cds_pos_summ = kbc_cds_pos %>% dplyr::select(c("#Uploaded_variation", "REF.x", "ALT.x", "major_allele", "minor_allele","Allele.y", "Nucleotide Change", "cds_in_allele", "pos",  "Consequence", "Amino_acids", "SIFT", "PolyPhen", "AFR_AF", "global", "major_af", "minor_af", "LWK"))
##########################End of Annotation###########################


##########333Phylo
gtype_phylo = fastaIN
gtype_phylo[,-1] <- as.data.frame(lapply(gtype_phylo[,-1], function(x) ifelse(x == "A", "A/A", ifelse(x == "G", "G/G", ifelse(x == "C", "C/C", ifelse(x == "T", "T/T", ifelse(x == "R", "A/G", ifelse(x == "Y", "C/T", ifelse(x == "S", "G/C", ifelse(x == "W", "A/T", NA))))))))))
gtype_phylo <- makeGenotypes(gtype_phylo)
#gtype_phylo_allele_cnt = gtype_phylo %>% map_at(vars(-1), ~allele.count(., allele.name = allele.names(.)[1]))%>% data.frame() %>% mutate_if(~length(unique(.))==1, ~ recode(.,`2` = 0))
gtype_phylo_allele_cnt = gtype_phylo %>% map_at(vars(-1), ~allele.count(., allele.name = allele.names(.)[1]))%>% data.frame() 

phylo_metadata = right_join(kbc_metadata,gtype_phylo_allele_cnt, by = c("serial" = "V0")) %>% mutate_at(vars(ct_ethnic, ethnicity_grouped), ~replace_na(., "MISSING"))

# phylo_matrix = fastaIN %>% dplyr::select(-c(1:12)) %>% mutate_all(., as.factor) %>% mutate_all(., fct_infreq) %>% mutate_all(., as.numeric) %>% mutate_all(., ~recode(., `1` = 0, `2`=1, `3`=2)) %>% data.matrix()
phylo_matrix = phylo_metadata %>% dplyr::select(-c(1:12)) %>% data.matrix() %>% reference()
gdist1 = Gdist(phylo_matrix, method = 6)

obj = df2genind(phylo_matrix, ploidy = 2, sep = "/t")

#grp = find.clusters(obj, max.n = 50, n.pca = 200, scale = F)
# plot(grp)

hier_cluster_m1 = (hclust(gdist1, method = 'ward.D'))

cutree.8 = cutree(hier_cluster_m1, k = 23)

summary(as.factor(cutree(hier_cluster_m1, k = 23)))

dend = as.dendrogram(hclust(gdist1, method = 'ward.D'))

dend = color_branches(dend)
labels(dend) = as.character(unlist(phylo_metadata[,"ct_ethnic"]))[order.dendrogram(dend)]


circlize_dendrogram(dend, dend_track_height = 0.85)


##########End of Phylo


##############Resume sanger seq analysis####################################################3
#create an empty dataframe, to be used in subsequent loop
freqStats_v1 <- data.frame(freq = as.numeric(),
                           percentage = as.numeric(),
                           base = as.character(),
                           locus = as.character())


#compute sample size and Func_SnpFrequency for each nucleotide per locus
for(i in 2:dim(fastaIN_v1)[2]) {
  freqStats_v2 <- Func_SnpFreq(fastaIN_v1[,i])
  freqStats_v2 <- as.data.frame(freqStats_v2)
  freqStats_v2 <- freqStats_v2 %>% mutate(base = row.names(freqStats_v2),
                                          locus = colnames(fastaIN_v1[i]))
  freqStats_v1 <- rbind(freqStats_v1, freqStats_v2)
}

#remove V's from locus
freqStats_v1 <- freqStats_v1 %>% mutate(locus = gsub("V", "", locus),
                                        locus = as.numeric(locus)
)

#filter all locus with 100% freq
consensus_snpFrequency <- freqStats_v1 %>% filter(percentage != 100)
write.csv(t(consensus_snpFrequency), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/consensus_snpFrequency.csv", quote = F, row.names = F)

##non-secretor alleles
erythrogene_nonsecretor_alleles <- c("461", "418", "335", "882", "311", "901", "602", "691", "604")
consensus_snpFrequency_non_sec <- consensus_snpFrequency[which(consensus_snpFrequency$locus %in% erythrogene_nonsecretor_alleles),]
write.csv(t(consensus_snpFrequency_non_sec), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/consensus_snpFrequency_non_sec.csv", quote = F, row.names = F)



##Splitting consensus into 2 alleles
allele1_file <- fastaIN_v1 %>% mutate(allele_name = paste(V0, "_A1", sep = ""))
allele1_file[,c(2:22)] <- lapply(allele1_file[,c(2:22)], function(x) ifelse(x == "Y", "C", (ifelse(x == "R", "A", (ifelse(x == "S", "C", (ifelse(x == "W", "A", x))))))))

allele2_file <- fastaIN_v1 %>% mutate(allele_name = paste(V0, "_A2", sep = ""))
allele2_file[,c(2:22)] <- lapply(allele2_file[,c(2:22)], function(x) ifelse(x == "Y", "T", (ifelse(x == "R", "G", (ifelse(x == "S", "G", (ifelse(x == "W", "T", x))))))))

bi_allelic_file <- rbind(allele1_file, allele2_file)

##merging the nucleotides into 1 column and naming the unique haplotypes for bi-allelic file
uniq_haplotypes <- distinct(bi_allelic_file[,c(2:22)]) %>% unite(merged_haps, 1:21, sep = "", remove = T) %>% mutate(hap_name = paste("haplotype", row_number()))
named_haplotypes <- merge(bi_allelic_file %>% unite(merged_haps, 2:22, sep = "", remove = FALSE), uniq_haplotypes, by = "merged_haps", all = T)
named_haplotypes <- named_haplotypes[order(named_haplotypes$hap_name),]

#create an empty dataframe, to be used in subsequent loop
freqStats_v1 <- data.frame(freq = as.numeric(),
                           percentage = as.numeric(),
                           base = as.character(),
                           locus = as.character())

#compute sample size and Func_SnpFrequency for each nucleotide per locus
for(i in 2:(dim(bi_allelic_file)[2]-1)) {
  freqStats_v2 <- Func_SnpFreq(bi_allelic_file[,i])
  freqStats_v2 <- as.data.frame(freqStats_v2)
  freqStats_v2 <- freqStats_v2 %>% mutate(base = row.names(freqStats_v2),
                                          locus = colnames(bi_allelic_file[i]))
  freqStats_v1 <- rbind(freqStats_v1, freqStats_v2)
}

#remove V's from locus
freqStats_v1 <- freqStats_v1 %>% mutate(locus = gsub("V", "", locus),
                                        locus = as.numeric(locus)
)

#filter all locus with 100% freq
bi_allele_snpFrequency <- freqStats_v1 %>% filter(percentage != 100)
write.csv(t(bi_allele_snpFrequency), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/fut2_bi_allele_snpFrequency.csv", quote = F, row.names = F)

##non-secretor alleles
bi_allele_snpFrequency_non_sec <- bi_allele_snpFrequency[which(bi_allele_snpFrequency$locus %in% erythrogene_nonsecretor_alleles),]
write.csv(t(bi_allele_snpFrequency_non_sec), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/bi_allele_snpFrequency_non_sec.csv", quote = F, row.names = F)

##with MAF greater than 5%
write.csv(bi_allele_snpFrequency[bi_allele_snpFrequency$percentage > 5 & bi_allele_snpFrequency$percentage < 95,], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/bi_allele_snpFrequency_common.csv", quote = F, row.names = F )

########End of altered code for biallelic file
##Displaying and calculating LD between 461 and 772 positions
fastaIN_v1 %>% group_by(V461, V772) %>% summarise(n = n())
gtype_object_461_472 <- fastaIN_v1[,c("V0", "V461", "V772")]
gtype_object_461_472[,c("V461", "V772")] <- lapply(gtype_object_461_472[,c("V461", "V772")], function(x) ifelse(x == "A", "A/A", ifelse(x == "G", "G/G", "G/A")))
gtype_object_461_472 <- makeGenotypes(gtype_object_461_472[,2:3])
LD(gtype_object_461_472)

##Comparing sequencing results with genotyping results
rs601338_gtypes_comparison_seq <- fastaIN_v1[,c("V0", "V461")]
rs601338_gtypes_comparison_seq[,c("V461")] <- ifelse(rs601338_gtypes_comparison_seq[,c("V461")] == "G", "GG", ifelse(rs601338_gtypes_comparison_seq[,c("V461")] == "A", "AA", "GA"))
write.csv(rs601338_gtypes_comparison_seq, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/rs601338_gtypes_comparison_seq.csv", quote = F, row.names = F)

rs601338_taqman_genotyped = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/genotyping results for correlation with sequences/kbc_rs601338_11102019_calls.csv", stringsAsFactors = F, header = T, na.strings =c("", "Undetermined", "Outlier"))
rs601338_gtypes_comparison_seq_taqman = merge(rs601338_gtypes_comparison_seq, rs601338_taqman_genotyped[,c("Sample.Name", "Call")], by.x = "V0", by.y  = "Sample.Name")
rs601338_gtypes_comparison_seq_taqman[,"Call"] <- ifelse(rs601338_gtypes_comparison_seq_taqman[,"Call"] == "GG", "AA", ifelse(rs601338_gtypes_comparison_seq_taqman[,"Call"] == "AA", "GG", "GA"))


#remove temp files
#rm(fastaIN, fastaIN_names, fastaIN_seqs, fastaIN_v1, freqStats_v1, freqStats_v2, colN, i, all_identical, Func_SnpFreq)

##troubleshooting the wrongly genotyped individuals
uniq_consensus_haplotypes <- distinct(fastaIN_v1[,c(2:22)]) %>% unite(merged_haps_het, 1:21, sep = "", remove = T) %>% mutate(hap_name = paste("haplotype", row_number()))
named_consensus_haplotypes <- merge(fastaIN_v1 %>% unite(merged_haps_het, 2:22, sep = "", remove = F), uniq_consensus_haplotypes, by = "merged_haps_het", all = T)

rs601338_taqman_vs_sequencing_all = merge(named_consensus_haplotypes, rs601338_taqman_genotyped[,c("Sample.Name", "Call")], by.x = "V0", by.y  = "Sample.Name")
rs601338_taqman_vs_sequencing_all$Call = ifelse(rs601338_taqman_vs_sequencing_all[,"Call"] == "GG", "AA", ifelse(rs601338_taqman_vs_sequencing_all[,"Call"] == "AA", "GG", "GA"))

##misgenotyped heterozygotes and GG hom
hets_genotyped_as_non_secretors = rs601338_taqman_vs_sequencing_all[rs601338_taqman_vs_sequencing_all[,"V461"] == "R" & rs601338_taqman_vs_sequencing_all[,"Call"] == "AA", ]
GGhom_genotyped_as_het_secretors = rs601338_taqman_vs_sequencing_all[rs601338_taqman_vs_sequencing_all[,"V461"] == "G" & rs601338_taqman_vs_sequencing_all[,"Call"] == "GA", ]

hets_genotyped_as_non_secretors_tr = hets_genotyped_as_non_secretors[,c("V0", "V73","V461","V513")]



##FUT3######################FUT3#######################FUT3
##############################################################################################################
# An R.script for generating SNP frequencies from sequence alignments
# It takes an alignment from CLC exported in Phylip format
# change the input file at the load data step to match your input file
##############################################################################################################
#clear environment

fastaIN <- fread("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/fut3 phylip from clc/KBC FUT3 CDS Consensus alignment.phy", skip =1, header = F, blank.lines.skip = T) 

#split sequence column into multiple columns
fastaIN_seqs <- fastaIN[,tstrsplit(V3, "")]  #require data.table

#select sequence names
fastaIN_names <- fastaIN[,.(V1)]

#bind the two datasets
fastaIN <- cbind(fastaIN_names, fastaIN_seqs)

#remame column 1 as it matches column 2
names(fastaIN)[1] <- c("V0")

#convert to data.frame
fastaIN <- as.data.frame(fastaIN)

#-------------------------------------------------------------------------------------------------------------
#replace all gaps with NA so that these gaps are not evaluated 
fastaIN[fastaIN == "-"] <- NA

#create new dataframe with just the sequence IDs, to be used in the subsequent loop
fastaIN_v1 <- as.data.frame(fastaIN[1])

#keep only columns that represent haplotypes (those that dont match)
for(colN in 2:dim(fastaIN)[2]) {
  if(all_identical(fastaIN[,colN]) == F) {
    fastaIN_v1 <- as.data.frame(cbind(fastaIN_v1, fastaIN[colN]))
  }
}


#create an empty dataframe, to be used in subsequent loop
freqStats_v1 <- data.frame(freq = as.numeric(),
                           percentage = as.numeric(),
                           base = as.character(),
                           locus = as.character())


#compute sample size and Func_SnpFrequency for each nucleotide per locus
for(i in 2:dim(fastaIN_v1)[2]) {
  freqStats_v2 <- Func_SnpFreq(fastaIN_v1[,i])
  freqStats_v2 <- as.data.frame(freqStats_v2)
  freqStats_v2 <- freqStats_v2 %>% mutate(base = row.names(freqStats_v2),
                                          locus = colnames(fastaIN_v1[i]))
  freqStats_v1 <- rbind(freqStats_v1, freqStats_v2)
}

#remove V's from locus
freqStats_v1 <- freqStats_v1 %>% mutate(locus = gsub("V", "", locus),
                                        locus = as.numeric(locus)
)

#filter all locus with 100% freq
consensus_snpFrequency <- freqStats_v1 %>% filter(percentage != 100)

write.csv(t(consensus_snpFrequency), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/consensus_snpFrequency.csv", quote = F, row.names = F)

##non-secretor alleles
erythrogene_lewis_neg_alleles <- c("59", "508", "202", "314", "1067", "13", "484", "667", "13", "484", "808", "858", "47", "974", "975", "732", "1007", "61", "980", "882", "667", "1067")
consensus_snpFrequency_lewis_neg <- consensus_snpFrequency[which(consensus_snpFrequency$locus %in% erythrogene_lewis_neg_alleles),]
write.csv(t(consensus_snpFrequency_lewis_neg), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/consensus_snpFrequency_non_sec.csv", quote = F, row.names = F)

##Splitting consensus into 2 alleles
allele1_file <- fastaIN_v1 %>% mutate(allele_name = paste(V0, "_A1", sep = ""))
allele1_file[,c(2:26)] <- lapply(allele1_file[,c(2:26)], function(x) ifelse(x == "Y", "C", (ifelse(x == "R", "A", (ifelse(x == "S", "C", (ifelse(x == "W", "A", (ifelse(x == "K", "G", x))))))))))

allele2_file <- fastaIN_v1 %>% mutate(allele_name = paste(V0, "_A2", sep = ""))
allele2_file[,c(2:26)] <- lapply(allele2_file[,c(2:26)], function(x) ifelse(x == "Y", "T", (ifelse(x == "R", "G", (ifelse(x == "S", "G", (ifelse(x == "W", "T", (ifelse(x == "K", "T", x))))))))))

bi_allelic_file <- rbind(allele1_file, allele2_file)

##merging the nucleotides into 1 column and naming the unique haplotypes for bi-allelic file
uniq_haplotypes <- distinct(bi_allelic_file[,c(2:26)]) %>% unite(merged_haps, 1:25, sep = "", remove = T) %>% mutate(hap_name = paste("haplotype", row_number()))
named_haplotypes <- merge(bi_allelic_file %>% unite(merged_haps, 2:26, sep = "", remove = FALSE), uniq_haplotypes, by = "merged_haps", all = T)
named_haplotypes <- named_haplotypes[order(named_haplotypes$hap_name),]

#create an empty dataframe, to be used in subsequent loop
freqStats_v1 <- data.frame(freq = as.numeric(),
                           percentage = as.numeric(),
                           base = as.character(),
                           locus = as.character())

#compute sample size and Func_SnpFrequency for each nucleotide per locus
for(i in 2:(dim(bi_allelic_file)[2]-1)) {
  freqStats_v2 <- Func_SnpFreq(bi_allelic_file[,i])
  freqStats_v2 <- as.data.frame(freqStats_v2)
  freqStats_v2 <- freqStats_v2 %>% mutate(base = row.names(freqStats_v2),
                                          locus = colnames(bi_allelic_file[i]))
  freqStats_v1 <- rbind(freqStats_v1, freqStats_v2)
}

#remove V's from locus
freqStats_v1 <- freqStats_v1 %>% mutate(locus = gsub("V", "", locus),
                                        locus = as.numeric(locus)
)


#filter all locus with 100% freq
bi_allele_snpFrequency <- freqStats_v1 %>% filter(percentage != 100)
write.csv(t(bi_allele_snpFrequency), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut3/fut3_bi_allele_snpFrequency.csv", quote = F, row.names = F)

##lewis_negative alleles
bi_allele_snpFrequency_lewis_neg <- bi_allele_snpFrequency[which(bi_allele_snpFrequency$locus %in% erythrogene_lewis_neg_alleles),]
write.csv(t(bi_allele_snpFrequency_lewis_neg), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut3/bi_allele_snpFrequency_non_sec.csv", quote = F, row.names = F)

##with MAF greater than 5%
write.csv(bi_allele_snpFrequency[bi_allele_snpFrequency$percentage > 5 & bi_allele_snpFrequency$percentage < 95,], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut3/bi_allele_snpFrequency_common.csv", quote = F, row.names = F )

##writing out file with common SNPs
write.csv(fastaIN_v1[,c(paste("V",c(unique(bi_allele_snpFrequency[bi_allele_snpFrequency$percentage > 5 & bi_allele_snpFrequency$percentage < 95, "locus"])), sep = ""))], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut3/Haplotypes_of_common_snps.csv", quote = F, row.names = F )


########End of altered code for biallelic file

##Displaying and calculating LD between 461 and 772 positions
#fastaIN_v1 %>% group_by(V461, V772) %>% summarise(n = n())
#gtype_object_461_472 <- fastaIN_v1[,c("V0", "V461", "V772")]
#gtype_object_461_472[,c("V461", "V772")] <- lapply(gtype_object_461_472[,c("V461", "V772")], function(x) ifelse(x == "A", "A/A", ifelse(x == "G", "G/G", "G/A")))
#gtype_object_461_472 <- makeGenotypes(gtype_object_461_472[,2:3])
#LD(gtype_object_461_472)

##Checking HWE for all markers
fut3_gtype_object_all <-  fastaIN_v1
fut3_gtype_object_all[,-1] <- as.data.frame(lapply(fut3_gtype_object_all[,-1], function(x) ifelse(x == "A", "A/A", ifelse(x == "G", "G/G", ifelse(x == "C", "C/C", ifelse(x == "T", "T/T", ifelse(x == "R", "A/G", ifelse(x == "Y", "C/T", ifelse(x == "S", "G/C", ifelse(x == "W", "A/T", ifelse(x == "K", "G/T", NA)))))))))))
fut3_gtype_object_all <- makeGenotypes(fut3_gtype_object_all)
lapply(fut3_gtype_object_all[,-1], function(x)HWE.exact(unlist(x)))
HWE.exact(fut3_gtype_object_all[,3])

#lapply(gtype_object_all[,-1], function(x)HWE.exact(unlist(x)))


##Comparing sequencing results with genotyping results
rs28362458_59_gtypes_comparison_seq <- fastaIN_v1[,c("V0", "V13", "V59")]
rs28362458_59_gtypes_comparison_seq[,c("V13", "V59")] <- lapply(rs28362458_59_gtypes_comparison_seq[,c("V13", "V59")], function(x) ifelse(x == "G", "GG", ifelse(x == "A", "AA", ifelse(x == "R", "GA", ifelse(x == "T", "TT", ifelse(x == "K", "TG", x)))))) 
write.csv(rs28362458_59_gtypes_comparison_seq, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/SEQUENCING/sequencing results/r_processed_seq_results/fut2/rs28362458_59_gtypes_comparison_seq.csv", quote = F, row.names = F)


rs28362458_taqman_genotyped = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/genotyping results for correlation with sequences/rota_28362458_kbc_calls.csv", stringsAsFactors = F, header = T)
rs28362459_taqman_genotyped = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/KBC RSV household/lab work/genotyping results for correlation with sequences/rota_rs28362459_kbc_calls.csv", stringsAsFactors = F, header = T)

rs28362458_59_gtypes_comparison_seq = merge(rs28362458_59_gtypes_comparison_seq, rs28362458_taqman_genotyped[,c("serial", "rs28362458")], by.x = "V0", by.y  = "serial")
rs28362458_59_gtypes_comparison_seq = merge(rs28362458_59_gtypes_comparison_seq, rs28362459_taqman_genotyped[,c("serial", "rs28362459")], by.x = "V0", by.y  = "serial")

#rs601338_gtypes_comparison_seq_taqman[,"Call"] <- ifelse(rs601338_gtypes_comparison_seq_taqman[,"Call"] == "GG", "AA", ifelse(rs601338_gtypes_comparison_seq_taqman[,"Call"] == "AA", "GG", "GA"))


#remove temp files
#rm(fastaIN, fastaIN_names, fastaIN_seqs, fastaIN_v1, freqStats_v1, freqStats_v2, colN, i, all_identical, Func_SnpFreq)

##troubleshooting the wrongly genotyped individuals
uniq_consensus_haplotypes <- distinct(fastaIN_v1[,c(2:26)]) %>% unite(merged_haps_het, 1:25, sep = "", remove = T) %>% mutate(hap_name = paste("haplotype", row_number()))
named_consensus_haplotypes <- merge(fastaIN_v1 %>% unite(merged_haps_het, 2:26, sep = "", remove = F), uniq_consensus_haplotypes, by = "merged_haps_het", all = T)

#rs601338_taqman_vs_sequencing_all = merge(named_consensus_haplotypes, rs601338_taqman_genotyped[,c("Sample.Name", "Call")], by.x = "V0", by.y  = "Sample.Name")
#rs601338_taqman_vs_sequencing_all$Call = ifelse(rs601338_taqman_vs_sequencing_all[,"Call"] == "GG", "AA", ifelse(rs601338_taqman_vs_sequencing_all[,"Call"] == "AA", "GG", "GA"))

##misgenotyped heterozygotes and GG hom
#hets_genotyped_as_non_secretors = rs601338_taqman_vs_sequencing_all[rs601338_taqman_vs_sequencing_all[,"V461"] == "R" & rs601338_taqman_vs_sequencing_all[,"Call"] == "AA", ]
#GGhom_genotyped_as_het_secretors = rs601338_taqman_vs_sequencing_all[rs601338_taqman_vs_sequencing_all[,"V461"] == "G" & rs601338_taqman_vs_sequencing_all[,"Call"] == "GA", ]


