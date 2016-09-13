library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

###############
#Similar to the way that R displays vectors, [[1]] means that R is showing the first element of a list
#To extract an element from a list, you have to use [[]] double square brackets. 

pangram <- "The quick brown fox jumps over the lazy dog"
strsplit(pangram, " ")[[1]][1]

###############


chromHMM_file = fread("/home/surya/Desktop/scripts/data/chrom_impute_chromHMM_data/E008_25_imputed12marks_dense.bed", sep = "\t")
head(chromHMM_file)

chromHMM_file$states <- str_replace(chromHMM_file$V4, "\\d+_", "") 
chromHMM_file$elm_lgth <- chromHMM_file$V3 - chromHMM_file$V2

state_levels <- chromHMM_file$states %>% as.factor %>% levels

elem_sum_table <- chromHMM_file %>% group_by(states) %>% summarize(elem_sum = sum(elm_lgth))

elem_sum_table

fishers_enrichment_file <- fread("/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/fishers_enrichment_result/hyper_chromHMM_fishers_enrichment_result_sorted.txt", sep = "\t", header = T)
#fishers_enrichment_file <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/TF_fishers_enrichment_result.txt", sep = "\t", header = T)
read_significant_DMR <- fread("/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/hyper_dmrs.txt", sep = "\t")
#read_significant_DMR <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/DMR_data.bed", sep = "\t")

names(read_siginificant_DMR) <- c('chr','start','end','CpG','group1','group2','group_name','group1_merged_percent','group2_merged_percent' )
#names(read_significant_DMR) <- c('chr','start','end','q_val','diff','CpG','group1','group2')
bed_g1_g2 <- read_significant_DMR

fishers_enrichment_file %>% head
fishers_enrichment_file$chrom_state <- str_replace(fishers_enrichment_file$ChromHMM_state, "\\d+_", "") 
fishers_enrichment_file


-log10(fisher_counts$Enrichment_pval) %>% hist(breaks = 25)



fisher_counts <- merge(fishers_enrichment_file, elem_sum_table, by.x = "chrom_state", by.y = "states")
fisher_counts$norm_sig_dmr_hits <- with(fisher_counts, Significant_DMR_hits/(elem_sum/1000000))
fisher_counts$norm_sig_dmr_hits
fisher_counts$chrom_state
fisher_counts %>% select(chrom_state, norm_sig_dmr_hits, Enrichment_pval) %>% head(25)
fisher_counts %>% names
names(fisher_counts) 



# TF_sig_hits_barplot <- ggplot(fisher_counts, aes(x=chrom_state, y = norm_sig_dmr_hits, fill= -log10(Enrichment_pval)))+
#   geom_bar(stat = "identity") + xlab("Chromatin State in H9 cells") + ylab("Enrichment % of Hypermethylated DMRs in H9 ES cells") + coord_flip() + theme_cfg_1 
# TF_sig_hits_barplot + scale_fill_gradientn(colors = topo.colors(10), limits = c(0, 60))

fisher_counts_thresh <- fisher_counts
fisher_counts_thresh$log_Enrichment_pvalue <- -log10(fisher_counts_thresh$Enrichment_pval)

fisher_counts_thresh$log_Enrichment_pvalue[1] <- 30

pdf(file="/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/fishers_enrichment_result/hyper_sig_DMR_chromHMM_dist.pdf")

TF_sig_hits_barplot <- ggplot(fisher_counts_thresh, aes(x=chrom_state, y = norm_sig_dmr_hits, fill= log_Enrichment_pvalue))+
  geom_bar(stat = "identity") + xlab("Chromatin State in H9 Cells") + ylab("Hyper-methylated DMRs per Mb of chromatin state in H9 ES Cells") + coord_flip() + theme_cfg_1 
 TF_sig_hits_barplot + scale_fill_gradientn(name="-log10(Fishers_pvalue)", colors = topo.colors(5), limits = c(0, 30))  

dev.off()




#Hypomethylated_DMRs
chromHMM_file = fread("/home/surya/Desktop/scripts/data/chrom_impute_chromHMM_data/E008_25_imputed12marks_dense.bed", sep = "\t")
head(chromHMM_file)

chromHMM_file$states <- str_replace(chromHMM_file$V4, "\\d+_", "") 
chromHMM_file$elm_lgth <- chromHMM_file$V3 - chromHMM_file$V2

state_levels <- chromHMM_file$states %>% as.factor %>% levels

elem_sum_table <- chromHMM_file %>% group_by(states) %>% summarize(elem_sum = sum(elm_lgth))

elem_sum_table

fishers_enrichment_file <- fread("/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/fishers_enrichment_result/hypo_chromHMM_fishers_enrichment_result_sorted.txt", sep = "\t", header = T)
#fishers_enrichment_file <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/TF_fishers_enrichment_result.txt", sep = "\t", header = T)
read_significant_DMR <- fread("/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/hypo_dmrs.txt", sep = "\t")
#read_significant_DMR <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/DMR_data.bed", sep = "\t")

names(read_siginificant_DMR) <- c('chr','start','end','CpG','group1','group2','group_name','group1_merged_percent','group2_merged_percent' )
#names(read_significant_DMR) <- c('chr','start','end','q_val','diff','CpG','group1','group2')
bed_g1_g2 <- read_significant_DMR

fishers_enrichment_file %>% head
fishers_enrichment_file$chrom_state <- str_replace(fishers_enrichment_file$ChromHMM_state, "\\d+_", "") 
fishers_enrichment_file


-log10(fisher_counts$Enrichment_pval) %>% hist(breaks = 25)



fisher_counts <- merge(fishers_enrichment_file, elem_sum_table, by.x = "chrom_state", by.y = "states")
fisher_counts$norm_sig_dmr_hits <- with(fisher_counts, Significant_DMR_hits/(elem_sum/1000000))
fisher_counts$norm_sig_dmr_hits
fisher_counts$chrom_state
fisher_counts %>% select(chrom_state, norm_sig_dmr_hits, Enrichment_pval) %>% head(25)
fisher_counts %>% names
names(fisher_counts) 



# TF_sig_hits_barplot <- ggplot(fisher_counts, aes(x=chrom_state, y = norm_sig_dmr_hits, fill= -log10(Enrichment_pval)))+
#   geom_bar(stat = "identity") + xlab("Chromatin State in H9 cells") + ylab("Enrichment % of Hypermethylated DMRs in H9 ES cells") + coord_flip() + theme_cfg_1 
# TF_sig_hits_barplot + scale_fill_gradientn(colors = topo.colors(10), limits = c(0, 60))

fisher_counts_thresh <- fisher_counts
fisher_counts_thresh$log_Enrichment_pvalue <- -log10(fisher_counts_thresh$Enrichment_pval)
fisher_counts_thresh[fisher_counts_thresh$Significant_DMR_hits ==0,]
fisher_counts_thresh$log_Enrichment_pvalue

pdf(file="/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/fishers_enrichment_result/hypo_sig_DMR_chromHMM_dist.pdf")

TF_sig_hits_barplot <- ggplot(fisher_counts_thresh, aes(x=chrom_state, y = norm_sig_dmr_hits, fill= log_Enrichment_pvalue))+
  geom_bar(stat = "identity") + xlab("Chromatin State in H9 Cells") + ylab("Hypo-methylated DMRs per Mb of chromatin state in H9 ES Cells") + coord_flip() + theme_cfg_1 
 TF_sig_hits_barplot + scale_fill_gradientn(name="-log10(Fishers_pvalue)", colors = topo.colors(5))  

dev.off()






############################
############################

#Here below are just a rough test:
#In fact, are the codes produced for the genomics_DMR_project 3:

###########################
############################




chromHMM_file = fread("/home/surya/Desktop/scripts/data/chrom_impute_chromHMM_data/E008_25_imputed12marks_dense.bed", sep = "\t")
head(chromHMM_file)

fishers_enrichment_file <- fread("/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/fishers_enrichment_result/hyper_chromHMM_fishers_enrichment_result_sorted.txt", sep = "\t", header = T)
#fishers_enrichment_file <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/TF_fishers_enrichment_result.txt", sep = "\t", header = T)
read_significant_DMR <- fread("/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/hyper_dmrs.txt", sep = "\t")
#read_significant_DMR <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/DMR_data.bed", sep = "\t")

names(read_siginificant_DMR) <- c('chr','start','end','CpG','group1','group2','group_name','group1_merged_percent','group2_merged_percent' )
#names(read_significant_DMR) <- c('chr','start','end','q_val','diff','CpG','group1','group2')
bed_g1_g2 <- read_significant_DMR


Total_DMR_count <- nrow(read_significant_DMR)
Total_DMR_count
fishers_enrichment_file$sig_percent_hits <- with(fishers_enrichment_file, (Significant_DMR_hits/Total_DMR_count) * 100) 


## make other extra anlysis plots
theme_cfg <- theme(
  #axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

theme_cfg_1 <-
theme(panel.background = element_rect(colour = "black"),
			axis.text=element_text(color="black",size=10),
	        axis.title=element_text(color="black",size=12,face="bold"),
	        plot.title = element_text(size = rel(1.3), colour = "darkblue", face="bold"))


pdf(file="/home/surya/Desktop/scripts/data/metilene_dmr_data/Project_3_DMR/Fisher_test_DMR_plot.pdf")
#sig_pval_indices <- fishers_enrichment_file$Enrichment_pval < 0.0001
#sig_hits_df <- fishers_enrichment_file[sig_pval_indices,]
sig_hits_df <- fishers_enrichment_file
sig_hits_df$log_Enrichment_pvalue <- -log10(sig_hits_df$Enrichment_pval)
max_limit <- max(sig_hits_df$Enrichment_pval)
min(sig_hits_df$Enrichment_pval)



-log10(sig_hits_df$Enrichment_pval)

sig_hits_df$log_Enrichment_pvalue %>% head
sig_hits_df$Enrichment_pval %>% head

sig_hits_df
TF_sig_hits_barplot <- ggplot(sig_hits_df, aes(x=ChromHMM_state, y = sig_percent_hits, fill= log_Enrichment_pvalue))+
  geom_bar(stat = "identity") + xlab("Chromatin State in H9 cells") + ylab("Enrichment % of Hypermethylated DMRs in H9 ES cells") + coord_flip() + theme_cfg_1 
TF_sig_hits_barplot

dev.off()



pdf("/home/surya/Desktop/scripts/data/metilene_dmr_data/Project_3_DMR/Ratio_DMR_hits_plot.pdf")
sig_hits_df$ratio_sig_to_background <- 0 
sig_hits_df$ratio_sig_to_background <- sig_hits_df$Significant_DMR_hits / sig_hits_df$Background_dmr_Hits

ratio_dmr_barplot <- ggplot(sig_hits_df, aes(x=TF_name, y = ratio_sig_to_background, fill= Enrichment_pval))+
  geom_density(stat = "identity", colour = "red") + xlab("Transcription Factor") + ylab(" Ratio(Significant_DMR_hits/Background_DMR_Hits)") + coord_flip() + theme_cfg_1
ratio_dmr_barplot

dev.off()



hist_pval<- hist(log_pval, main = "Fishers pvalue dist. for TF Enrichment",
	xlab= "Log(Fishers pvalue for Enrichment)",
	ylab="Transcription Factor count", 
    border="blue", 
    col="green",
    las=1, 
    breaks=20)
hist_pval


pdf("/home/surya/Desktop/scripts/data/metilene_dmr_data/Project_3_DMR/Pvalue_bonferroni_plot.pdf")
par(mfrow=c(1,2))
log_pval<-log(fishers_enrichment_file$Enrichment_pval)
hist_pval<- hist(log_pval, main = "Fishers pvalue dist. for TF Enrichment",
	xlab= "Log(Fishers pvalue for Enrichment)",
	ylab="Transcription Factor count", 
    border="blue", 
    col="green",
    las=1, 
    breaks=20)
hist_pval

log_bonferroni <- log(fishers_enrichment_file$Enrichment_Bonferonni_Adj)
hist_bonferroni <- hist(log_bonferroni, main = "Bonferroni correction dist. for TF Enrichment",
	xlab= "Log(Bonferroni Adjusted pvalue)",
	ylab="Transcription Factor count", 
    border="blue", 
    col="steelblue",
    las=1, 
    breaks=20)
hist_bonferroni

dev.off()


##Further analysis considering (<0.25 as hypomethylated, >= 0.25  & <=0.75 as intermediate DMR:
#and > 0.75 as Hypermethylated DMRs
read_significant_DMR$H9_ESC_DMR_STATE <- ""
read_significant_DMR$SM_CELL_DMR_STATE <- ""
read_significant_DMR

#For ESC_DMR
#Just finding the indices for marking the state of DMRs:
hypo_indices  <- which(read_significant_DMR$group1 < 0.25)
read_significant_DMR[hypo_indices,]
intermed_indices  <- which(read_significant_DMR$group1 >= 0.25 & read_significant_DMR$group1 <= 0.75)
read_significant_DMR[intermed_indices,]
hyper_indices  <- which(read_significant_DMR$group1 > 0.75)
read_significant_DMR[hyper_indices,]


#Naming the DMRs:
hypo_indices  <- which(read_significant_DMR$group1 < 0.25)
read_significant_DMR[hypo_indices,]$H9_ESC_DMR_STATE <- "Hypo DMRs"
intermed_indices  <- which(read_significant_DMR$group1 >= 0.25 & read_significant_DMR$group1 <= 0.75)
read_significant_DMR[intermed_indices,]$H9_ESC_DMR_STATE <- "Inter DMRs"
hyper_indices  <- which(read_significant_DMR$group1 > 0.75)
read_significant_DMR[hyper_indices,]$H9_ESC_DMR_STATE <- "Hyper DMRs"
read_significant_DMR


#For SM_cell_state_DMRs:
#Just finding the indices for marking the state of DMRs:
hypo_indices  <- which(read_significant_DMR$group2 < 0.25)
read_significant_DMR[hypo_indices,]
intermed_indices  <- which(read_significant_DMR$group2 >= 0.25 & read_significant_DMR$group2 <= 0.75)
read_significant_DMR[intermed_indices,]
hyper_indices  <- which(read_significant_DMR$group2 > 0.75)
read_significant_DMR[hyper_indices,]


#Naming the DMRs:
hypo_indices  <- which(read_significant_DMR$group2 < 0.25)
read_significant_DMR[hypo_indices,]$SM_CELL_DMR_STATE <- "Hypo DMRs"
read_significant_DMR %>% 	head
intermed_indices  <- which(read_significant_DMR$group2 >= 0.25 & read_significant_DMR$group2 <= 0.75)
read_significant_DMR[intermed_indices,]$SM_CELL_DMR_STATE <- "Inter DMRs"
hyper_indices  <- which(read_significant_DMR$group2 > 0.75)
read_significant_DMR[hyper_indices,]$SM_CELL_DMR_STATE <- "Hyper DMRs"
read_significant_DMR[which(read_significant_DMR$SM_CELL_DMR_STATE == ""),]
nrow(read_significant_DMR)


pdf(file="/home/surya/Desktop/scripts/data/metilene_dmr_data/Project_3_DMR/Piechart_ESC_vs_SM_cell.pdf")
p1 <- ggplot(read_significant_DMR, aes(x=factor(1), fill=H9_ESC_DMR_STATE))+
  geom_bar(width = 1)+ coord_polar("y") + theme_cfg
p1

p2 <- ggplot(read_significant_DMR, aes(x=factor(1), fill=SM_CELL_DMR_STATE))+
  geom_bar(width = 1) + coord_polar("y") + theme_cfg
p2
dev.off()

pdf(file="/home/surya/Desktop/scripts/data/metilendde_dmr_data/Project_3_DMR/test.pdf")
grid.arrange(p3, arrangeGrob(p1,p2),ncol=2)
dev.off()



## make other extra anlysis plots
theme_cfg <-
theme(panel.background = element_rect(colour = "black"),
			axis.text=element_text(color="black",size=10),
	        axis.title=element_text(color="black",size=12,face="bold"),
	        plot.title = element_text(size = rel(1.3), colour = "darkblue", face="bold"))


## Plot statistics
pdf(file="/home/surya/Desktop/scripts/data/metilene_dmr_data/Project_3_DMR/ESC_vs_SM_cell_methylation_diff.pdf")
#difference histogram
p1 <- ggplot(read_significant_DMR, aes(x=diff)) + geom_histogram(binwidth=0.04, fill='darkgrey', color='red') + xlab("Mean Methylation Difference ( ESC - SM cell )") + ylab("DMR count") + scale_x_continuous(limits=c(-1,1)) + theme_cfg
p1
dev.off()

#q_val vs difference
pdf(file="/home/surya/Desktop/scripts/data/metilene_dmr_data/Project_3_DMR/ESC_vs_SM_cell_qval.pdf")
p2 <- ggplot(read_significant_DMR, aes(x=diff, y=q_val)) + geom_point(alpha=.5,color="blue") + scale_y_log10() + xlab("Mean Methylation Difference ( ESC - SM cell )") + ylab("q-value") + theme_cfg
p2
dev.off()

melt_bed_g1_g2 <- bed_g1_g2
names(melt_bed_g1_g2) <- c('chr','start','end','q_val','diff','CpG','H9_ESC','SM_Cell')
merged_df <- melt(melt_bed_g1_g2, id.vars = c('chr','start','end','q_val','diff','CpG'))
names(merged_df) <- c('chr','start','end','q_val','diff','CpG', 'Cell_Type', 'value')
cdat <- ddply(merged_df, "Cell_Type", summarise, value.mean=mean(value))

#Mean methylation diff on 2 cell types:
pdf(file="/home/surya/Desktop/scripts/data/metilene_dmr_data/Project_3_DMR/ESC_vs_SM_cell_methylation_density.pdf")
p3 <-ggplot2::ggplot(merged_df, aes(x=value, y=..count../sum(..count..), group=Cell_Type, colour=Cell_Type, fill=Cell_Type)) +
    ggplot2::geom_density(alpha= 0.3) + geom_vline(data=cdat, aes(xintercept=value.mean,  colour=Cell_Type),
               linetype="dashed", size=1) + xlab("Mean methylation distribution for DMRs") + ylab("Density") + theme_cfg
p3
dev.off()

