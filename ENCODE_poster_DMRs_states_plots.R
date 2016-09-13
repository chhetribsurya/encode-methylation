
#args <- commandArgs(TRUE)
#setwd("/home/surya/Desktop/scripts/data")
#getwd()

library(plyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

#write.table(df, file="foo.bed_g2_g3", quote=F, sep="\t", row.names=F, col.names=F)

#######################
#Similar to the way that R displays vectors, [[1]] means that R is showing the first element of a list
#To extract an element from a list, you have to use [[]] double square brackets. 
pangram <- "The quick brown fox jumps over the lazy dog"
strsplit(pangram, " ")[[1]][1]
###############################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




bed_g1_g2 <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/merged_processed_DMR/analyse_DMR_output_data/merged_processed_DMR/human_ES_cell_final_dmrs.bed", sep="\t")
names(bed_g1_g2) <- c('chr','start','end','CpG','group1','group2','group_name','group1_merged_percent','group2_merged_percent' )
bed_g1_g2$H9_ESC_DMR_STATE <- ""
bed_g1_g2$SM_CELL_DMR_STATE <- ""
bed_g1_g2

#For ESC_DMR
#Just finding the indices for marking the state of DMRs:
hypo_indices  <- which(bed_g1_g2$group1 < 0.25)
bed_g1_g2[hypo_indices,]
intermed_indices  <- which(bed_g1_g2$group1 >= 0.25 & bed_g1_g2$group1 <= 0.75)
bed_g1_g2[intermed_indices,]
hyper_indices  <- which(bed_g1_g2$group1 > 0.75)
bed_g1_g2[hyper_indices,]


#Naming the DMRs:
hypo_indices  <- which(bed_g1_g2$group1 < 0.25)
bed_g1_g2[hypo_indices,]$H9_ESC_DMR_STATE <- "Hypo DMRs"
intermed_indices  <- which(bed_g1_g2$group1 >= 0.25 & bed_g1_g2$group1 <= 0.75)
bed_g1_g2[intermed_indices,]$H9_ESC_DMR_STATE <- "Inter DMRs"
hyper_indices  <- which(bed_g1_g2$group1 > 0.75)
bed_g1_g2[hyper_indices,]$H9_ESC_DMR_STATE <- "Hyper DMRs"
bed_g1_g2[which(bed_g1_g2$H9_ESC_DMR_STATE == ""),]
bed_g1_g2
nrow(bed_g1_g2)


#For SM_cell_state_DMRs:
#Just finding the indices for marking the state of DMRs:
hypo_indices  <- which(bed_g1_g2$group2 < 0.25)
bed_g1_g2[hypo_indices,]
intermed_indices  <- which(bed_g1_g2$group2 >= 0.25 & bed_g1_g2$group2 <= 0.75)
bed_g1_g2[intermed_indices,]
hyper_indices  <- which(bed_g1_g2$group2 > 0.75)
bed_g1_g2[hyper_indices,]


#Naming the DMRs:
hypo_indices  <- which(bed_g1_g2$group2 < 0.25)
bed_g1_g2[hypo_indices,]$SM_CELL_DMR_STATE <- "Hypo DMRs"
intermed_indices  <- which(bed_g1_g2$group2 >= 0.25 & bed_g1_g2$group2 <= 0.75)
bed_g1_g2[intermed_indices,]$SM_CELL_DMR_STATE <- "Inter DMRs"
hyper_indices  <- which(bed_g1_g2$group2 > 0.75)
bed_g1_g2[hyper_indices,]$SM_CELL_DMR_STATE <- "Hyper DMRs"
bed_g1_g2[which(bed_g1_g2$SM_CELL_DMR_STATE == ""),]
bed_g1_g2
nrow(bed_g1_g2)

dir_path="/home/surya/Desktop/ENCODE_poster/parsed_data_R_output"
#Countings of DMRs:
bed_g1_g2
Total_dmr_count <- nrow(bed_g1_g2)
Total_dmr_count

hyper_dmrs <- bed_g1_g2[which(bed_g1_g2$H9_ESC_DMR_STATE == "Hyper DMRs"),]
hyper_dmrs
write.table(hyper_dmrs, paste(dir_path,"hyper_dmrs.txt",sep="/"),row.names=F, col.names= F, quote=FALSE, sep="\t")
hyper_dmr_count <- nrow(hyper_dmrs)
hyper_dmr_count

hypo_dmrs <- bed_g1_g2[which(bed_g1_g2$H9_ESC_DMR_STATE == "Hypo DMRs"),]
hypo_dmrs
write.table(hypo_dmrs, paste(dir_path,"hypo_dmrs.txt",sep="/"),row.names=F, col.names= F, quote=FALSE, sep="\t")
hypo_dmr_count <- nrow(hypo_dmrs)
hypo_dmr_count

inter_dmrs <- bed_g1_g2[which(bed_g1_g2$H9_ESC_DMR_STATE == "Inter DMRs"),]
inter_dmrs
write.table(inter_dmrs, paste(dir_path,"inter_dmrs.txt",sep="/"),row.names=F, col.names= F, quote=FALSE, sep="\t")
inter_dmr_count <- nrow(inter_dmrs)
inter_dmr_count


########
########

DMR_count_df <- data.frame(H9_ESC_DMR_STATE= c("Hyper DMRs", "Inter DMRs", "Hypo DMRs"), DMR_COUNT=c(hyper_dmr_count,inter_dmr_count,hypo_dmr_count ))
DMR_count_df1 <- DMR_count_df %>% arrange(DMR_COUNT)
DMR_count_df1

pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/Encode_poster_fig_1.pdf")

p1 <- ggplot(DMR_count_df1, aes(x="", y=DMR_COUNT, fill=H9_ESC_DMR_STATE)) +
  geom_bar(stat="identity", width = 0.25) + xlab("H9-Human ES Cell DMR States") +ylab("DMR Counts") 

ggplot_build(p1)$data

hypo_color <- "#00BA38"
inter_color <-"#619CFF"
hyper_color <-"#F8766D" 

p1

dev.off()


theme_cfg <- theme(
  #axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )



###### For hyper-transitioning analysis ( hyper to hypo or inter DMRS)#########

bed_g1_g2_hyper_hypo <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Hyper DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Hypo DMRs"]
bed_g1_g2_hyper_hypo
hyper_hypo_count <- nrow(bed_g1_g2_hyper_hypo)

bed_g1_g2_hyper_inter <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Hyper DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Inter DMRs"]
bed_g1_g2_hyper_inter
hyper_inter_count <- nrow(bed_g1_g2_hyper_inter)

bed_g1_g2_hyper_hyper <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Hyper DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Hyper DMRs"]
bed_g1_g2_hyper_hyper
nrow(bed_g1_g2_hyper_hyper)

hyper_transition_df <- data.frame(DMR_transition=c("hyper_to_hypo", "hyper_to_inter"), DMR_COUNT=c(hyper_hypo_count, hyper_inter_count))
hyper_transition_df


pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/Encode_poster_fig_2.pdf")

p2 <- ggplot(hyper_transition_df, aes(x="", y=DMR_COUNT))+
  geom_bar(stat = "identity", width=0.25, fill=c(hypo_color,inter_color))+ xlab("H9-hESC hyper-methylated region transition") +
  ylab("DMR Counts")+ theme(legend.position="none") 
  #+ theme(axis.title.y = element_text(size=14, face="bold")) + 
  #theme(axis.title.x = element_text(size=14, face="bold"))
p2

dev.off()


###### For inter-transitioning analysis ( inter to hypo or inter DMRS)#########

bed_g1_g2_inter_hypo <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Inter DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Hypo DMRs"]
bed_g1_g2_inter_hypo
inter_hypo_count <- nrow(bed_g1_g2_inter_hypo)

bed_g1_g2_inter_inter <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Inter DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Inter DMRs"]
bed_g1_g2_inter_inter %>% head
inter_inter_count <- nrow(bed_g1_g2_inter_inter)

bed_g1_g2_inter_hyper <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Inter DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Hyper DMRs"]
bed_g1_g2_inter_hyper
inter_hyper_count <- nrow(bed_g1_g2_inter_hyper)

inter_transition_df <- data.frame(DMR_transition=c("inter_to_hyper", "inter_to_inter", "inter_to_hypo"), DMR_COUNT=c(inter_hyper_count,inter_inter_count, inter_hypo_count))
inter_transition_df

pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/Encode_poster_fig_3.pdf")

p3 <- ggplot(inter_transition_df, aes(x="", y=DMR_COUNT))+
  geom_bar(stat = "identity", width=0.25, fill=c(hyper_color,inter_color,hypo_color))+ xlab("H9-hESC intermediate-methylated region transition") +
  ylab("DMR Counts")+ theme(legend.position="none") 
  #+ theme(axis.title.y = element_text(size=14, face="bold")) + 
  #theme(axis.title.x = element_text(size=14, face="bold"))
p3

dev.off()



###### For hypo-transitioning analysis ( hypo to hypo or inter DMRS)#########

bed_g1_g2_hypo_hypo <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Hypo DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Hypo DMRs"]
bed_g1_g2_hypo_hypo
hypo_hypo_count <- nrow(bed_g1_g2_hypo_hypo)

bed_g1_g2_hypo_inter <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Hypo DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Inter DMRs"]
bed_g1_g2_hypo_inter
hypo_inter_count <- nrow(bed_g1_g2_hypo_inter)

bed_g1_g2_hypo_hyper <- bed_g1_g2[bed_g1_g2$H9_ESC_DMR_STATE == "Hypo DMRs" & bed_g1_g2$SM_CELL_DMR_STATE == "Hyper DMRs"]
bed_g1_g2_hypo_hyper
hypo_hyper_count <- nrow(bed_g1_g2_hypo_hyper)

hypo_transition_df <- data.frame(DMR_transition=c("hypo_to_hyper", "hypo_to_inter", "hypo_to_hypo"), DMR_COUNT=c(hypo_hyper_count, hypo_inter_count, hypo_hypo_count))
hypo_transition_df


pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/Encode_poster_fig_4.pdf")

p4 <- ggplot(hypo_transition_df, aes(x="", y=DMR_COUNT))+
  geom_bar(stat = "identity", width=0.25, fill=c(hyper_color,inter_color,hypo_color))+ xlab("H9-hESC hypo-methylated region transition") +
  ylab("DMR Counts")+ theme(legend.position="none") 
  #+ theme(axis.title.y = element_text(size=14, face="bold")) + 
  #theme(axis.title.x = element_text(size=14, face="bold"))
p4

dev.off()




#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################





bed_g2_g3 <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/merged_processed_DMR/analyse_DMR_output_data/output_g2_g3_3_sorted_qval.0.05.out", sep="\t")
bed_g2_g3
names(bed_g2_g3) <- c('chr','start','end','q_val','diff','CpG','group1','group2')
bed_g2_g3$SM_CELL_DMR_STATE <- ""
bed_g2_g3$Hepatocytes_DMR_STATE <- ""
bed_g2_g3

#For ESC_DMR
#Just finding the indices for marking the state of DMRs:
hypo_indices  <- which(bed_g2_g3$group1 < 0.25)
bed_g2_g3[hypo_indices,]
intermed_indices  <- which(bed_g2_g3$group1 >= 0.25 & bed_g2_g3$group1 <= 0.75)
bed_g2_g3[intermed_indices,]
hyper_indices  <- which(bed_g2_g3$group1 > 0.75)
bed_g2_g3[hyper_indices,]


#Naming the DMRs:
hypo_indices  <- which(bed_g2_g3$group1 < 0.25)
bed_g2_g3[hypo_indices,]$SM_CELL_DMR_STATE <- "Hypo DMRs"
bed_g2_g3
intermed_indices  <- which(bed_g2_g3$group1 >= 0.25 & bed_g2_g3$group1 <= 0.75)
bed_g2_g3[intermed_indices,]$SM_CELL_DMR_STATE <- "Inter DMRs"
hyper_indices  <- which(bed_g2_g3$group1 > 0.75)
bed_g2_g3[hyper_indices,]$SM_CELL_DMR_STATE <- "Hyper DMRs"
bed_g2_g3


#For SM_cell_state_DMRs:
#Just finding the indices for marking the state of DMRs:
hypo_indices  <- which(bed_g2_g3$group2 < 0.25)
bed_g2_g3[hypo_indices,]
intermed_indices  <- which(bed_g2_g3$group2 >= 0.25 & bed_g2_g3$group2 <= 0.75)
bed_g2_g3[intermed_indices,]
hyper_indices  <- which(bed_g2_g3$group2 > 0.75)
bed_g2_g3[hyper_indices,]


#Naming the DMRs:
hypo_indices  <- which(bed_g2_g3$group2 < 0.25)
bed_g2_g3[hypo_indices,]$Hepatocytes_DMR_STATE <- "Hypo DMRs"
bed_g2_g3
intermed_indices  <- which(bed_g2_g3$group2 >= 0.25 & bed_g2_g3$group2 <= 0.75)
bed_g2_g3[intermed_indices,]$Hepatocytes_DMR_STATE <- "Inter DMRs"
hyper_indices  <- which(bed_g2_g3$group2 > 0.75)
bed_g2_g3[hyper_indices,]$Hepatocytes_DMR_STATE <- "Hyper DMRs"
bed_g2_g3
bed_g2_g3[which(bed_g2_g3$Hepatocytes_DMR_STATE == ""),]
nrow(bed_g2_g3)


theme_cfg <- theme(
  #axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

bed_g2_g3

pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/ENCODE_poster_Piechart_SM_CELL_vs_Hepatocytes_cell_1.pdf")
p1 <- ggplot(bed_g2_g3, aes(x=factor(1), fill=SM_CELL_DMR_STATE))+
  geom_bar(width = 1)+ coord_polar("y") + theme_cfg
p1
dev.off()

pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/ENOCDE_poster_Piechart_SM_CELL_vs_Hepatocytes_cell_2.pdf")
p2 <- ggplot(bed_g2_g3, aes(x=factor(1), fill=Hepatocytes_DMR_STATE))+
  geom_bar(width = 1) + coord_polar("y") + theme_cfg
p2

dev.off()


## Plot statistics
#dir.create("~/Desktop/R_dir")
bed_g2_g3 <- fread("/home/surya/Desktop/scripts/data/metilene_dmr_data/merged_processed_DMR/analyse_DMR_output_data/output_g2_g3_3_sorted_qval.0.05.out", sep="\t")

pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/final_smooth_vs_hepatocytes.pdf")
#difference histogram
p1 <- ggplot(bed_g2_g3, aes(x=diff)) + geom_histogram(binwidth=0.04, fill='darkgrey', color='red') + xlab("Mean Methylation Difference ( SM_Cell - Hepatocytes )") + ylab("DMR count") + scale_x_continuous(limits=c(-1,1)) + theme_cfg
p1
#q_val vs difference
p2 <- ggplot(bed_g2_g3, aes(x=diff, y=q_val)) + geom_point(alpha=.5,color="blue") + scale_y_log10() + xlab("Mean Methylation Difference ( SM cell - Hepatocytes )") + ylab("q-value") + theme_cfg
p2
dev.off()

melt_bed_g1_g2 <- bed_g2_g3
names(melt_bed_g1_g2) <- c('chr','start','end','q_val','diff','CpG','SM_Cell','Hepatocytes')
melt_bed_g1_g2

merged_df <- melt(bed_g2_g3, id.vars = c('chr','start','end','q_val','diff','CpG'))
names(merged_df) <- c('chr','start','end','q_val','diff','CpG', 'Cell_Type', 'value')
merged_df

#cdat <- ddply(merged_df, "Cell_Type", summarise, value.mean=mean(value))
#cdat

mean_line <- merged_df %>% group_by(Cell_Type) %>% summarise(value_mean=mean(value))
mean_line

pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/density_plot_smooth_vs_hepatocytes.pdf")
p3 <-ggplot(merged_df, aes(x=value, y=..count../sum(..count..), group=Cell_Type, colour=Cell_Type, fill=Cell_Type)) +
     geom_density(alpha= 0.3) + geom_vline(data=mean_line, aes(xintercept=value_mean,  colour=Cell_Type),
               linetype="dashed", size=1) + xlab("Mean methylation distribution for DMRs") + ylab("Density") + theme_cfg
p3

dev.off()

pdf(file="~/Desktop/ENCODE_poster/DMR_distribution_plot/Combined_SM_cell_vs_Hepatocytes_cell.pdf")
grid.arrange(p3, arrangeGrob(p1,p2),ncol=2)
dev.off()



###############################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

