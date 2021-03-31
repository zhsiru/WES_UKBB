setwd("C:\\Users\\Sirui\\Desktop\\WORKS\\regeneron")
setwd("C:\\Users\\sirui.zhou\\work\\regeneron-wes")

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(CMplot) 
library(readxl)
library(cowplot)
library(coloc)
library(hrbrthemes)
library(viridis)
library(gapminder)

burden <- read.table("RGC_eBMD__RINT_finemapadj_cohort_UKB_Freeze_300.burden.tsv",header=T)
burden2 <- burden %>% separate(Name, into = paste0('Gene', 1:2), sep = '[()]')
burden2 <- burden2[order(burden2$Pval),]
burden2$MAC <- burden2$Num_Cases-burden2$Cases_Ref
###remove M1.01##
burden2 <- burden2[which(burden2$Alt != "M1.01"), ]
burden2_lowestp <- burden2 %>% 
  group_by(Gene2) %>% 
  filter(Alt == first(Alt))

burden2 %>%
  ggplot( aes(x=AAF, color=Alt, fill=Alt)) +
  geom_density(alpha=0.6) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none"
  ) +
  ylab("") +
  xlab("Alternative Allele Frequency")

burden2_M1 <- burden2[ which(burden2$Alt == "M1.1"), ]
#burden2_M1_1 <- burden2[ which(burden2$Alt == "M1.01"), ]
burden2_M2 <- burden2[ which(burden2$Alt == "M2.1"), ]
burden2_M2_1 <- burden2[ which(burden2$Alt == "M2.01"), ]
burden2_M3 <- burden2[ which(burden2$Alt == "M3.1"), ]
burden2_M3_1 <- burden2[ which(burden2$Alt == "M3.01"), ]
burden2_M4 <- burden2[ which(burden2$Alt == "M4.1"), ]
burden2_M4_1 <- burden2[ which(burden2$Alt == "M4.01"), ]

burden2_M1$P_FDR <- p.adjust(burden2_M1$Pval, method = "BH")
#burden2_M1_1$P_FDR <- p.adjust(burden2_M1_1$Pval, method = "BH")
burden2_M2$P_FDR <- p.adjust(burden2_M2$Pval, method = "BH")
burden2_M2_1$P_FDR <- p.adjust(burden2_M2_1$Pval, method = "BH")
burden2_M3$P_FDR <- p.adjust(burden2_M3$Pval, method = "BH")
burden2_M3_1$P_FDR <- p.adjust(burden2_M3_1$Pval, method = "BH")
burden2_M4$P_FDR <- p.adjust(burden2_M4$Pval, method = "BH")
burden2_M4_1$P_FDR <- p.adjust(burden2_M4_1$Pval, method = "BH")



burden2_M1_sig <- burden2_M1[which(burden2_M1$P_FDR < 0.01), ]
#burden2_M1_1_sig <- burden2_M1_1[which(burden2_M1_1$P_FDR < 0.01), ]
burden2_M2_sig <- burden2_M2[which(burden2_M2$P_FDR < 0.01), ]
burden2_M2_1_sig <- burden2_M2_1[which(burden2_M2_1$P_FDR < 0.01), ]
burden2_M3_sig <- burden2_M3[which(burden2_M3$P_FDR < 0.01), ]
burden2_M3_1_sig <- burden2_M3_1[which(burden2_M3_1$P_FDR < 0.01), ]
burden2_M4_sig <- burden2_M4[which(burden2_M4$P_FDR < 0.01), ]
burden2_M4_1_sig <- burden2_M4_1[which(burden2_M4_1$P_FDR < 0.01), ]




burden2_lowestp$P_FDR <- p.adjust(burden2_lowestp$Pval, method = "BH")
burden_sig_FDR <- burden2_lowestp[which(burden2_lowestp$P_FDR < 0.05), ]

burden2_lowestp$P_BONF <- p.adjust(burden2_lowestp$Pval, method = "bonferroni")
burden_sig_BONF <- burden2_lowestp[which(burden2_lowestp$P_BONF < 0.05), ]



GWAS <- read.table("RGC_eBMD__RINT_finemapadj_cohort_UKB_Freeze_300.single.res.tsv",header=T)

CMplot(burden2_M3[,c(1,3,4,13)], plot.type="m",LOG10=TRUE,col=c("#bdc9e1","#045a8d"),highlight=burden2_M3_sig,
       highlight.col=c("#e34a33"),highlight.cex=2,highlight.pch=c(18), highlight.text=burden2_M3_sig,      
       highlight.text.col=c("#810f7c"),threshold=0.00002,threshold.lty=2,file="jpg",memo="M3_FDR0.01",   
       amplify=FALSE,file.output=T,verbose=TRUE,width=14,height=6)

JM <- read.table("JM_4204_genes.txt",header=T)
new <- read.table("43_new_genes.txt",header=T)

burden_sig_novel <- merge(burden_sig, JM, by = "Gene1")

write.table(burden_sig_2, file="burden_sig_0.05.txt",col.names=T,row.names=F,quote=F,sep="\t")



EI <- read.table("eBMD_EI.txt",header=T)
Fenland <- read.table("Fenland_BMD_protein2.txt",header=T)
Fenland <- Fenland[order(Fenland$pval_ALL),]
Fenland <- Fenland %>% 
  group_by(Gene1) %>% 
  filter(beta_ALL == first(beta_ALL))
pQTL_Sun <- read.table("Sun.Em.cis.IV.indp.final.qc.mr.txt",header=T)
pQTL_Sun <- pQTL_Sun[order(pQTL_Sun$pval_Sun),]
pQTL_Sun <- pQTL_Sun %>% 
  group_by(Gene1) %>% 
  filter(b_Sun == first(b_Sun))
pQTL_Em <- read.table("Em.cis.IV.single.final.qc.mr.txt",header=T)
pQTL_Em <- pQTL_Em[order(pQTL_Em$pval_EM),]
pQTL_Em <- pQTL_Em %>% 
  group_by(Gene1) %>% 
  filter(b_EM == first(b_EM))
COLOC <- read.table("Sun.MR.coloc.txt",header=T)
COLOC <- COLOC[order(COLOC$H4),]
COLOC <- COLOC %>% 
  group_by(Gene1) %>% 
  filter(H4 == last(H4))

eQTL <- read.table("BMD_GTEx_v7_eQTLGen_all.txt",header = T)
eQTL <- eQTL[order(eQTL$p_SMR),]
eQTL_lowestp <- eQTL %>% 
  group_by(Gene2) %>% 
  filter(Tissue == first(Tissue))

eQTLgen <- read.table("bmd_cis-eQTL_eQTLGen.smr",header=T)
eQTLgen <- eQTLgen[order(eQTLgen$p_SMR),]
eQTLgen <- eQTLgen %>% 
  group_by(Gene2) %>% 
  filter(b_SMR == first(b_SMR))


dat <- merge(burden2_lowestp, EI, by="Gene1", all = TRUE)
dat <- merge(dat, Fenland, by = "Gene1", all = TRUE)
dat <- merge(dat, pQTL_Sun, by = "Gene1", all = TRUE)
dat <- merge(dat, pQTL_Em, by = "Gene1", all = TRUE)
dat <- merge(dat, COLOC, by = "Gene1", all = TRUE)
dat <- merge(dat, eQTL_lowestp, by = "Gene2", all = TRUE)
#dat <- merge(dat, eQTLgen, by = "Gene2", all = TRUE)
write.table(dat, file="master_data.txt",col.names=T,row.names=F,quote=F,sep="\t")

###positive control and calculate enrichment###

pc <- unique(read.table("positive_control.txt",header=F))

dat$pc <- match(dat$Gene1, pc[,1], nomatch = 0, incomparables = NULL)

dat$pc=ifelse(dat$pc == "0",0,1)

dat$FDR_sig=ifelse(dat$P_FDR < 0.05,1,0)
dat$BONF_sig=ifelse(dat$P_BONF < 0.05,1,0)

dat$EI_sig=ifelse(dat$EI > 0.75,1,0)
dat$EI_FDR_sig=ifelse(dat$EI > 0.75 & dat$P_FDR < 0.05,1,0)
dat$EI_BONF_sig=ifelse(dat$EI > 0.75 & dat$P_BONF < 0.05,1,0)

dat1_table=table(dat$FDR_sig,dat$pc)
dat2_table=table(dat$BONF_sig,dat$pc)
dat3_table=table(dat$EI_sig,dat$pc)
dat4_table=table(dat$EI_FDR_sig,dat$pc)
dat5_table=table(dat$EI_BONF_sig,dat$pc)

fisher.test(dat1_table)
fisher.test(dat2_table)
fisher.test(dat3_table)
fisher.test(dat4_table)
fisher.test(dat5_table)


#####

CADM1 <- dat[which(dat$Gene1 == "CADM1"), ]
CD109 <- dat[which(dat$Gene1 == "CD109"), ]

dat_WES <- dat[which(dat$Pval < 0.001), ]
dat_WES_eMR <- dat_WES[which(dat_WES$p_SMR < 0.000001), ]
dat_WES_eMR_EI <- dat_WES_eMR[which(dat_WES_eMR$EI > 0.8 ), ]
dat_WES_eMR_EI_coloc <- dat_WES_eMR_EI[which(dat_WES_eMR_EI$p_HEIDI > 0.05), ]
dat_WES_eMR_coloc <- dat_WES_eMR[which(dat_WES_eMR$p_HEIDI > 0.05), ]

dat_WES_Fenland <- dat_WES[which(dat_WES$pval_ALL < 0.001), ]

dat_WES <- dat[which(dat$Pval < 0.05), ]
dat_WES_pMR <- dat_WES[which(dat_WES$pval_Sun < 0.05), ]
dat_WES_pMR_coloc <- unique(dat_WES_pMR[which(dat_WES_pMR$H4 > 0.5), ])

###New scatter

dat_WES <- dat[which(dat$P_FDR < 0.01), ]

dat_WES$P_SMR_FDR <- p.adjust(dat_WES$p_SMR, method = "BH")
dat_WES$pval_ALL_FDR <- p.adjust(dat_WES$pval_ALL, method = "BH")

dat_WES_eMR <- dat_WES[which(dat_WES$P_SMR_FDR < 0.01), ]
dat_WES_eMR_EI <- dat_WES_eMR[which(dat_WES_eMR$EI > 0.8 ), ]

write.table(dat_WES, file="dat_WES_37.txt",col.names=T,row.names=F,quote=F,sep="\t")

###EI vs WES####
EI <- dat %>% filter( EI > 0.6) %>% select (Gene1)
dat$Pval_log <- -log10(dat$Pval)
dat$p_SMR_log <- -log10(dat$p_SMR)
dat$abs_Effect <- abs(dat$Effect)
dat$abs_beta <- abs(dat$b_SMR)
dat2 <- dat[which(dat$Pval < 0.01), ]




# Most basic bubble plot
dat2 %>%
  arrange(desc(abs_Effect)) %>%
  #mutate(country = factor(country, country)) %>%
  ggplot(aes(x=EI, y=Pval_log, size=abs_Effect, 
             fill=abs_beta,
             )) +
  geom_point(alpha=0.9, shape=21, color="black") +
  scale_size(range = c(0.1, 10),
             breaks = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
             #labels = c("250", "500", "750", "1000", "1250"),
             name="absolute effect size for burden") +
  #scale_fill_viridis(discrete=F, option="A") +
  scale_fill_gradientn(colours=c("#bcbddc","#807dba","#3f007d","#810f7c","#e7298a", "#ce1256"),
                       #values=c(0.01,0.05,0.1,0.3,0.5,1),
                       na.value="grey", 
                       guide="colourbar",
                       name="absolute effect size for eMR",limits=c(0,0.3),
                       breaks=c(0,0.05,0.1,0.2,0.3), 
                       labels=c(0,0.05,0.1,0.2,0.3),
                       )+
  #theme_ipsum() +
  #theme(legend.position="bottom") +
  ylab("-log10 pvalue for burden") +
  xlab("EI") 
  #theme(legend.position = "left")


#####pMR#####
EI <- dat_WES_pMR %>% filter( EI > 0.6) %>% select (Gene1)
Fenland <- dat_WES_pMR %>% filter( pval_ALL < 0.001) %>% select (Gene1)	
ggscatter(dat_WES_pMR, x = "b_SMR", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#225ea8", size = 3, alpha = 0.8, 
          title = "Effect Size of pMR and WES Burden (p<0.05)",
          #title = "Effect Size of pMR and WES Burden (p<0.05) in colocalized genes",
          xlab = "pMR Beta in Sun", ylab = "Burdern Beta of gene with lowest P", 
          label.select = Fenland$Gene1,
          label="Gene1", 
          font.label = c(11, "bold", "#045a8d"), 
          repel = TRUE) +
  stat_cor(method = "pearson")

#####eMR#####

EI <- dat_WES_eMR %>% filter( EI > 0.8) %>% select (Gene1)
Fenland <- dat_WES_eMR %>% filter( pval_ALL_FDR < 0.01) %>% select (Gene1)
pQTL <- dat_WES_eMR %>% filter( pval_Sun < 0.01) %>% select (Gene1)


p1 <- ggscatter(dat_WES_eMR, x = "b_SMR", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#fb6a4a", size = 3, alpha = 0.8, 
          title = "A. FDR < 0.01 (Hightlight EI)",
          #title = "Effect Size of eMR (p<1e-8) and WES Burden (p<5e-4) in colocalized genes",
          xlab = "eMR Beta of gene with lowest P", ylab = "Burdern Beta of gene with lowest P",
          label.select = EI$Gene1,
          label="Gene1", 
          font.label = c(11, "bold", "#045a8d"), 
          repel = TRUE) + stat_cor(method = "pearson", label.x = -0.15, label.y = -0.28)

p2 <- ggscatter(dat_WES_eMR_EI, x = "b_SMR", y = "Effect",
                add = "reg.line",
                conf.int = TRUE, color = "#045a8d", size = 3, alpha = 0.8, 
                title = "B. FDR < 0.01 and EI > 0.8",
                #title = "Effect Size of eMR (p<1e-8) and WES Burden (p<5e-4) in colocalized genes",
                xlab = "eMR Beta of gene with lowest P", ylab = "Burdern Beta of gene with lowest P",
                label="Gene1", 
                font.label = c(11, "bold", "#c51b8a"), 
                repel = TRUE) + stat_cor(method = "pearson", label.x = -0.15, label.y = -0.28)

p3 <- ggscatter(dat_WES_eMR, x = "b_SMR", y = "Effect",
                add = "reg.line",
                conf.int = TRUE, color = "#fb6a4a", size = 3, alpha = 0.8, 
                title = "C. FDR < 0.01 (Highlight Fenland)",
                #title = "Effect Size of eMR (p<1e-8) and WES Burden (p<5e-4) in colocalized genes",
                xlab = "eMR Beta of gene with lowest P", ylab = "Burdern Beta of gene with lowest P",
                label="Gene1", 
                label.select = Fenland$Gene1,
                font.label = c(11, "bold", "#005a32"), 
                repel = TRUE) + stat_cor(method = "pearson", label.x = -0.15, label.y = -0.28)

p4 <- ggscatter(dat_WES_eMR, x = "b_SMR", y = "Effect",
                add = "reg.line",
                conf.int = TRUE, color = "#fb6a4a", size = 3, alpha = 0.8, 
                title = "D. FDR < 0.01 (Highlight pQTL)",
                #title = "Effect Size of eMR (p<1e-8) and WES Burden (p<5e-4) in colocalized genes",
                xlab = "eMR Beta of gene with lowest P", ylab = "Burdern Beta of gene with lowest P",
                label="Gene1", 
                label.select = pQTL$Gene1,
                font.label = c(11, "bold", "#4a1486"), 
                repel = TRUE) + stat_cor(method = "pearson", label.x = -0.15, label.y = -0.28)

grid.arrange(p1, p2, p3, p4, ncol =2)

#####eMR+EI#####
ggscatter(dat_WES_eMR_EI_coloc, x = "b_SMR", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#ff7f00", size = 3, alpha = 0.8, 
          #title = "Effect Size of eMR (p<1e-6) and WES Burden (p<1e-3) in EI genes (>0.8)",
          title = "Effect Size of eMR (p<1e-6) and WES Burden (p<1e-3) in colocalized and EI genes (>0.8)",
          xlab = "eMR Beta of gene with lowest P", ylab = "Burdern Beta of gene with lowest P", 
          label="Gene1", 
          font.label = c(11, "bold", "#3f007d"), 
          repel = TRUE) +
  stat_cor(method = "pearson")




#####Fenland#####
ggscatter(dat_WES_Fenland, x = "beta_ALL", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#ff7f00", size = 3, alpha = 0.8, 
          #title = "Effect Size of eMR (p<1e-6) and WES Burden (p<1e-3) in EI genes (>0.8)",
          title = "Effect Size of Fenland and WES Burden (p<1e-3) genes",
          xlab = "Fenland tBMD beta", ylab = "Burdern Beta of gene with lowest P", 
          label="Gene1", 
          font.label = c(11, "bold", "#3f007d"), 
          repel = TRUE) +
  stat_cor(method = "pearson")


strict <- dat[which(dat$Pval < 0.000003), ]
strict <- strict[which(strict$p_SMR < 8.4e-6), ]
strict <- strict[which(strict$EI > 0.8 ), ]
strict <- strict[which(strict$pval_ALL < 0.001), ]


dat %>% select(1, 3)


old <- read.table("WES_old.txt",header=T)

old <- old[order(old$Pval_old),]

old <- old %>% 
  group_by(Gene1) %>% 
  filter(Alt_old == first(Alt_old))

test <- merge(dat, old, by="Gene2")
test$logp_old <- -log10(test$Pval_old)
test$logp_new <- -log10(test$Pval)
test_2 <- test[which(test$Pval < 0.05 & test$Pval_old < 0.05), ]
test_3 <- test[which(test$Pval < 0.000003), ]
test_4 <- test[which(test$Pval_old < 0.000003), ]
sig <- unique(rbind(test_3, test_4))
CADM1 <- test[which(test$Gene1.x == "CADM1"), ]

###correlation between old and new WES

ggscatter(sig, x = "Effect_old", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#984ea3", size = 3, alpha = 0.8, 
          title = "Effect Size of Burden new and old in GW significant genes",
          xlab = "Burdern Beta of gene with lowest P (Old)", ylab = "Burdern Beta of gene with lowest P (New)", 
          label="Gene1.x", 
          font.label = c(10, "bold", "#3f007d"), 
          repel = TRUE) +
  stat_cor(method = "pearson", label.x = -0.05, label.y = -0.25)


ggscatter(test_4, x = "logp_old", y = "logp_new",
          add = "reg.line",
          conf.int = TRUE, color = "#33a02c", size = 3, alpha = 0.8, 
          title = "Pvalue of significant genes in old data",
          xlab = "-log 10 p value old", ylab = "-log 10 p value new", 
          label="Gene1.x", 
          font.label = c(10, "bold", "#3f007d"), 
          repel = TRUE) +
  stat_cor(method = "pearson")



MR_xWES_lowestp_sig1 <- MR_xWES_lowestp[which(MR_xWES_lowestp$Pval < 0.005), ]
MR_xWES_lowestp_sig2 <- MR_xWES_lowestp_sig1[which(MR_xWES_lowestp_sig1$p_SMR < 0.005), ]
MR_xWES_lowestp_sig2_Heidi <- MR_xWES_lowestp_sig2[which(MR_xWES_lowestp_sig2$p_HEIDI > 0.05), ]

ggscatter(MR_xWES_lowestp_sig2, x = "b_SMR", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#225ea8", size = 3, alpha = 0.8, 
          title = "Effect Size of eQTLgen MR and WES Burden (p<0.05)",
          xlab = "eQTLgen MR Beta", ylab = "Burdern Beta of gene with lowest P", 
          repel = TRUE) +
  stat_cor(method = "pearson")

ggscatter(MR_xWES_lowestp_sig2_Heidi, x = "b_SMR", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#7a0177", size = 3, alpha = 0.8, 
          title = "eQTLgen MR/Heidi and WES Burden (p<0.05)",
          xlab = "eQTLgen MR (Heidi colocalized) Beta", ylab = "Burdern Beta of gene with lowest P", 
          #label="GENE1", 
          #font.label = c(14, "bold", "#3f007d"), 
          repel = TRUE) +
  stat_cor(method = "pearson")



source("plot_SMR.r") 
SMRData = ReadSMRData("CADM1_plot_ENSG_GTEx.ENSG00000182985.12.txt")
SMRLocusPlot(data=SMRData, smr_thresh=8.4e-6, heidi_thresh=0.05, plotWindow=500, max_anno_probe=1)
SMREffectPlot(data=SMRData) 


write.table(MR_xWES_lowestp, file="MR_xWES_lowestp.txt",col.names=T,row.names=F,quote=F,sep="\t")

write.table(MR_xWES_lowestp_sig1_coloc, file="MR_xWES_lowestp_sig_coloc.txt",col.names=T,row.names=F,quote=F,sep="\t")

MR_xWES_lowestp_sig1 <- MR_xWES_lowestp[which(MR_xWES_lowestp$Pval < 0.5), ]

MR_xWES_lowestp_sig1_Sun <- MR_xWES_lowestp_sig1[which(MR_xWES_lowestp_sig1$MR.Pval.Sun < 0.05), ]

MR_xWES_lowestp_sig1_coloc <- MR_xWES_lowestp_sig1_Sun[which(MR_xWES_lowestp_sig1_Sun$COLOC.H4 > 0.5), ]


ggscatter(MR_xWES_lowestp_sig1_Sun, x = "MR.beta.Sun", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#756bb1", size = 3, alpha = 0.6, 
          title = "Effect Size of MR (Sun et al) and WES Burden Test Significant Genes",
          xlab = "MR Beta", ylab = "Burdern Beta of gene with lowest P", 
          repel = TRUE) +
  stat_cor(method = "pearson", label.x = -0.1, label.y = -0.5)

ggscatter(MR_xWES_lowestp_sig1_coloc, x = "MR.beta.Sun", y = "Effect",
          add = "reg.line",
          conf.int = TRUE, color = "#756bb1", size = 3, alpha = 0.6, 
          title = "Effect Size of MR Colocalized and WES Burden Test Significant Genes",
          xlab = "MR Beta", ylab = "Burdern Beta of gene with lowest P", label="GENE1", 
          font.label = c(14, "bold", "#f03b20"), 
          repel = TRUE) +
  stat_cor(method = "pearson", label.x = -0.05, label.y = -0.25)




library(ggplot2)
library(grid)
library(gridExtra)
require(scales)

test=data.frame(read.table("Coloc.txt",header=T,sep = "\t"))
p = ggplot(data=test,
           aes(x = Group, y = Beta, ymin = LCI, ymax = UCI ), width=0.1, cex=1.2)+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group), yintercept =0, linetype=2)+
  xlab('')+ ylab("Beta (95% Confidence Interval)") +
  labs(color='Source') + 
  geom_errorbar(aes(ymin=LCI, ymax=UCI,col=Group),width=0.1,cex=1.2)+ 
  facet_wrap(~GENE1,strip.position="left",nrow=4,scales = "free_y") +
  theme_classic() + 
  theme(plot.title=element_text(size=14,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"), strip.text.y = element_text(size=10, face="bold")) +
  scale_y_continuous(breaks=seq(-0.5,0.19,0.08))+
  #scale_y_continuous(trans = log10_trans(),
  #                  breaks = trans_breaks("log10", function(x) 10^x, n=8),
  #                 labels = trans_format("log10", math_format(10^.x))) +
  coord_flip()
p + scale_color_manual(values=c('#1c9099','#e34a33','#31a354')) + theme(legend.position = c(0.1, 0.2)) + guides(colour = guide_legend(reverse=TRUE))





MR_coloc <- read.table("MR.coloc.txt",header=T)
MR <- read.table("MR.txt",header=T)

MR_coloc <- MR_coloc[order(MR_coloc$Pval),]

MR_coloc_2 <- MR_coloc %>% 
  group_by(GENE1) %>% 
  filter(Effect == first(Effect))

MR_coloc_3 <- unique(MR_coloc_2[, c("GENE1", "P", "Beta", "Effect", "Pval")])

MR <- MR[order(MR$Pval),]

MR_2 <- MR %>% 
  group_by(GENE1) %>% 
  filter(Effect == first(Effect))

MR_3 <- unique(MR_2[, c("GENE1", "P", "Beta", "Effect", "Pval")])



ggscatter(MR_coloc_3, x = "Beta", y = "Effect",
          add = "reg.line",
          conf.int = FALSE, color = "#756bb1", size = 3, alpha = 0.6, 
          cor.coef = T, cor.method = "pearson",
          title = "Effect Size of MR Colocalized and WES Burden Test Significant Genes",
          xlab = "MR Beta", ylab = "Burdern Beta of gene with lowest P", label="GENE1", font.label = c(14, "bold", "#f03b20"), repel = TRUE)


ggscatter(MR_3, x = "Beta", y = "Effect",
          add = "reg.line",
          conf.int = FALSE, color = "#756bb1", size = 3, alpha = 0.6, 
          cor.coef = T, cor.method = "pearson",
          title = "Effect Size of MR and WES Burden Test Significant Genes",
          xlab = "MR Beta", ylab = "Burdern Beta of gene with lowest P", repel = TRUE)


###coloc using susie###

CADM1_E <- read.table("CADM1_Thyroid.coloc",header=T)
CADM1_O <- read.table("CADM1_eBMD.coloc",header=T)
CADM1_E$eaf <- as.numeric(CADM1_E$eaf)
CADM1_O$eaf <- as.numeric(CADM1_O$eaf)
CADM1_E$snp <- as.character(CADM1_E$snp)
CADM1_O$snp <- as.character(CADM1_O$snp)
CADM1_E$maf=ifelse(CADM1_E$eaf>0.5,1-CADM1_E$eaf,CADM1_E$eaf)
CADM1_O$maf=ifelse(CADM1_O$eaf>0.5,1-CADM1_O$eaf,CADM1_O$eaf)
CADM1_E <- CADM1_E[order(CADM1_E$p),]
CADM1_E <- CADM1_E %>% 
  group_by(snp) %>% 
  filter(p == first(p))
CADM1_O <- CADM1_O[order(CADM1_O$p),]
CADM1_O <- CADM1_O %>% 
  group_by(snp) %>% 
  filter(p == first(p))

SNP <- merge(CADM1_E, CADM1_O, by="snp")
write.table(SNP$snp,file="SNP.txt",col.names=F,row.names=F,quote=F,sep=" ") 

coloc.res <- coloc.abf(dataset1=list(pvalues=CADM1_E$p, 
                                     beta=CADM1_E$beta, 
                                     varbeta=CADM1_E$varbeta, 
                                     snp=CADM1_E$snp, 
                                     MAF=CADM1_E$maf, 
                                     N=714,
                                     type="quant"), 
                       dataset2=list(beta=CADM1_O$beta, 
                                     varbeta=CADM1_O$varbeta, 
                                     snp=CADM1_O$snp, 
                                     MAF=CADM1_O$maf, 
                                     N=426824,  
                                     type="quant"), 
                       p1 = 1e-04, 
                       p2 = 1e-04, 
                       p12 = 1e-05)

LD <- read.table("CADM1_eur.ld")
colnames(LD) <- SNP$snp
rownames(LD) <- SNP$snp
CADM1_E <- SNP[c(1:6)]
CADM1_O <- SNP[c(1,7:11)]

DE <- list(pvalues=CADM1_E$p.x, 
           beta=CADM1_E$beta.x, 
           varbeta=CADM1_E$varbeta.x, 
           snp=CADM1_E$snp, 
           MAF=CADM1_E$maf.x, 
           N=714,
           type="quant")

DO <- list(beta=CADM1_O$beta.y, 
           varbeta=CADM1_O$varbeta.y, 
           snp=CADM1_O$snp, 
           MAF=CADM1_O$maf.y, 
           N=426824,  
           type="quant")

str(DO)

DE$LD <- as.matrix(LD)
DO$LD <- as.matrix(LD)

pos <- read.table("CADM1_1kg_eur.bim")
pos <- pos[c(2,4)]
colnames(pos) <- c("snp","position")

DE$position <- pos$position
DO$position <- pos$position

SO<- runsusie(DO,nref=503,check_R = FALSE)
SE<- runsusie(DE,nref=503,check_R = FALSE)

susie.res=coloc.susie(SO,SE)

print(susie.res$summary)
sensitivity(susie.res,"H4 > 0.9",row=6
            ,dataset1=DO,dataset2=DE)


