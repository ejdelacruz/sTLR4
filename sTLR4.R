#'---
#' author: "Erin dela Cruz"
#' -----------------------------------------------------------------------------------------------
#' 
#' **Program**: sTLR4.R
#' 
#' **Purpose**: Analyze results from TLR4 CLIA and correlation to [BVAB]
#' 
#' **Last Edited**: Jan 17 2020
#' 
#' -----------------------------------------------------------------------------------------------
#' 
#+ setup-libs, include=FALSE
#' **Prep**: Load necessary libraries and clear data from memory
rm(list = ls())
# Tidyverse Libs
library(tidyverse)
# For epidemiologic analysis
library(epitools)
# ggplot add-ons
library(ggthemes)
library(cowplot)

#' For use on PC
fp <- "H:/code/sTLR4_pub/"

#' Load sTLR4 CLIA measurements, exclude age > 50
tlr4.conc <- read.csv(paste0(fp, "indata/sTLR4_samples.csv")) %>%
  filter(age<50) %>%
  mutate(mc.phase = ifelse(days.since.lmp<14, "Proliferative", NA)) %>%
  mutate(mc.phase = ifelse(days.since.lmp>14, "Secretory", mc.phase)) %>%
  mutate(hc.mens = mc.phase) %>%
  mutate(hc.mens = ifelse(contraception_horm==1, "Hormonal Contraception", hc.mens))
tlr4.conc$hc.mens <- factor(tlr4.conc$hc.mens, c("Proliferative", "Secretory", "Hormonal Contraception"))

#' Subsets
bvn.nohc <- tlr4.conc %>% filter(nugentdx==0 & contraception_horm==0)
bvp.nohc <- tlr4.conc %>% filter(nugentdx==1 & contraception_horm==0)
bvn <- tlr4.conc %>% filter(nugentdx==0)
bvp <- tlr4.conc %>% filter(nugentdx==1)
bvp.nohc <- tlr4.conc %>% filter(nugentdx==1)
nohc <- tlr4.conc %>% filter(contraception_horm==0)

#' Graphs
theme_set(theme_tufte())
#' --------------------------------------------------------------------
#' Distribution of [sTLR4] in endocervical cytobrush supernatant
#' LOD=0.0625 ng/mL
#' --------------------------------------------------------------------
summary.tlr4.stats <- tlr4.conc %>%
  summarise(avg.tlr4 = mean(tlr4.ngmL), median.tlr4 = median(tlr4.ngmL), sd.tlr4 = sd(tlr4.ngmL), N=n())

#' --------------------------------------------------------------------
#' Testing patient characteristics associated with [sTLR4] in endocervical
#' cytobrush supernatant: (1) phase of menstrual cycle, (2) age, 
#' (3) hormonal contraeptive use, (4) BV, microbiological or clinical, 
#' (5) TLR4 deficiency (by SNP), and (6) endocervical ectopy
#' --------------------------------------------------------------------

age.cat <- ggplot(tlr4.conc, aes(x = factor(age_cat), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  scale_y_continuous("[sTLR4] (ng/mL)") +
  scale_x_discrete("Age Category (years)", labels=c("18-20", "21-30", "31-40", "41-49"))
kruskal.test(tlr4.conc$tlr4.ngmL~tlr4.conc$age_cat)

cyclephase <- ggplot(subset(nohc, is.na(mc.phase)==FALSE), aes(x = factor(mc.phase), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  labs(x="Menstrual Cycle Phase") + 
  scale_y_continuous("[sTLR4] (ng/mL)") +
  labs(subtitle="No hormonal contraception")
wilcox.test(nohc$tlr4.ngmL~nohc$mc.phase)

hc <- ggplot(subset(tlr4.conc, is.na(hc.mens)==FALSE), aes(x = factor(contraception_horm, ordered=FALSE), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  labs(x="Hormonal Contraception") + 
  scale_y_continuous("[sTLR4] (ng/mL)") +
  scale_x_discrete(labels=c("No", "Yes"))
wilcox.test(tlr4.conc$tlr4.ngmL~tlr4.conc$contraception_horm)


cyclephase.hc <- ggplot(subset(tlr4.conc, is.na(hc.mens)==FALSE), aes(x = factor(hc.mens), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  labs(x="Menstrual Cycle Phase") + 
  scale_y_continuous("[sTLR4] (ng/mL)")
kruskal.test(tlr4.conc$tlr4.ngmL~tlr4.conc$hc.mens)

hc.prolif <- tlr4.conc %>%
  filter(hc.mens!="Secretory" & is.na(hc.mens)==FALSE)
wilcox.test(hc.prolif$tlr4.ngmL~hc.prolif$hc.mens)
hc.secretory <- tlr4.conc %>%
  filter(hc.mens!="Proliferative" & is.na(hc.mens)==FALSE)
wilcox.test(hc.secretory$tlr4.ngmL~hc.secretory$hc.mens)
prolif.secretory <- tlr4.conc %>%
  filter(hc.mens!="Hormonal Contraception" & is.na(hc.mens)==FALSE)
wilcox.test(prolif.secretory$tlr4.ngmL~prolif.secretory$hc.mens)

nugent <- ggplot(tlr4.conc, aes(x = factor(nugentdx), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  labs(x="Microbiological (Nugent) Diagnosis") + 
  scale_y_continuous("[sTLR4] (ng/mL)") +
  scale_x_discrete(labels=c("BV-", "BV+"))
wilcox.test(tlr4.conc$tlr4.ngmL~tlr4.conc$nugentdx)

amsel <- ggplot(tlr4.conc, aes(x = factor(amsel), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  labs(x="Clinical (Amsel) Diagnosis") + 
  scale_y_continuous("[sTLR4] (ng/mL)") +
  scale_x_discrete(labels=c("BV-", "BV+"))
wilcox.test(tlr4.conc$tlr4.ngmL~tlr4.conc$amsel)

tlr4.conc$nugentdx <- factor(tlr4.conc$nugentdx, labels = c("BV-negative", "BV-positive"))
tlr4 <- ggplot(tlr4.conc, aes(x = factor(TLR4_896deficient), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  facet_grid(~nugentdx) +
  labs(x="TLR4 rs4986790") + 
  scale_y_continuous("[sTLR4] (ng/mL)") +
  scale_x_discrete(labels=c("AA", "AG/GG"))
wilcox.test(tlr4.conc$tlr4.ngmL~tlr4.conc$TLR4_896deficient)
wilcox.test(bvn$tlr4.ngmL~bvn$TLR4_896deficient)
wilcox.test(bvp$tlr4.ngmL~bvp$TLR4_896deficient)


ectopy <- ggplot(tlr4.conc, aes(x = factor(ect.any), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  scale_y_continuous("[sTLR4] (ng/mL)") +
  scale_x_discrete("Endocervical Ectopy", labels=c("No Ectopy", "Ectopy >25%"))
wilcox.test(tlr4.conc$tlr4.ngmL, tlr4.conc$ect.any)

ect.allcats <- ggplot(tlr4.conc, aes(x = factor(ECTOPY), y = tlr4.ngmL)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.75) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  scale_y_continuous("[sTLR4] (ng/mL)") +
  scale_x_discrete("Estimated Extent of Ectopy", labels = c("0%", "25%", "50%", "75%", "100%"))

#'----------------------------------------------------------------
#' Correlation between sTLR4 and [G-negative BVAB]
#'----------------------------------------------------------------
#' Load qPCR
qpcr.long <- read.csv(paste0(fp, "indata/qPCR_data.csv")) %>%
  right_join(., tlr4.conc, by="patientID")
#' Replace with full assay name
qpcr.long$assay <- as.character(qpcr.long$assay)
qpcr.long$assay[qpcr.long$assay=="leptosneath"] <- "Sneathia"
qpcr.long$assay[qpcr.long$assay=="pamnii"] <- "P. amnii"
qpcr.long$assay[qpcr.long$assay=="ptimon"] <- "P. timonensis"
#' Create separate datasets for each qPCR assay
leptosneath <- qpcr.long[which(qpcr.long$assay == "Sneathia"),]
pamnii <- qpcr.long[which(qpcr.long$assay == "P. amnii"),]
ptimon <- qpcr.long[which(qpcr.long$assay == "P. timonensis"),]

g.pamnii <- ggplot(subset(pamnii, copies_per_swab>lod), aes(x=copies_per_swab, y=tlr4.ngmL)) +
  scale_x_log10("Prevotella amnii\n16S rRNA Copies Per Swab", breaks=c(1e3, 1e5, 1e7, 1e9, 1e11)) +
  scale_y_sqrt("[sTLR4] (ng/mL)") + 
  geom_point()

g.ptimon <- ggplot(subset(ptimon, copies_per_swab>lod), aes(x=copies_per_swab, y=tlr4.ngmL)) +
  scale_x_log10("Prevotella timonensis\n16S rRNA Copies Per Swab", breaks=c(1e3, 1e5, 1e7, 1e9, 1e11)) +
  scale_y_sqrt("[sTLR4] (ng/mL)") + 
  geom_point()

g.sneathia <- ggplot(subset(leptosneath, copies_per_swab>lod), aes(x=copies_per_swab, y=tlr4.ngmL)) +
  scale_x_log10("Sneathia spp.\n16S rRNA Copies Per Swab", breaks=c(1e3, 1e5, 1e7, 1e9, 1e11)) +
  scale_y_sqrt("[sTLR4] (ng/mL)") + 
  geom_point()



#'----------------------------------------------------------------
#'Measurements on primary cell supernatants and longitudinal sampling
#'on single enrollee
#'----------------------------------------------------------------
d18.02.02 <- read.csv(paste0(fp,"indata/2018-02-02_TLR4CLIA.csv"),
                         header=TRUE, check.names=TRUE, as.is=TRUE,
                         stringsAsFactors = FALSE)
pri.cell <- d18.02.02 %>% filter(Sample.Name != "4241" & Sample.Name != "Saline") %>%
  rename(studyID=Sample.Name, filt.ind=Day, rep=collection.date) %>%
  mutate(filt.ind = ifelse(studyID=="Media", "Unfiltered", filt.ind)) %>%
  group_by(studyID, filt.ind, rep) %>%
  mutate(tlr4.ng=tlr4.pgmL/1e3)
pri.cell.sups <- ggplot(pri.cell, aes(x = factor(studyID), y = tlr4.ng, shape=filt.ind)) +
  geom_jitter(width=0.1, size = 3) +
  labs(x="", caption = "* Below level of detection (0.0625 ng/mL)") + 
  scale_y_continuous("[sTLR4] (ng/mL)") +
  scale_x_discrete(labels=c("Fm Media*", "HCE2 Supernatant")) +
  theme(legend.title=element_blank(),
        legend.position = "right")
#' Summary statistics for primary cells
sups <- d18.02.02 %>% filter(Sample.Name!="4241" & Sample.Name!="Saline") %>%
  rename(filt = Day, batch = collection.date) %>%
  mutate(tlr4.ngmL=tlr4.pgmL/1e3) %>%
  group_by(filt, batch) %>%
  summarise(mean.ngmL = mean(tlr4.ngmL), sd = sd(tlr4.ngmL), N=n(), se=sd/sqrt(N))

bvr01.pt1a <- d18.02.02 %>% filter(Sample.Name=="pt1a") %>%
  rename(studyID=Sample.Name, days.since.lmp=Day) %>%
  mutate(studyID = as.numeric(studyID)) %>%
  mutate(days.since.lmp = as.numeric(days.since.lmp)) %>%
  mutate(collection.date = as.Date(collection.date, "%m/%d/%y")) %>%
  group_by(studyID, days.since.lmp, collection.date) %>%
  mutate(tlr4.ng=tlr4.pgmL/1e3) %>%
  summarise(tlr4.ngmL = mean(tlr4.ng))

g.pt1a <- ggplot(bvr01.4241, aes(x=days.since.lmp, y=tlr4.ngmL)) +
  geom_rect(aes(xmin=13, xmax=15, ymin=0, ymax=Inf), alpha=0.5, fill="grey75") +
  geom_point() + geom_line() +
  scale_y_continuous("Vaginal swab supernatant\n[sTLR4] (ng/mL)") +
  scale_x_continuous("Days since last menstrual period")

#'----------------------------------------------------------------
#' AMASS FIGURES FOR PUBLICATION
#'----------------------------------------------------------------

#' Main Figure
main.fig <- plot_grid(pri.cell.sups, nugent, tlr4, 
                      ectopy, cyclephase, g.pt1a, nrow=2, ncol=3, labels="AUTO")

#' Supplemental figure
fig.s1 <- plot_grid(amsel, age.cat, NULL,
                    g.pamnii, g.ptimon, g.sneathia, 
                    ect.allcats, hc, cyclephase.hc,
                    labels=c("A", "B", "", "C", "D","E", "F", "G", "H"), nrow=3, ncol=3)