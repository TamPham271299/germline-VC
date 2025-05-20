library(tidyr)
library(dplyr)

# SNV
setwd("/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/WES_Korean/vep")

cons_type <- c("missense", "LoF", "syno", "utr")

for (c in cons_type) {
    print(c)
    df <- read.table(paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".tsv"), header = T, sep = "\t", na.strings=".")
    df <- df %>% 
        separate_rows(ALT, AC, AF, sep=",") %>%
        filter(ALT==Allele)
    df <- as.data.frame(df)
    df$AF <- as.numeric(df$AF)
    df$AC <- as.numeric(df$AC)
    df <- df[df$AF <= 0.05,]
    df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)

    ## SVA
    df.sva <- df
    df.sva$RC <- df.sva$AN - df.sva$AC
    df.sva$gnomADe_RC <- df.sva$gnomADe_AN - df.sva$gnomADe_AC
    df.sva$gnomADe_RC_eas <- df.sva$gnomADe_AN_eas - df.sva$gnomADe_AC_eas

    ### EAS gnomeADe
    t.gnomADe.eas <- apply(df.sva[, c("AC", "RC", "gnomADe_AC_eas", "gnomADe_RC_eas")], 1, function(x) {
        t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
        return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
    })
    t.gnomADe.eas <- as.data.frame(t(t.gnomADe.eas))
    colnames(t.gnomADe.eas)[2] <- "OR"
    adj.pVal.gnomADe.eas <- p.adjust(t.gnomADe.eas$pVal, method = "bonferroni")
    EF.gnomADe.eas <- (t.gnomADe.eas$OR-1)/t.gnomADe.eas$OR

    ### ALL gnomeADe
    t.gnomADe <- apply(df.sva[, c("AC", "RC", "gnomADe_AC", "gnomADe_RC")], 1, function(x) {
        t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
        return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
    })
    t.gnomADe <- as.data.frame(t(t.gnomADe))
    colnames(t.gnomADe)[2] <- "OR"
    adj.pVal.gnomADe <- p.adjust(t.gnomADe$pVal, method = "bonferroni")
    EF.gnomADe <- (t.gnomADe$OR-1)/t.gnomADe$OR

    df.res <- cbind(df.sva,
                    p_value_EA_gnomADe=t.gnomADe.eas$pVal, 
                    FDR_EA_gnomADe=adj.pVal.gnomADe.eas, 
                    OR_EA_gnomADe=t.gnomADe.eas$OR,
                    CI_lower_EA_gnomADe=t.gnomADe.eas$CI_lower,
                    CI_upper_EA_gnomADe=t.gnomADe.eas$CI_upper,
                    EF_EA_gnomADe=EF.gnomADe.eas,
                    p_value_ALL_gnomADe=t.gnomADe$pVal, 
                    FDR_ALL_gnomADe=adj.pVal.gnomADe,
                    OR_ALL_gnomADe=t.gnomADe$OR,
                    CI_lower_ALL_gnomADe=t.gnomADe$CI_lower,
                    CI_upper_ALL_gnomADe=t.gnomADe$CI_upper,
                    EF_ALL_gnomADe=EF.gnomADe,
                    IS.SIG=rep("ns", nrow(df.sva)))

    df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe > 0.05] <- "EA"
    df.res$IS.SIG[df.res$OR_ALL_gnomADe > 1 & df.res$FDR_ALL_gnomADe < 0.05 & df.res$FDR_EA_gnomADe > 0.05] <- "ALL"
    df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$OR_ALL_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe < 0.05] <- "EA/ALL"

    df.out <- cbind(df.res[, c("Pos", "SYMBOL", "Gene", "AC", "AN", "AF")],
                    `AC - gnomADe (EA)`=df.res$gnomADe_AC_eas, `AN - gnomADe (EA)`=df.res$gnomADe_AN_eas, `AF - gnomADe (EA)`=sprintf("%.2e", df.res[, c("gnomADe_AF_eas")]),
                    `OR (CI) - gnomADe (EA)`=paste0(round(df.res$OR_EA_gnomADe,2), " (", round(df.res$CI_lower_EA_gnomADe,2), "-", round(df.res$CI_upper_EA_gnomADe,2), ")"),
                    `p-value - gnomADe (EA)`=sprintf("%.2e", df.res$p_value_EA_gnomADe), `FDR - gnomADe (EA)`=sprintf("%.2e", df.res$FDR_EA_gnomADe), `EF - gnomADe (EA)`=round(df.res$EF_EA_gnomADe,2),
                    `AC - gnomADe (ALL)`=df.res$gnomADe_AC, `AN - gnomADe (ALL)`=df.res$gnomADe_AN, `AF - gnomADe (ALL)`=sprintf("%.2e", df.res[, c("gnomADe_AF")]),
                    `OR (CI) - gnomADe (ALL)`=paste0(round(df.res$OR_ALL_gnomADe,2), " (", round(df.res$CI_lower_ALL_gnomADe,2), "-", round(df.res$CI_upper_ALL_gnomADe,2), ")"),
                    `p-value - gnomADe (ALL)`=sprintf("%.2e", df.res$p_value_ALL_gnomADe), `FDR - gnomADe (ALL)`=sprintf("%.2e", df.res$FDR_ALL_gnomADe), `EF - gnomADe (ALL)`=round(df.res$EF_ALL_gnomADe,2),
                    is.sig=df.res$IS.SIG,
                    df.res[, c("Consequence", "EXON", "HGVSc", "HGVSp", "SIFT", "PolyPhen", "CADD_RAW", "CADD_PHRED", "CLIN_SIG", "gnomADe", "VAR_SYNONYMS")],
                    df.res[, c(60:75)])

    df.out <- df.out[order(factor(df.out$is.sig, levels=c("EA/ALL", "EA", "ALL", "ns"))),]
    
    write.table(df.res, paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".SVA.tsv"), sep="\t", quote=F, row.names=F)
    write.csv(df.out, paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".SVA.csv"))

    ## GBA
    df.gba <- df
    df.gba <- df.gba %>%
        group_by(Gene, SYMBOL) %>%
        summarise(No.of.variants=n(),
                    Pos=paste(Pos, collapse=";"),
                    Consequence=paste(Consequence, collapse=";"),
                    EXON=paste(EXON, collapse=";"),
                    HGVSc=paste(HGVSc, collapse=";"),
                    HGVSp=paste(HGVSp, collapse=";"),
                    SIFT=paste(SIFT, collapse=";"),
                    PolyPhen=paste(PolyPhen, collapse=";"),
                    CADD_RAW=paste(CADD_RAW, collapse=";"),
                    CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                    CLIN_SIG=paste(CLIN_SIG, collapse=";"),
                    gnomADe=paste(gnomADe, collapse=";"),
                    VAR_SYNONYMS=paste(VAR_SYNONYMS, collapse=";"),
                    AC=sum(AC),
                    AN=sum(AN),
                    RC=sum(AN) - sum(AC),
                    gnomADe_AC_eas=sum(gnomADe_AC_eas),
                    gnomADe_AN_eas=sum(gnomADe_AN_eas),
                    gnomADe_RC_eas=sum(gnomADe_AN_eas) - sum(gnomADe_AC_eas),
                    gnomADe_AC=sum(gnomADe_AC),
                    gnomADe_AN=sum(gnomADe_AN),
                    gnomADe_RC=sum(gnomADe_AN) - sum(gnomADe_AC), .groups="drop")
    df.gba <- as.data.frame(df.gba)
    # m <- match(df.gba$Gene, df$Gene)
    # df.gba$SYMBOL <- df$SYMBOL[m]

    samples <- colnames(df)[60:75]
    df.gba.2 <- df[, c("Gene", "SYMBOL", "Pos", samples)] %>%
        pivot_longer(cols=all_of(samples), names_to="sampleID", values_to="GT") %>%
        filter(!(GT %in% c("0/0", "0|0", "./.")))
    df.gba.var <- df.gba.2 %>%
        group_by(Gene, SYMBOL, Pos) %>%
        summarise(sampleID=paste(sampleID, "-", GT, collapse=","), .groups="drop")
    df.gba.var <- df.gba.var %>%
        group_by(Gene, SYMBOL) %>%
        summarise(variant_summary=paste0(Pos, " (", sampleID, ")", collapse="; "), .groups="drop")

    df.gba.sam <- df.gba.2 %>%
        group_by(Gene, SYMBOL, sampleID) %>%
        summarise(Pos=paste(Pos, collapse=","),
                    n=n(),
                    .groups="drop",
                    )
    df.gba.sam.1 <- df.gba.sam %>%
        group_by(Gene, SYMBOL) %>%
        summarise(sample_summary=paste0(sampleID, " (", Pos, ")", collapse="; "), .groups="drop")
    df.gba.sam.2 <- df.gba.sam %>%
        group_by(Gene, SYMBOL, n) %>%
        summarise(variant_stat=n(), .groups="drop")
    df.gba.sam.2 <- df.gba.sam.2 %>%
        group_by(Gene, SYMBOL) %>%
        summarize(variant_stat=paste0(variant_stat, " (", n, ")", collapse="; "), .groups="drop")

    m <- match(df.gba.var$Gene, df.gba$Gene)
    df.gba$variant_summary <- df.gba.var$variant_summary[m]

    m <- match(df.gba.sam.1$Gene, df.gba$Gene)
    df.gba$sample_summary <- df.gba.sam.1$sample_summary[m]

    m <- match(df.gba.sam.2$Gene, df.gba$Gene)
    df.gba$variant_stat <- df.gba.sam.2$variant_stat[m]  

    ### EAS gnomeADe
    t.gnomADe.eas <- apply(df.gba[, c("AC", "RC", "gnomADe_AC_eas", "gnomADe_RC_eas")], 1, function(x) {
        t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
        return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
    })
    t.gnomADe.eas <- as.data.frame(t(t.gnomADe.eas))
    colnames(t.gnomADe.eas)[2] <- "OR"
    adj.pVal.gnomADe.eas <- p.adjust(t.gnomADe.eas$pVal, method = "bonferroni")
    EF.gnomADe.eas <- (t.gnomADe.eas$OR-1)/t.gnomADe.eas$OR

    ### ALL gnomeADe
    t.gnomADe <- apply(df.gba[, c("AC", "RC", "gnomADe_AC", "gnomADe_RC")], 1, function(x) {
        t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
        return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
    })
    t.gnomADe <- as.data.frame(t(t.gnomADe))
    colnames(t.gnomADe)[2] <- "OR"
    adj.pVal.gnomADe <- p.adjust(t.gnomADe$pVal, method = "bonferroni")
    EF.gnomADe <- (t.gnomADe$OR-1)/t.gnomADe$OR

    df.res <- cbind(df.gba, 
                    p_value_EA_gnomADe=t.gnomADe.eas$pVal, 
                    FDR_EA_gnomADe=adj.pVal.gnomADe.eas, 
                    OR_EA_gnomADe=t.gnomADe.eas$OR,
                    CI_lower_EA_gnomADe=t.gnomADe.eas$CI_lower,
                    CI_upper_EA_gnomADe=t.gnomADe.eas$CI_upper,
                    EF_EA_gnomADe=EF.gnomADe.eas,   
                    p_value_ALL_gnomADe=t.gnomADe$pVal, 
                    FDR_ALL_gnomADe=adj.pVal.gnomADe,
                    OR_ALL_gnomADe=t.gnomADe$OR,
                    CI_lower_ALL_gnomADe=t.gnomADe$CI_lower,
                    CI_upper_ALL_gnomADe=t.gnomADe$CI_upper,
                    EF_ALL_gnomADe=EF.gnomADe,
                    IS.SIG=rep("ns", nrow(df.gba)))

    df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe > 0.05] <- "EA"
    df.res$IS.SIG[df.res$OR_ALL_gnomADe > 1 & df.res$FDR_ALL_gnomADe < 0.05 & df.res$FDR_EA_gnomADe > 0.05] <- "ALL"
    df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$OR_ALL_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe < 0.05] <- "EA/ALL"

    df.out <- cbind(df.res[, c("SYMBOL", "Gene", "No.of.variants", "variant_summary", "sample_summary", "variant_stat", "AC", "AN")],
                    `AC - gnomADe (EA)`=df.res$gnomADe_AC_eas, `AN - gnomADe (EA)`=df.res$gnomADe_AN_eas,
                    `OR (CI) - gnomADe (EA)`=paste0(round(df.res$OR_EA_gnomADe,2), " (", round(df.res$CI_lower_EA_gnomADe,2), "-", round(df.res$CI_upper_EA_gnomADe,2), ")"),
                    `p-value - gnomADe (EA)`=sprintf("%.2e", df.res$p_value_EA_gnomADe), `FDR - gnomADe (EA)`=sprintf("%.2e", df.res$FDR_EA_gnomADe), `EF - gnomADe (EA)`=round(df.res$EF_EA_gnomADe,2),
                    `AC - gnomADe (ALL)`=df.res$gnomADe_AC, `AN - gnomADe (ALL)`=df.res$gnomADe_AN,
                    `OR (CI) - gnomADe (ALL)`=paste0(round(df.res$OR_ALL_gnomADe,2), " (", round(df.res$CI_lower_ALL_gnomADe,2), "-", round(df.res$CI_upper_ALL_gnomADe,2), ")"),
                    `p-value - gnomADe (ALL)`=sprintf("%.2e", df.res$p_value_ALL_gnomADe), `FDR - gnomADe (ALL)`=sprintf("%.2e", df.res$FDR_ALL_gnomADe), `EF - gnomADe (ALL)`=round(df.res$EF_ALL_gnomADe,2),
                    is.sig=df.res$IS.SIG,
                    df.res[, c("Pos", "Consequence", "EXON", "HGVSc", "HGVSp", "SIFT", "PolyPhen", "CADD_RAW", "CADD_PHRED", "CLIN_SIG", "gnomADe", "VAR_SYNONYMS")])
    
    df.out <- df.out[order(factor(df.out$is.sig, levels=c("EA/ALL", "EA", "ALL", "ns"))),]

    write.table(df.res, paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".GBA.tsv"), sep="\t", quote=F, row.names=F)
    write.csv(df.out, paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".GBA.csv"))

    # ### MD genes
    # g <- c("FAM136A", "DTNA", "PRKCB", "COCH", "DPT", "SEMA3D", "TECTA" , "GUSB", "SLC6A7", "STRC", "HMX2", "PIP4P1", "OTOG", "LSAMP", "MYO7A", "ADGRV1", "CDH23", "PCDH15", "USH1C", "SHROOM2")
}

# indels
c <- "frameshift"
df <- read.table(paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".tsv"), header = T, sep = "\t", na.strings=".")

df.1 <- df[grep(",", df$ALT),]
df.1 <- df.1 %>% 
    separate_rows(ALT, AC, AF, sep=",")
df.1 <- as.data.frame(df.1)
df.1 <- df.1[c(1,4,5),]
df <- df[-grep(",", df$ALT),]
df <- rbind(df, df.1)
df$AF <- as.numeric(df$AF)
df$AC <- as.numeric(df$AC)
df <- df[df$AF <= 0.05,]
df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)

## SVA
df.sva <- df
df.sva$RC <- df.sva$AN - df.sva$AC
df.sva$gnomADe_RC <- df.sva$gnomADe_AN - df.sva$gnomADe_AC
df.sva$gnomADe_RC_eas <- df.sva$gnomADe_AN_eas - df.sva$gnomADe_AC_eas

### EAS gnomeADe
t.gnomADe.eas <- apply(df.sva[, c("AC", "RC", "gnomADe_AC_eas", "gnomADe_RC_eas")], 1, function(x) {
    t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
    return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
})
t.gnomADe.eas <- as.data.frame(t(t.gnomADe.eas))
colnames(t.gnomADe.eas)[2] <- "OR"
adj.pVal.gnomADe.eas <- p.adjust(t.gnomADe.eas$pVal, method = "bonferroni")
EF.gnomADe.eas <- (t.gnomADe.eas$OR-1)/t.gnomADe.eas$OR

### ALL gnomeADe
t.gnomADe <- apply(df.sva[, c("AC", "RC", "gnomADe_AC", "gnomADe_RC")], 1, function(x) {
    t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
    return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
})
t.gnomADe <- as.data.frame(t(t.gnomADe))
colnames(t.gnomADe)[2] <- "OR"
adj.pVal.gnomADe <- p.adjust(t.gnomADe$pVal, method = "bonferroni")
EF.gnomADe <- (t.gnomADe$OR-1)/t.gnomADe$OR

df.res <- cbind(df.sva,
                p_value_EA_gnomADe=t.gnomADe.eas$pVal, 
                FDR_EA_gnomADe=adj.pVal.gnomADe.eas, 
                OR_EA_gnomADe=t.gnomADe.eas$OR,
                CI_lower_EA_gnomADe=t.gnomADe.eas$CI_lower,
                CI_upper_EA_gnomADe=t.gnomADe.eas$CI_upper,
                EF_EA_gnomADe=EF.gnomADe.eas,
                p_value_ALL_gnomADe=t.gnomADe$pVal, 
                FDR_ALL_gnomADe=adj.pVal.gnomADe,
                OR_ALL_gnomADe=t.gnomADe$OR,
                CI_lower_ALL_gnomADe=t.gnomADe$CI_lower,
                CI_upper_ALL_gnomADe=t.gnomADe$CI_upper,
                EF_ALL_gnomADe=EF.gnomADe,
                IS.SIG=rep("ns", nrow(df.sva)))

df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe > 0.05] <- "EA"
df.res$IS.SIG[df.res$OR_ALL_gnomADe > 1 & df.res$FDR_ALL_gnomADe < 0.05 & df.res$FDR_EA_gnomADe > 0.05] <- "ALL"
df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$OR_ALL_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe < 0.05] <- "EA/ALL"

df.out <- cbind(df.res[, c("Pos", "SYMBOL", "Gene", "AC", "AN", "AF")],
                `AC - gnomADe (EA)`=df.res$gnomADe_AC_eas, `AN - gnomADe (EA)`=df.res$gnomADe_AN_eas, `AF - gnomADe (EA)`=sprintf("%.2e", df.res[, c("gnomADe_AF_eas")]),
                `OR (CI) - gnomADe (EA)`=paste0(round(df.res$OR_EA_gnomADe,2), " (", round(df.res$CI_lower_EA_gnomADe,2), "-", round(df.res$CI_upper_EA_gnomADe,2), ")"),
                `p-value - gnomADe (EA)`=sprintf("%.2e", df.res$p_value_EA_gnomADe), `FDR - gnomADe (EA)`=sprintf("%.2e", df.res$FDR_EA_gnomADe), `EF - gnomADe (EA)`=round(df.res$EF_EA_gnomADe,2),
                `AC - gnomADe (ALL)`=df.res$gnomADe_AC, `AN - gnomADe (ALL)`=df.res$gnomADe_AN, `AF - gnomADe (ALL)`=sprintf("%.2e", df.res[, c("gnomADe_AF")]),
                `OR (CI) - gnomADe (ALL)`=paste0(round(df.res$OR_ALL_gnomADe,2), " (", round(df.res$CI_lower_ALL_gnomADe,2), "-", round(df.res$CI_upper_ALL_gnomADe,2), ")"),
                `p-value - gnomADe (ALL)`=sprintf("%.2e", df.res$p_value_ALL_gnomADe), `FDR - gnomADe (ALL)`=sprintf("%.2e", df.res$FDR_ALL_gnomADe), `EF - gnomADe (ALL)`=round(df.res$EF_ALL_gnomADe,2),
                is.sig=df.res$IS.SIG,
                df.res[, c("Consequence", "EXON", "HGVSc", "HGVSp", "SIFT", "PolyPhen", "CADD_RAW", "CADD_PHRED", "CLIN_SIG", "gnomADe", "VAR_SYNONYMS")],
                df.res[, c(60:75)])

df.out <- df.out[order(factor(df.out$is.sig, levels=c("EA/ALL", "EA", "ALL", "ns"))),]

write.table(df.res, paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".SVA.tsv"), sep="\t", quote=F, row.names=F)
write.csv(df.out, paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".SVA.csv"))

## GBA
df.gba <- df
df.gba <- df.gba %>%
    group_by(Gene, SYMBOL) %>%
    summarise(No.of.variants=n(),
                Pos=paste(Pos, collapse=";"),
                Consequence=paste(Consequence, collapse=";"),
                EXON=paste(EXON, collapse=";"),
                HGVSc=paste(HGVSc, collapse=";"),
                HGVSp=paste(HGVSp, collapse=";"),
                SIFT=paste(SIFT, collapse=";"),
                PolyPhen=paste(PolyPhen, collapse=";"),
                CADD_RAW=paste(CADD_RAW, collapse=";"),
                CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                CLIN_SIG=paste(CLIN_SIG, collapse=";"),
                gnomADe=paste(gnomADe, collapse=";"),
                VAR_SYNONYMS=paste(VAR_SYNONYMS, collapse=";"),
                AC=sum(AC),
                AN=sum(AN),
                RC=sum(AN) - sum(AC),
                gnomADe_AC_eas=sum(gnomADe_AC_eas),
                gnomADe_AN_eas=sum(gnomADe_AN_eas),
                gnomADe_RC_eas=sum(gnomADe_AN_eas) - sum(gnomADe_AC_eas),
                gnomADe_AC=sum(gnomADe_AC),
                gnomADe_AN=sum(gnomADe_AN),
                gnomADe_RC=sum(gnomADe_AN) - sum(gnomADe_AC), .groups="drop")
df.gba <- as.data.frame(df.gba)
# m <- match(df.gba$Gene, df$Gene)
# df.gba$SYMBOL <- df$SYMBOL[m]

samples <- colnames(df)[60:75]
df.gba.2 <- df[, c("Gene", "SYMBOL", "Pos", samples)] %>%
    pivot_longer(cols=all_of(samples), names_to="sampleID", values_to="GT") %>%
    filter(!(GT %in% c("0/0", "0|0", "./.")))
df.gba.var <- df.gba.2 %>%
    group_by(Gene, SYMBOL, Pos) %>%
    summarise(sampleID=paste(sampleID, "-", GT, collapse=","), .groups="drop")
df.gba.var <- df.gba.var %>%
    group_by(Gene, SYMBOL) %>%
    summarise(variant_summary=paste0(Pos, " (", sampleID, ")", collapse="; "), .groups="drop")

df.gba.sam <- df.gba.2 %>%
    group_by(Gene, SYMBOL, sampleID) %>%
    summarise(Pos=paste(Pos, collapse=","),
                n=n(),
                .groups="drop",
                )
df.gba.sam.1 <- df.gba.sam %>%
    group_by(Gene, SYMBOL) %>%
    summarise(sample_summary=paste0(sampleID, " (", Pos, ")", collapse="; "), .groups="drop")
df.gba.sam.2 <- df.gba.sam %>%
    group_by(Gene, SYMBOL, n) %>%
    summarise(variant_stat=n(), .groups="drop")
df.gba.sam.2 <- df.gba.sam.2 %>%
    group_by(Gene, SYMBOL) %>%
    summarize(variant_stat=paste0(variant_stat, " (", n, ")", collapse="; "), .groups="drop")

m <- match(df.gba.var$Gene, df.gba$Gene)
df.gba$variant_summary <- df.gba.var$variant_summary[m]

m <- match(df.gba.sam.1$Gene, df.gba$Gene)
df.gba$sample_summary <- df.gba.sam.1$sample_summary[m]

m <- match(df.gba.sam.2$Gene, df.gba$Gene)
df.gba$variant_stat <- df.gba.sam.2$variant_stat[m]  

### EAS gnomeADe
t.gnomADe.eas <- apply(df.gba[, c("AC", "RC", "gnomADe_AC_eas", "gnomADe_RC_eas")], 1, function(x) {
    t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
    return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
})
t.gnomADe.eas <- as.data.frame(t(t.gnomADe.eas))
colnames(t.gnomADe.eas)[2] <- "OR"
adj.pVal.gnomADe.eas <- p.adjust(t.gnomADe.eas$pVal, method = "bonferroni")
EF.gnomADe.eas <- (t.gnomADe.eas$OR-1)/t.gnomADe.eas$OR

### ALL gnomeADe
t.gnomADe <- apply(df.gba[, c("AC", "RC", "gnomADe_AC", "gnomADe_RC")], 1, function(x) {
    t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
    return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
})
t.gnomADe <- as.data.frame(t(t.gnomADe))
colnames(t.gnomADe)[2] <- "OR"
adj.pVal.gnomADe <- p.adjust(t.gnomADe$pVal, method = "bonferroni")
EF.gnomADe <- (t.gnomADe$OR-1)/t.gnomADe$OR

df.res <- cbind(df.gba, 
                p_value_EA_gnomADe=t.gnomADe.eas$pVal, 
                FDR_EA_gnomADe=adj.pVal.gnomADe.eas, 
                OR_EA_gnomADe=t.gnomADe.eas$OR,
                CI_lower_EA_gnomADe=t.gnomADe.eas$CI_lower,
                CI_upper_EA_gnomADe=t.gnomADe.eas$CI_upper,
                EF_EA_gnomADe=EF.gnomADe.eas,   
                p_value_ALL_gnomADe=t.gnomADe$pVal, 
                FDR_ALL_gnomADe=adj.pVal.gnomADe,
                OR_ALL_gnomADe=t.gnomADe$OR,
                CI_lower_ALL_gnomADe=t.gnomADe$CI_lower,
                CI_upper_ALL_gnomADe=t.gnomADe$CI_upper,
                EF_ALL_gnomADe=EF.gnomADe,
                IS.SIG=rep("ns", nrow(df.gba)))

df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe > 0.05] <- "EA"
df.res$IS.SIG[df.res$OR_ALL_gnomADe > 1 & df.res$FDR_ALL_gnomADe < 0.05 & df.res$FDR_EA_gnomADe > 0.05] <- "ALL"
df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$OR_ALL_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe < 0.05] <- "EA/ALL"

df.out <- cbind(df.res[, c("SYMBOL", "Gene", "No.of.variants", "variant_summary", "sample_summary", "variant_stat", "AC", "AN")],
                `AC - gnomADe (EA)`=df.res$gnomADe_AC_eas, `AN - gnomADe (EA)`=df.res$gnomADe_AN_eas,
                `OR (CI) - gnomADe (EA)`=paste0(round(df.res$OR_EA_gnomADe,2), " (", round(df.res$CI_lower_EA_gnomADe,2), "-", round(df.res$CI_upper_EA_gnomADe,2), ")"),
                `p-value - gnomADe (EA)`=sprintf("%.2e", df.res$p_value_EA_gnomADe), `FDR - gnomADe (EA)`=sprintf("%.2e", df.res$FDR_EA_gnomADe), `EF - gnomADe (EA)`=round(df.res$EF_EA_gnomADe,2),
                `AC - gnomADe (ALL)`=df.res$gnomADe_AC, `AN - gnomADe (ALL)`=df.res$gnomADe_AN,
                `OR (CI) - gnomADe (ALL)`=paste0(round(df.res$OR_ALL_gnomADe,2), " (", round(df.res$CI_lower_ALL_gnomADe,2), "-", round(df.res$CI_upper_ALL_gnomADe,2), ")"),
                `p-value - gnomADe (ALL)`=sprintf("%.2e", df.res$p_value_ALL_gnomADe), `FDR - gnomADe (ALL)`=sprintf("%.2e", df.res$FDR_ALL_gnomADe), `EF - gnomADe (ALL)`=round(df.res$EF_ALL_gnomADe,2),
                is.sig=df.res$IS.SIG,
                df.res[, c("Pos", "Consequence", "EXON", "HGVSc", "HGVSp", "SIFT", "PolyPhen", "CADD_RAW", "CADD_PHRED", "CLIN_SIG", "gnomADe", "VAR_SYNONYMS")])

df.out <- df.out[order(factor(df.out$is.sig, levels=c("EA/ALL", "EA", "ALL", "ns"))),]

write.table(df.res, paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".GBA.tsv"), sep="\t", quote=F, row.names=F)
write.csv(df.out, paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".GBA.csv"))


## Dealing with 'match v2'
library(tidyr)
library(dplyr)

# SNV
setwd("/home/mdnl/Documents/Tam_14Mar25/MD/variant_calling/WES_Korean/vep")

cons_type <- c("missense", "LoF")

for (c in cons_type) {
    print(c)
    df <- read.table(paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.tsv"), header = T, sep = "\t", na.strings=".")
    df <- df %>% 
        separate_rows(ALT, AC, AF, sep=",") %>%
        filter(ALT==Allele)
    df <- as.data.frame(df)
    df$AF <- as.numeric(df$AF)
    df$AC <- as.numeric(df$AC)
    df <- df[df$AF <= 0.05,]
    df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)

    ## SVA
    df.sva <- df
    df.sva$RC <- df.sva$AN - df.sva$AC
    df.sva$gnomADe_RC <- df.sva$gnomADe_AN - df.sva$gnomADe_AC
    df.sva$gnomADe_RC_eas <- df.sva$gnomADe_AN_eas - df.sva$gnomADe_AC_eas

    ### EAS gnomeADe
    t.gnomADe.eas <- apply(df.sva[, c("AC", "RC", "gnomADe_AC_eas", "gnomADe_RC_eas")], 1, function(x) {
        t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
        return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
    })
    t.gnomADe.eas <- as.data.frame(t(t.gnomADe.eas))
    colnames(t.gnomADe.eas)[2] <- "OR"
    adj.pVal.gnomADe.eas <- p.adjust(t.gnomADe.eas$pVal, method = "bonferroni")
    EF.gnomADe.eas <- (t.gnomADe.eas$OR-1)/t.gnomADe.eas$OR

    ### ALL gnomeADe
    t.gnomADe <- apply(df.sva[, c("AC", "RC", "gnomADe_AC", "gnomADe_RC")], 1, function(x) {
        t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
        return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
    })
    t.gnomADe <- as.data.frame(t(t.gnomADe))
    colnames(t.gnomADe)[2] <- "OR"
    adj.pVal.gnomADe <- p.adjust(t.gnomADe$pVal, method = "bonferroni")
    EF.gnomADe <- (t.gnomADe$OR-1)/t.gnomADe$OR

    df.res <- cbind(df.sva,
                    p_value_EA_gnomADe=t.gnomADe.eas$pVal, 
                    FDR_EA_gnomADe=adj.pVal.gnomADe.eas, 
                    OR_EA_gnomADe=t.gnomADe.eas$OR,
                    CI_lower_EA_gnomADe=t.gnomADe.eas$CI_lower,
                    CI_upper_EA_gnomADe=t.gnomADe.eas$CI_upper,
                    EF_EA_gnomADe=EF.gnomADe.eas,
                    p_value_ALL_gnomADe=t.gnomADe$pVal, 
                    FDR_ALL_gnomADe=adj.pVal.gnomADe,
                    OR_ALL_gnomADe=t.gnomADe$OR,
                    CI_lower_ALL_gnomADe=t.gnomADe$CI_lower,
                    CI_upper_ALL_gnomADe=t.gnomADe$CI_upper,
                    EF_ALL_gnomADe=EF.gnomADe,
                    IS.SIG=rep("ns", nrow(df.sva)))

    df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe > 0.05] <- "EA"
    df.res$IS.SIG[df.res$OR_ALL_gnomADe > 1 & df.res$FDR_ALL_gnomADe < 0.05 & df.res$FDR_EA_gnomADe > 0.05] <- "ALL"
    df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$OR_ALL_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe < 0.05] <- "EA/ALL"

    df.out <- cbind(df.res[, c("Pos", "SYMBOL", "Gene", "AC", "AN", "AF")],
                    `AC - gnomADe (EA)`=df.res$gnomADe_AC_eas, `AN - gnomADe (EA)`=df.res$gnomADe_AN_eas, `AF - gnomADe (EA)`=sprintf("%.2e", df.res[, c("gnomADe_AF_eas")]),
                    `OR (CI) - gnomADe (EA)`=paste0(round(df.res$OR_EA_gnomADe,2), " (", round(df.res$CI_lower_EA_gnomADe,2), "-", round(df.res$CI_upper_EA_gnomADe,2), ")"),
                    `p-value - gnomADe (EA)`=sprintf("%.2e", df.res$p_value_EA_gnomADe), `FDR - gnomADe (EA)`=sprintf("%.2e", df.res$FDR_EA_gnomADe), `EF - gnomADe (EA)`=round(df.res$EF_EA_gnomADe,2),
                    `AC - gnomADe (ALL)`=df.res$gnomADe_AC, `AN - gnomADe (ALL)`=df.res$gnomADe_AN, `AF - gnomADe (ALL)`=sprintf("%.2e", df.res[, c("gnomADe_AF")]),
                    `OR (CI) - gnomADe (ALL)`=paste0(round(df.res$OR_ALL_gnomADe,2), " (", round(df.res$CI_lower_ALL_gnomADe,2), "-", round(df.res$CI_upper_ALL_gnomADe,2), ")"),
                    `p-value - gnomADe (ALL)`=sprintf("%.2e", df.res$p_value_ALL_gnomADe), `FDR - gnomADe (ALL)`=sprintf("%.2e", df.res$FDR_ALL_gnomADe), `EF - gnomADe (ALL)`=round(df.res$EF_ALL_gnomADe,2),
                    is.sig=df.res$IS.SIG,
                    df.res[, c("Consequence", "EXON", "HGVSc", "HGVSp", "SIFT", "PolyPhen", "CADD_RAW", "CADD_PHRED", "CLIN_SIG", "gnomADe", "VAR_SYNONYMS")],
                    df.res[, c(60:75)])

    df.out <- df.out[order(factor(df.out$is.sig, levels=c("EA/ALL", "EA", "ALL", "ns"))),]
    
    write.table(df.res, paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.SVA.tsv"), sep="\t", quote=F, row.names=F)
    write.csv(df.out, paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.SVA.csv"))

    ## GBA
    df.gba <- df
    df.gba <- df.gba %>%
        group_by(Gene, SYMBOL) %>%
        summarise(No.of.variants=n(),
                    Pos=paste(Pos, collapse=";"),
                    Consequence=paste(Consequence, collapse=";"),
                    EXON=paste(EXON, collapse=";"),
                    HGVSc=paste(HGVSc, collapse=";"),
                    HGVSp=paste(HGVSp, collapse=";"),
                    SIFT=paste(SIFT, collapse=";"),
                    PolyPhen=paste(PolyPhen, collapse=";"),
                    CADD_RAW=paste(CADD_RAW, collapse=";"),
                    CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                    CLIN_SIG=paste(CLIN_SIG, collapse=";"),
                    gnomADe=paste(gnomADe, collapse=";"),
                    VAR_SYNONYMS=paste(VAR_SYNONYMS, collapse=";"),
                    AC=sum(AC),
                    AN=sum(AN),
                    RC=sum(AN) - sum(AC),
                    gnomADe_AC_eas=sum(gnomADe_AC_eas),
                    gnomADe_AN_eas=sum(gnomADe_AN_eas),
                    gnomADe_RC_eas=sum(gnomADe_AN_eas) - sum(gnomADe_AC_eas),
                    gnomADe_AC=sum(gnomADe_AC),
                    gnomADe_AN=sum(gnomADe_AN),
                    gnomADe_RC=sum(gnomADe_AN) - sum(gnomADe_AC), .groups="drop")
    df.gba <- as.data.frame(df.gba)
    # m <- match(df.gba$Gene, df$Gene)
    # df.gba$SYMBOL <- df$SYMBOL[m]

    samples <- colnames(df)[60:75]
    df.gba.2 <- df[, c("Gene", "SYMBOL", "Pos", samples)] %>%
        pivot_longer(cols=all_of(samples), names_to="sampleID", values_to="GT") %>%
        filter(!(GT %in% c("0/0", "0|0", "./.")))
    df.gba.var <- df.gba.2 %>%
        group_by(Gene, SYMBOL, Pos) %>%
        summarise(sampleID=paste(sampleID, "-", GT, collapse=","), .groups="drop")
    df.gba.var <- df.gba.var %>%
        group_by(Gene, SYMBOL) %>%
        summarise(variant_summary=paste0(Pos, " (", sampleID, ")", collapse="; "), .groups="drop")

    df.gba.sam <- df.gba.2 %>%
        group_by(Gene, SYMBOL, sampleID) %>%
        summarise(Pos=paste(Pos, collapse=","),
                    n=n(),
                    .groups="drop",
                    )
    df.gba.sam.1 <- df.gba.sam %>%
        group_by(Gene, SYMBOL) %>%
        summarise(sample_summary=paste0(sampleID, " (", Pos, ")", collapse="; "), .groups="drop")
    df.gba.sam.2 <- df.gba.sam %>%
        group_by(Gene, SYMBOL, n) %>%
        summarise(variant_stat=n(), .groups="drop")
    df.gba.sam.2 <- df.gba.sam.2 %>%
        group_by(Gene, SYMBOL) %>%
        summarize(variant_stat=paste0(variant_stat, " (", n, ")", collapse="; "), .groups="drop")

    m <- match(df.gba.var$Gene, df.gba$Gene)
    df.gba$variant_summary <- df.gba.var$variant_summary[m]

    m <- match(df.gba.sam.1$Gene, df.gba$Gene)
    df.gba$sample_summary <- df.gba.sam.1$sample_summary[m]

    m <- match(df.gba.sam.2$Gene, df.gba$Gene)
    df.gba$variant_stat <- df.gba.sam.2$variant_stat[m]  

    ### EAS gnomeADe
    t.gnomADe.eas <- apply(df.gba[, c("AC", "RC", "gnomADe_AC_eas", "gnomADe_RC_eas")], 1, function(x) {
        t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
        return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
    })
    t.gnomADe.eas <- as.data.frame(t(t.gnomADe.eas))
    colnames(t.gnomADe.eas)[2] <- "OR"
    adj.pVal.gnomADe.eas <- p.adjust(t.gnomADe.eas$pVal, method = "bonferroni")
    EF.gnomADe.eas <- (t.gnomADe.eas$OR-1)/t.gnomADe.eas$OR

    ### ALL gnomeADe
    t.gnomADe <- apply(df.gba[, c("AC", "RC", "gnomADe_AC", "gnomADe_RC")], 1, function(x) {
        t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
        return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
    })
    t.gnomADe <- as.data.frame(t(t.gnomADe))
    colnames(t.gnomADe)[2] <- "OR"
    adj.pVal.gnomADe <- p.adjust(t.gnomADe$pVal, method = "bonferroni")
    EF.gnomADe <- (t.gnomADe$OR-1)/t.gnomADe$OR

    df.res <- cbind(df.gba, 
                    p_value_EA_gnomADe=t.gnomADe.eas$pVal, 
                    FDR_EA_gnomADe=adj.pVal.gnomADe.eas, 
                    OR_EA_gnomADe=t.gnomADe.eas$OR,
                    CI_lower_EA_gnomADe=t.gnomADe.eas$CI_lower,
                    CI_upper_EA_gnomADe=t.gnomADe.eas$CI_upper,
                    EF_EA_gnomADe=EF.gnomADe.eas,   
                    p_value_ALL_gnomADe=t.gnomADe$pVal, 
                    FDR_ALL_gnomADe=adj.pVal.gnomADe,
                    OR_ALL_gnomADe=t.gnomADe$OR,
                    CI_lower_ALL_gnomADe=t.gnomADe$CI_lower,
                    CI_upper_ALL_gnomADe=t.gnomADe$CI_upper,
                    EF_ALL_gnomADe=EF.gnomADe,
                    IS.SIG=rep("ns", nrow(df.gba)))

    df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe > 0.05] <- "EA"
    df.res$IS.SIG[df.res$OR_ALL_gnomADe > 1 & df.res$FDR_ALL_gnomADe < 0.05 & df.res$FDR_EA_gnomADe > 0.05] <- "ALL"
    df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$OR_ALL_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe < 0.05] <- "EA/ALL"

    df.out <- cbind(df.res[, c("SYMBOL", "Gene", "No.of.variants", "variant_summary", "sample_summary", "variant_stat", "AC", "AN")],
                    `AC - gnomADe (EA)`=df.res$gnomADe_AC_eas, `AN - gnomADe (EA)`=df.res$gnomADe_AN_eas,
                    `OR (CI) - gnomADe (EA)`=paste0(round(df.res$OR_EA_gnomADe,2), " (", round(df.res$CI_lower_EA_gnomADe,2), "-", round(df.res$CI_upper_EA_gnomADe,2), ")"),
                    `p-value - gnomADe (EA)`=sprintf("%.2e", df.res$p_value_EA_gnomADe), `FDR - gnomADe (EA)`=sprintf("%.2e", df.res$FDR_EA_gnomADe), `EF - gnomADe (EA)`=round(df.res$EF_EA_gnomADe,2),
                    `AC - gnomADe (ALL)`=df.res$gnomADe_AC, `AN - gnomADe (ALL)`=df.res$gnomADe_AN,
                    `OR (CI) - gnomADe (ALL)`=paste0(round(df.res$OR_ALL_gnomADe,2), " (", round(df.res$CI_lower_ALL_gnomADe,2), "-", round(df.res$CI_upper_ALL_gnomADe,2), ")"),
                    `p-value - gnomADe (ALL)`=sprintf("%.2e", df.res$p_value_ALL_gnomADe), `FDR - gnomADe (ALL)`=sprintf("%.2e", df.res$FDR_ALL_gnomADe), `EF - gnomADe (ALL)`=round(df.res$EF_ALL_gnomADe,2),
                    is.sig=df.res$IS.SIG,
                    df.res[, c("Pos", "Consequence", "EXON", "HGVSc", "HGVSp", "SIFT", "PolyPhen", "CADD_RAW", "CADD_PHRED", "CLIN_SIG", "gnomADe", "VAR_SYNONYMS")])
    
    df.out <- df.out[order(factor(df.out$is.sig, levels=c("EA/ALL", "EA", "ALL", "ns"))),]

    write.table(df.res, paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.GBA.tsv"), sep="\t", quote=F, row.names=F)
    write.csv(df.out, paste0("cohort.snp.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.GBA.csv"))

    # ### MD genes
    # g <- c("FAM136A", "DTNA", "PRKCB", "COCH", "DPT", "SEMA3D", "TECTA" , "GUSB", "SLC6A7", "STRC", "HMX2", "PIP4P1", "OTOG", "LSAMP", "MYO7A", "ADGRV1", "CDH23", "PCDH15", "USH1C", "SHROOM2")
}

# indels
c <- "frameshift"
df <- read.table(paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.tsv"), header = T, sep = "\t", na.strings=".")

df.1 <- df[grep(",", df$ALT),]
df.1 <- df.1 %>% 
    separate_rows(ALT, AC, AF, sep=",")
df.1 <- as.data.frame(df.1)
df.1 <- df.1[c(4),]
df <- df[-grep(",", df$ALT),]
df <- rbind(df, df.1)
df$AF <- as.numeric(df$AF)
df$AC <- as.numeric(df$AC)
df <- df[df$AF <= 0.05,]
df$Pos <- paste0(df$CHROM, ":", df$POS, ":", df$REF, ">", df$ALT)

## SVA
df.sva <- df
df.sva$RC <- df.sva$AN - df.sva$AC
df.sva$gnomADe_RC <- df.sva$gnomADe_AN - df.sva$gnomADe_AC
df.sva$gnomADe_RC_eas <- df.sva$gnomADe_AN_eas - df.sva$gnomADe_AC_eas

### EAS gnomeADe
t.gnomADe.eas <- apply(df.sva[, c("AC", "RC", "gnomADe_AC_eas", "gnomADe_RC_eas")], 1, function(x) {
    t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
    return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
})
t.gnomADe.eas <- as.data.frame(t(t.gnomADe.eas))
colnames(t.gnomADe.eas)[2] <- "OR"
adj.pVal.gnomADe.eas <- p.adjust(t.gnomADe.eas$pVal, method = "bonferroni")
EF.gnomADe.eas <- (t.gnomADe.eas$OR-1)/t.gnomADe.eas$OR

### ALL gnomeADe
t.gnomADe <- apply(df.sva[, c("AC", "RC", "gnomADe_AC", "gnomADe_RC")], 1, function(x) {
    t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
    return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
})
t.gnomADe <- as.data.frame(t(t.gnomADe))
colnames(t.gnomADe)[2] <- "OR"
adj.pVal.gnomADe <- p.adjust(t.gnomADe$pVal, method = "bonferroni")
EF.gnomADe <- (t.gnomADe$OR-1)/t.gnomADe$OR

df.res <- cbind(df.sva,
                p_value_EA_gnomADe=t.gnomADe.eas$pVal, 
                FDR_EA_gnomADe=adj.pVal.gnomADe.eas, 
                OR_EA_gnomADe=t.gnomADe.eas$OR,
                CI_lower_EA_gnomADe=t.gnomADe.eas$CI_lower,
                CI_upper_EA_gnomADe=t.gnomADe.eas$CI_upper,
                EF_EA_gnomADe=EF.gnomADe.eas,
                p_value_ALL_gnomADe=t.gnomADe$pVal, 
                FDR_ALL_gnomADe=adj.pVal.gnomADe,
                OR_ALL_gnomADe=t.gnomADe$OR,
                CI_lower_ALL_gnomADe=t.gnomADe$CI_lower,
                CI_upper_ALL_gnomADe=t.gnomADe$CI_upper,
                EF_ALL_gnomADe=EF.gnomADe,
                IS.SIG=rep("ns", nrow(df.sva)))

df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe > 0.05] <- "EA"
df.res$IS.SIG[df.res$OR_ALL_gnomADe > 1 & df.res$FDR_ALL_gnomADe < 0.05 & df.res$FDR_EA_gnomADe > 0.05] <- "ALL"
df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$OR_ALL_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe < 0.05] <- "EA/ALL"

df.out <- cbind(df.res[, c("Pos", "SYMBOL", "Gene", "AC", "AN", "AF")],
                `AC - gnomADe (EA)`=df.res$gnomADe_AC_eas, `AN - gnomADe (EA)`=df.res$gnomADe_AN_eas, `AF - gnomADe (EA)`=sprintf("%.2e", df.res[, c("gnomADe_AF_eas")]),
                `OR (CI) - gnomADe (EA)`=paste0(round(df.res$OR_EA_gnomADe,2), " (", round(df.res$CI_lower_EA_gnomADe,2), "-", round(df.res$CI_upper_EA_gnomADe,2), ")"),
                `p-value - gnomADe (EA)`=sprintf("%.2e", df.res$p_value_EA_gnomADe), `FDR - gnomADe (EA)`=sprintf("%.2e", df.res$FDR_EA_gnomADe), `EF - gnomADe (EA)`=round(df.res$EF_EA_gnomADe,2),
                `AC - gnomADe (ALL)`=df.res$gnomADe_AC, `AN - gnomADe (ALL)`=df.res$gnomADe_AN, `AF - gnomADe (ALL)`=sprintf("%.2e", df.res[, c("gnomADe_AF")]),
                `OR (CI) - gnomADe (ALL)`=paste0(round(df.res$OR_ALL_gnomADe,2), " (", round(df.res$CI_lower_ALL_gnomADe,2), "-", round(df.res$CI_upper_ALL_gnomADe,2), ")"),
                `p-value - gnomADe (ALL)`=sprintf("%.2e", df.res$p_value_ALL_gnomADe), `FDR - gnomADe (ALL)`=sprintf("%.2e", df.res$FDR_ALL_gnomADe), `EF - gnomADe (ALL)`=round(df.res$EF_ALL_gnomADe,2),
                is.sig=df.res$IS.SIG,
                df.res[, c("Consequence", "EXON", "HGVSc", "HGVSp", "SIFT", "PolyPhen", "CADD_RAW", "CADD_PHRED", "CLIN_SIG", "gnomADe", "VAR_SYNONYMS")],
                df.res[, c(60:75)])

df.out <- df.out[order(factor(df.out$is.sig, levels=c("EA/ALL", "EA", "ALL", "ns"))),]

write.table(df.res, paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.SVA.tsv"), sep="\t", quote=F, row.names=F)
write.csv(df.out, paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.SVA.csv"))

## GBA
df.gba <- df
df.gba <- df.gba %>%
    group_by(Gene, SYMBOL) %>%
    summarise(No.of.variants=n(),
                Pos=paste(Pos, collapse=";"),
                Consequence=paste(Consequence, collapse=";"),
                EXON=paste(EXON, collapse=";"),
                HGVSc=paste(HGVSc, collapse=";"),
                HGVSp=paste(HGVSp, collapse=";"),
                SIFT=paste(SIFT, collapse=";"),
                PolyPhen=paste(PolyPhen, collapse=";"),
                CADD_RAW=paste(CADD_RAW, collapse=";"),
                CADD_PHRED=paste(CADD_PHRED, collapse=";"),
                CLIN_SIG=paste(CLIN_SIG, collapse=";"),
                gnomADe=paste(gnomADe, collapse=";"),
                VAR_SYNONYMS=paste(VAR_SYNONYMS, collapse=";"),
                AC=sum(AC),
                AN=sum(AN),
                RC=sum(AN) - sum(AC),
                gnomADe_AC_eas=sum(gnomADe_AC_eas),
                gnomADe_AN_eas=sum(gnomADe_AN_eas),
                gnomADe_RC_eas=sum(gnomADe_AN_eas) - sum(gnomADe_AC_eas),
                gnomADe_AC=sum(gnomADe_AC),
                gnomADe_AN=sum(gnomADe_AN),
                gnomADe_RC=sum(gnomADe_AN) - sum(gnomADe_AC), .groups="drop")
df.gba <- as.data.frame(df.gba)
# m <- match(df.gba$Gene, df$Gene)
# df.gba$SYMBOL <- df$SYMBOL[m]

samples <- colnames(df)[60:75]
df.gba.2 <- df[, c("Gene", "SYMBOL", "Pos", samples)] %>%
    pivot_longer(cols=all_of(samples), names_to="sampleID", values_to="GT") %>%
    filter(!(GT %in% c("0/0", "0|0", "./.")))
df.gba.var <- df.gba.2 %>%
    group_by(Gene, SYMBOL, Pos) %>%
    summarise(sampleID=paste(sampleID, "-", GT, collapse=","), .groups="drop")
df.gba.var <- df.gba.var %>%
    group_by(Gene, SYMBOL) %>%
    summarise(variant_summary=paste0(Pos, " (", sampleID, ")", collapse="; "), .groups="drop")

df.gba.sam <- df.gba.2 %>%
    group_by(Gene, SYMBOL, sampleID) %>%
    summarise(Pos=paste(Pos, collapse=","),
                n=n(),
                .groups="drop",
                )
df.gba.sam.1 <- df.gba.sam %>%
    group_by(Gene, SYMBOL) %>%
    summarise(sample_summary=paste0(sampleID, " (", Pos, ")", collapse="; "), .groups="drop")
df.gba.sam.2 <- df.gba.sam %>%
    group_by(Gene, SYMBOL, n) %>%
    summarise(variant_stat=n(), .groups="drop")
df.gba.sam.2 <- df.gba.sam.2 %>%
    group_by(Gene, SYMBOL) %>%
    summarize(variant_stat=paste0(variant_stat, " (", n, ")", collapse="; "), .groups="drop")

m <- match(df.gba.var$Gene, df.gba$Gene)
df.gba$variant_summary <- df.gba.var$variant_summary[m]

m <- match(df.gba.sam.1$Gene, df.gba$Gene)
df.gba$sample_summary <- df.gba.sam.1$sample_summary[m]

m <- match(df.gba.sam.2$Gene, df.gba$Gene)
df.gba$variant_stat <- df.gba.sam.2$variant_stat[m]  

### EAS gnomeADe
t.gnomADe.eas <- apply(df.gba[, c("AC", "RC", "gnomADe_AC_eas", "gnomADe_RC_eas")], 1, function(x) {
    t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
    return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
})
t.gnomADe.eas <- as.data.frame(t(t.gnomADe.eas))
colnames(t.gnomADe.eas)[2] <- "OR"
adj.pVal.gnomADe.eas <- p.adjust(t.gnomADe.eas$pVal, method = "bonferroni")
EF.gnomADe.eas <- (t.gnomADe.eas$OR-1)/t.gnomADe.eas$OR

### ALL gnomeADe
t.gnomADe <- apply(df.gba[, c("AC", "RC", "gnomADe_AC", "gnomADe_RC")], 1, function(x) {
    t_res <- fisher.test(matrix(x, nrow = 2, byrow=T))
    return(c(pVal=t_res$p.value, OR=t_res$estimate, CI_lower=t_res$conf.int[1], CI_upper=t_res$conf.int[2]))
})
t.gnomADe <- as.data.frame(t(t.gnomADe))
colnames(t.gnomADe)[2] <- "OR"
adj.pVal.gnomADe <- p.adjust(t.gnomADe$pVal, method = "bonferroni")
EF.gnomADe <- (t.gnomADe$OR-1)/t.gnomADe$OR

df.res <- cbind(df.gba, 
                p_value_EA_gnomADe=t.gnomADe.eas$pVal, 
                FDR_EA_gnomADe=adj.pVal.gnomADe.eas, 
                OR_EA_gnomADe=t.gnomADe.eas$OR,
                CI_lower_EA_gnomADe=t.gnomADe.eas$CI_lower,
                CI_upper_EA_gnomADe=t.gnomADe.eas$CI_upper,
                EF_EA_gnomADe=EF.gnomADe.eas,   
                p_value_ALL_gnomADe=t.gnomADe$pVal, 
                FDR_ALL_gnomADe=adj.pVal.gnomADe,
                OR_ALL_gnomADe=t.gnomADe$OR,
                CI_lower_ALL_gnomADe=t.gnomADe$CI_lower,
                CI_upper_ALL_gnomADe=t.gnomADe$CI_upper,
                EF_ALL_gnomADe=EF.gnomADe,
                IS.SIG=rep("ns", nrow(df.gba)))

df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe > 0.05] <- "EA"
df.res$IS.SIG[df.res$OR_ALL_gnomADe > 1 & df.res$FDR_ALL_gnomADe < 0.05 & df.res$FDR_EA_gnomADe > 0.05] <- "ALL"
df.res$IS.SIG[df.res$OR_EA_gnomADe > 1 & df.res$OR_ALL_gnomADe > 1 & df.res$FDR_EA_gnomADe < 0.05 & df.res$FDR_ALL_gnomADe < 0.05] <- "EA/ALL"

df.out <- cbind(df.res[, c("SYMBOL", "Gene", "No.of.variants", "variant_summary", "sample_summary", "variant_stat", "AC", "AN")],
                `AC - gnomADe (EA)`=df.res$gnomADe_AC_eas, `AN - gnomADe (EA)`=df.res$gnomADe_AN_eas,
                `OR (CI) - gnomADe (EA)`=paste0(round(df.res$OR_EA_gnomADe,2), " (", round(df.res$CI_lower_EA_gnomADe,2), "-", round(df.res$CI_upper_EA_gnomADe,2), ")"),
                `p-value - gnomADe (EA)`=sprintf("%.2e", df.res$p_value_EA_gnomADe), `FDR - gnomADe (EA)`=sprintf("%.2e", df.res$FDR_EA_gnomADe), `EF - gnomADe (EA)`=round(df.res$EF_EA_gnomADe,2),
                `AC - gnomADe (ALL)`=df.res$gnomADe_AC, `AN - gnomADe (ALL)`=df.res$gnomADe_AN,
                `OR (CI) - gnomADe (ALL)`=paste0(round(df.res$OR_ALL_gnomADe,2), " (", round(df.res$CI_lower_ALL_gnomADe,2), "-", round(df.res$CI_upper_ALL_gnomADe,2), ")"),
                `p-value - gnomADe (ALL)`=sprintf("%.2e", df.res$p_value_ALL_gnomADe), `FDR - gnomADe (ALL)`=sprintf("%.2e", df.res$FDR_ALL_gnomADe), `EF - gnomADe (ALL)`=round(df.res$EF_ALL_gnomADe,2),
                is.sig=df.res$IS.SIG,
                df.res[, c("Pos", "Consequence", "EXON", "HGVSc", "HGVSp", "SIFT", "PolyPhen", "CADD_RAW", "CADD_PHRED", "CLIN_SIG", "gnomADe", "VAR_SYNONYMS")])

df.out <- df.out[order(factor(df.out$is.sig, levels=c("EA/ALL", "EA", "ALL", "ns"))),]

write.table(df.res, paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.GBA.tsv"), sep="\t", quote=F, row.names=F)
write.csv(df.out, paste0("cohort.indel.90.0.hard_filter.MAF.0.05.rare.vep_filt.", c, ".v2.GBA.csv"))
