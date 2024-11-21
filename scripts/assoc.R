library(dplyr)
library(data.table)
library(gghighlight)
library(tidyr)
library(ggplot2)
library(ggsci)
library(AssotesteR)
library(vartools)
library(GenABEL)
library(metaSeq)
library(metap)
library(scran)
library(patchwork)
library(gridExtra)
library(stringr)

source('scripts/functions.R')  # Upload needed functions.


df_eur <- as_data_frame(fread('data/pre_gwas_datasets/eur.worst.tsv', na.strings = c('.')))  # Eur dataset after extracting from Hail's MatrixTable.
df_afr <- as_data_frame(fread('data/pre_gwas_datasets/afr.worst.tsv', na.strings = c('.')))
df_amr <- as_data_frame(fread('data/pre_gwas_datasets/amr.worst.tsv', na.strings = c('.')))

df_eur <- prep_tmp(df_eur)
df_afr <- prep_tmp(df_afr)
df_amr <- prep_tmp(df_amr)

n_case_eur <- df_eur$non_miss_case %>% max()  # Suppose max number of cases is the same as max number of non-empty genotypes.
n_cont_eur <- df_eur$non_miss_cont %>% max()  # Same for controls.
n_case_afr <- df_afr$non_miss_case %>% max()
n_cont_afr <- df_afr$non_miss_cont %>% max()
n_case_amr <- df_amr$non_miss_case %>% max()
n_cont_amr <- df_amr$non_miss_cont %>% max()

# Calculate MAC, MAF, variant_call_rate and hwe chisq pvalue
df_eur <- prepare_data(df_eur, n_case_eur, n_cont_eur)
df_afr <- prepare_data(df_afr, n_case_afr, n_cont_afr)
df_amr <- prepare_data(df_amr, n_case_amr, n_cont_amr)

# "star"-mutation (*) filter.
drop <- as_data_frame(fread('data/drop.tsv')) %>%  # obtained from cases.vcf
  dplyr::mutate(chr_pos=paste(V1,V2,sep=':'))

df_eur <- df_eur %>% dplyr::mutate(chr_pos=paste(chr,pos,sep=':')) %>% dplyr::filter(!chr_pos %in% drop$chr_pos)
df_afr <- df_afr %>% dplyr::mutate(chr_pos=paste(chr,pos,sep=':')) %>% dplyr::filter(!chr_pos %in% drop$chr_pos)
df_amr <- df_amr %>% dplyr::mutate(chr_pos=paste(chr,pos,sep=':')) %>% dplyr::filter(!chr_pos %in% drop$chr_pos)

# Apply corresponding filters
df_eur <- apply_filters(df_eur, call_rate_n=0.9, call_rate_s=0.9, hwe_pvalue_n=10^-7, hwe_pvalue_s=10^-7)
df_afr <- apply_filters(df_afr, call_rate_n=0.9, call_rate_s=0.9, hwe_pvalue_n=10^-7, hwe_pvalue_s=10^-7)
df_amr <- apply_filters(df_amr, call_rate_n=0.9, call_rate_s=0.9, hwe_pvalue_n=10^-7, hwe_pvalue_s=10^-7)

# Confirm matching quality
matching_confirm_eur <- df_eur %>% CVAS(., gnomad_nfe_af=0.01, conseq_filter=c('synonymous_variant'))
matching_confirm_afr <- df_afr %>% CVAS(., gnomad_afr_af=0.01, conseq_filter=c('synonymous_variant'))
matching_confirm_amr <- df_amr %>% CVAS(., gnomad_amr_af=0.01, conseq_filter=c('synonymous_variant'))

p_eur_conf <- lambdaGC(matching_confirm_eur$pvalue, plot=T,
         filter_massage = 'Matching confirmation\n\tgnomAD NFE AF >= 0.01\n\tconsequence: synonymous variants\n\ttest: linear regression', text_coord = c(0, 3))
p_afr_conf <- lambdaGC(matching_confirm_afr$pvalue, plot=T,
         filter_massage = 'Matching confirmation\n\tgnomAD AFR AF >= 0.01\n\tconsequence: synonymous variants\n\ttest: linear regression', text_coord = c(0, 3))
p_amr_conf <- lambdaGC(matching_confirm_amr$pvalue, plot=T,
                     filter_massage = 'Matching confirmation\n\tgnomAD AMR AF >= 0.01\n\tconsequence: synonymous variants\n\ttest: linear regression', text_coord = c(0, 3))
(p_eur_conf|p_afr_conf|p_amr_conf)

###########################################################
### CVAS, missense variants + LoF, gnomAF_POP_AF > 0.01 ###
###########################################################

### EUR CVAS ###
cvas_eur <- df_eur %>% CVAS(., gnomad_nfe_af=0.01,
                            conseq_filter=c('missense_variant', 'stop_gained', 'frameshift_variant',
                                            'splice_acceptor_variant', 'splice_donor_variant'))
cvas_eur %>% dplyr::select(chr, pos, pvalue,  gene_symbol, rsid) %>% arrange(pvalue) %>%  as.data.frame() %>% head(50)
lbl_eur <- cvas_eur %>% add_rownames() %>%
  filter(rowname %in% c(1)) %>% mutate(rsid = str_extract(rsid, "[^\\&$]+")) %>%
  unite(lbl, c('rsid', 'gene_symbol')) %>% dplyr::select(rowname, lbl)

p_cvas_eur <- lambdaGC(cvas_eur$pvalue, plot=T, lbl_eur,
         filter_massage = 'gnomAD NFE AF > 0.01\n\tconsequence: missense variants + LoF', text_coord = c(0, 3))

write.table(cvas_eur, 'data/assoc_results/eur.cvas.tsv',
            quote = F, row.names = F, sep='\t')

### AFR CVAS ###
cvas_afr <- df_afr %>%
  CVAS(., gnomad_afr_af=0.01,
       conseq_filter=c('missense_variant', 'stop_gained', 'frameshift_variant','splice_acceptor_variant', 'splice_donor_variant'))
cvas_afr %>% dplyr::select(chr, pos, pvalue,  gene_symbol, rsid) %>% arrange(pvalue) %>%  as.data.frame() %>% head(50)

lbl_afr <- cvas_afr %>% add_rownames() %>% filter(rowname %in% c()) %>%
  mutate(rsid = str_extract(rsid, "rs[^&]+")) %>%
  unite(lbl, c('rsid', 'gene_symbol')) %>% dplyr::select(rowname, lbl)
p_cvas_afr <- lambdaGC(cvas_afr$pvalue, plot=T, lbl_afr,
         filter_massage = 'CVAS\n\tgnomAD AFR AF >= 0.01\n\tconsequence: missense variants + LoF\n\ttest: linear regression', text_coord = c(0, 3))

write.table(cvas_afr, 'data/assoc_results/afr.cvas.tsv', quote = F, row.names = F, sep='\t')

# ### AMR CVAS ###
cvas_amr <- df_amr %>%
  CVAS(., gnomad_amr_af=0.01,
       conseq_filter=c('missense_variant', 'stop_gained', 'frameshift_variant','splice_acceptor_variant', 'splice_donor_variant'))
cvas_amr %>% dplyr::select(chr, pos, pvalue,  gene_symbol, rsid) %>% arrange(pvalue) %>%  as.data.frame() %>% head(50)

lbl_amr <- cvas_amr %>% add_rownames() %>% filter(rowname %in% c()) %>%
  mutate(rsid = str_extract(rsid, "rs[^&]+")) %>%
  unite(lbl, c('rsid', 'gene_symbol')) %>% dplyr::select(rowname, lbl)
lambdaGC(cvas_amr$pvalue, plot=T,
         filter_massage = 'CVAS\n\tgnomAD AMR AF >= 0.01\n\tconsequence: missense variants + LoF\n\ttest: linear regression', text_coord = c(0, 3))
write.table(cvas_amr, 'data/assoc_results/amr.cvas.tsv', quote = F, row.names = F, sep='\t')

###########################################################
### RVAS, missense variants + LoF, gnomAD_POP_AF < 0.01 ###
###########################################################

### EUR RVAS ###
rvas_eur <- df_eur %>% 
  RVAS(., n_case.=n_case_eur, n_cont.=n_cont_eur,
                   gnomad_nfe_af=0.01, conseq_filter=c('missense_variant', 'stop_gained', 'frameshift_variant',
                                                       'splice_acceptor_variant', 'splice_donor_variant'))

rvas_eur <- as_data_frame(fread('data/assoc_results/eur.rvas.tsv'))

rvas_eur %>% head(15) %>%
  dplyr::select(gene_symbol, Simes) %>% arrange(Simes) %>% as.data.frame()
  

lbl_eur <- rvas_eur %>% arrange(Simes) %>% add_rownames() %>%
  dplyr::filter(gene_symbol %in% c('APOL1', 'CR1', 'KANK1', 'IL36G', 'COL4A4')) %>% 
  mutate(lbl=gene_symbol) %>% dplyr::select(rowname, lbl)

lambdaGC(rvas_eur$Simes %>% na.omit(), plot = T, lbl_eur,
         filter_massage = 'RVAS\n\tgnomAD EUR AF < 0.01\n\tconsequence: missense variants + LoF\n\ttest: Simes\n\t(Fisher\'s exact test\n\tC-alpha, ASUM, WSS, KBAC)',
         hline = -log10(0.05/n_genes),
         text_coord = c(0, 7))

lbl_eur <- rvas_eur %>% arrange(`C-alpha`) %>% add_rownames() %>%
  filter(gene_symbol %in% c('APOL1', 'CR1', 'KANK1', 'IL36G', 'COL4A4')) %>%
  mutate(lbl=gene_symbol) %>% dplyr::select(rowname, lbl)
p_calpha <- lambdaGC(rvas_eur$`C-alpha` %>% na.omit(), plot = T, lbl_eur,
         hline = -log10(0.05/n_genes),
         text_coord = c(0, 7))

lbl_eur <- rvas_eur %>% arrange(fisher_carriers.pvalue) %>% add_rownames() %>%
  filter(gene_symbol %in% c('APOL1', 'CR1', 'KANK1', 'IL36G', 'COL4A4')) %>%
  mutate(lbl=gene_symbol) %>% dplyr::select(rowname, lbl)
p_fisher <- lambdaGC(rvas_eur$fisher_carriers.pvalue %>% na.omit(), plot = T, lbl_eur,
                     hline = -log10(0.05/n_genes),
                     text_coord = c(0, 7))

lbl_eur <- rvas_eur %>% arrange(WSS) %>% add_rownames() %>%
  filter(gene_symbol %in% c('APOL1', 'CR1', 'KANK1', 'IL36G', 'COL4A4')) %>%
  mutate(lbl=gene_symbol) %>% dplyr::select(rowname, lbl)
p_wss <- lambdaGC(rvas_eur$WSS %>% na.omit(), plot = T, lbl_eur,
                     hline = -log10(0.05/n_genes),
                     text_coord = c(0, 7))

lbl_eur <- rvas_eur %>% arrange(KBAC) %>% add_rownames() %>%
  filter(gene_symbol %in% c('APOL1', 'CR1', 'KANK1', 'IL36G', 'COL4A4')) %>%
  mutate(lbl=gene_symbol) %>% dplyr::select(rowname, lbl)
p_kbac <- lambdaGC(rvas_eur$KBAC %>% na.omit(), plot = T, lbl_eur,
                  hline = -log10(0.05/2482),
                  text_coord = c(0, 7))

lbl_eur <- rvas_eur %>% arrange(ASUM) %>% add_rownames() %>%
  filter(gene_symbol %in% c('APOL1', 'CR1', 'KANK1', 'IL36G', 'COL4A4')) %>%
  mutate(lbl=gene_symbol) %>% dplyr::select(rowname, lbl)
p_asum <- lambdaGC(rvas_eur$ASUM %>% na.omit(), plot = T, lbl_eur,
                   hline = -log10(0.05/2482),
                   text_coord = c(0, 7))

(p_calpha|p_fisher)/(p_wss|p_kbac|p_asum)+plot_annotation(tag_levels = 'A')


rvas_eur %>% write.table('data/assoc_results/eur.rvas.tsv',
            quote = F, row.names = F, sep='\t')

### AFR RVAS ###
rvas_afr <- df_afr %>%
  RVAS(., n_case.=n_case_afr, n_cont.=n_cont_afr,
       gnomad_afr_af=0.01, conseq_filter=c('missense_variant', 'stop_gained', 'frameshift_variant',
                                           'splice_acceptor_variant', 'splice_donor_variant'))

rvas_afr %>% dplyr::select(gene_symbol, Simes) %>% arrange(Simes) %>%  as.data.frame() %>% head(15)

lbl_afr <- rvas_afr %>% arrange(linear_regression) %>% add_rownames() %>%
  filter(rowname %in% c()) %>%
  mutate(lbl=gene_symbol) %>% dplyr::select(rowname, lbl)

lambdaGC(rvas_afr$Simes %>% na.omit(), plot = T, lbl_afr,
         filter_massage = 'RVAS\n\tgnomAD AFR AF < 0.01\n\tconsequence: missense variants + LoF\n\ttest: Simes\n\t(Fisher\'s exact test\n\tC-alpha, ASUM, WSS, KBAC)',
         text_coord = c(0, 7))

write.table(rvas_afr, 'data/assoc_results/afr.rvas.tsv',
            quote = F, row.names = F, sep='\t')

### AMR RVAS ###
rvas_amr <- df_amr %>%
  RVAS(., n_case.=n_case_amr, n_cont.=n_cont_amr,
       gnomad_amr_af=0.01, conseq_filter=c('missense_variant', 'stop_gained', 'frameshift_variant',
                                           'splice_acceptor_variant', 'splice_donor_variant'))

rvas_amr %>% dplyr::select(gene_symbol, Simes) %>% arrange(Simes) %>% as.data.frame() %>% head(50)

lbl_amr <- rvas_amr %>% arrange(linear_regression) %>% add_rownames() %>%
  filter(rowname %in% c()) %>%
  mutate(lbl=gene_symbol) %>% dplyr::select(rowname, lbl)

lambdaGC(rvas_amr$Simes %>% na.omit(), plot = T, lbl_amr,
         filter_massage = 'RVAS\n\tgnomAD AMR AF < 0.01\n\tconsequence: missense variants + LoF\n\ttest: Simes\n\t(Fisher\'s exact test\n\tC-alpha, ASUM, WSS, KBAC)',
         text_coord = c(0, 7))

write.table(rvas_amr, 'data/assoc_results/amr.rvas.tsv',
            quote = F, row.names = F, sep='\t')
