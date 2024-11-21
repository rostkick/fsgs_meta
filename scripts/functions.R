#### HWE ####
hwe_chisq <- function(AA_obs, Aa_obs, aa_obs) {
  sum_alleles=AA_obs+Aa_obs+aa_obs
  AF=(aa_obs+0.5*Aa_obs)/sum_alleles
  
  AA_exp=(1-AF)^2*sum_alleles
  Aa_exp=2*AF*(1-AF)*sum_alleles
  aa_exp=AF^2*sum_alleles
  
  hwe_chisq=((AA_obs-AA_exp)^2/AA_exp)+((Aa_obs-Aa_exp)^2/Aa_exp)+((aa_obs-aa_exp)^2/aa_exp)
  pvalue <- pchisq(hwe_chisq, df=1, lower.tail = F)
  return(pvalue)
}


#### CVAS linear regression ####
l_reg <- function(cases, controls){
  #` Linear regression function
  table <- data.frame(geno = c(0, 1, 2, 0, 1, 2), phe = c(0, 0, 0, 1, 1, 1))
  ws <- c(controls, cases)
  m <- lm(phe ~ geno, data = table, weights = ws)
  m$df.residual <- sum(ws) - 2
  ks <- summary(m)$coefficients
  list(beta = ks[2, 1], pvalue = ks[2, 4])
}

CVAS <- function(df,
                 gnomad_nfe_af=0, gnomad_afr_af=0, gnomad_amr_af=0, gnomad_eas_af=0, gnomad_sas_af=0,
                 conseq_filter=c('missense_variant')){
  result <- df %>% rowwise() %>%
      #` Apply filters
      dplyr::filter(gnomAD_NFE_AF >= gnomad_nfe_af,
                    gnomAD_AFR_AF >= gnomad_afr_af,
                    gnomAD_AMR_AF >= gnomad_amr_af,
                    gnomAD_EAS_AF >= gnomad_eas_af,
                    gnomAD_SAS_AF >= gnomad_sas_af,
                    Consequence %in% conseq_filter) %>% 
    do(data_frame(chr=.$chr, pos=.$pos, rsid=.$rsid, gene_symbol=.$gene_symbol,
                  maf=.$maf,
                  gnomAD_NFE_AF=.$gnomAD_NFE_AF,
                  gnomAD_AFR_AF=.$gnomAD_AFR_AF,
                  gnomAD_AMR_AF=.$gnomAD_AMR_AF,
                  gnomAD_EAS_AF=.$gnomAD_EAS_AF,
                  gnomAD_SAS_AF=.$gnomAD_SAS_AF,
                  hom_ref_case=.$hom_ref_case, het_case=.$het_case, hom_var_case=.$hom_var_case,
                  hom_ref_cont=.$hom_ref_cont, het_cont=.$het_cont, hom_var_cont=.$hom_var_cont,
                  call_rate_case=.$call_rate_case, call_rate_cont=.$call_rate_cont,
                  as_data_frame(l_reg(c(.$hom_ref_case, .$het_case, .$hom_var_case),
                                      c(.$hom_ref_cont, .$het_cont, .$hom_var_cont))))) %>% 
  arrange(pvalue)
  return(result)
}

#### RVAS functions ####
prep_tmp <- function(df){
  names(df) <- gsub(x = names(df), pattern = "control", replacement = "cont")
  names(df) <- gsub(x = names(df), pattern = "\\.", replacement = "_")
  df <- df %>% dplyr::rename(chr=`#CHROM`, pos=POS, rsid=ID, ref=REF, alt=ALT, gene_symbol=SYMBOL)
  return(df)
}

prepare_data <- function(df, n_s, n_n){
  df <- df %>% 
    mutate(gnomAD_AF=replace_na(gnomAD_AF, 0),
           gnomAD_NFE_AF=replace_na(gnomAD_NFE_AF, 0),
           gnomAD_AFR_AF=replace_na(gnomAD_AFR_AF, 0),
           gnomAD_AMR_AF=replace_na(gnomAD_AMR_AF, 0),
           gnomAD_EAS_AF=replace_na(gnomAD_EAS_AF, 0),
           gnomAD_SAS_AF=replace_na(gnomAD_SAS_AF, 0)) %>% 
    mutate(minor_hom_case=if_else(hom_ref_case < hom_var_case, hom_ref_case, hom_var_case),
           minor_hom_cont=if_else(hom_ref_cont < hom_var_cont, hom_ref_cont, hom_var_cont),
           major_hom_case=if_else(hom_ref_case > hom_var_case, hom_ref_case, hom_var_case),
           major_hom_cont=if_else(hom_ref_cont > hom_var_cont, hom_ref_cont, hom_var_cont),
           MAC_case=het_case+minor_hom_case*2,
           MAC_cont=het_cont+minor_hom_cont*2,
           maf=((het_cont + het_case) + (hom_var_case + hom_var_cont)*2)/
             ((hom_ref_cont+hom_ref_case)*2+(het_cont+het_case)+(hom_var_cont+hom_var_case)*2),
           maf_case=pmin((hom_var_case*2+het_case)/(2*(hom_ref_case+het_case+hom_var_case)),
                         (hom_ref_case*2+het_case)/(2*(hom_ref_case+het_case+hom_var_case))),
           maf_cont=pmin((hom_var_cont*2+het_cont)/(2*(hom_ref_cont+het_cont+hom_var_cont)),
                        (hom_ref_cont*2+het_cont)/(2*(hom_ref_cont+het_cont+hom_var_cont))),
           call_rate_cont = (major_hom_cont+minor_hom_cont+het_cont)/n_n,
           call_rate_case= (major_hom_case+minor_hom_case+het_case)/n_s,
           hwe_pvalue_cont = hwe_chisq(hom_ref_cont, het_cont, hom_var_cont),
           hwe_pvalue_case = hwe_chisq(hom_ref_case, het_case, hom_var_case))
  return(df)
}

apply_filters <- function(df, call_rate_n=0.9, call_rate_s=0.9, hwe_pvalue_n=10^-4, hwe_pvalue_s=10^-4){
  df <- df %>%
    # Apply variant call rate filter for cases and controls separatly.
    dplyr::filter(call_rate_cont > call_rate_n &
             call_rate_case > call_rate_s) %>%
    # Apply HWE filter
    dplyr::filter(!(is.na(hwe_pvalue_n) & is.na(hwe_pvalue_s)) & (hwe_pvalue_case > hwe_pvalue_s & hwe_pvalue_cont > hwe_pvalue_n)) %>% 
    dplyr::filter(MAC_case + MAC_cont > 0)  # at least 1 alternative allele
  return(df)
}

fisher <- function(a,b,c,d){
  #' Wrapper function based on stats::fisher.test()
  data <- matrix(c(a,b,c,d),ncol=2)
  f <- fisher.test(data)
  return(f)
}

fisher_carriers <- function(df, n_case, n_cont){
  #` Fisher's exact test on carriers contingency table
  df %>%
    group_by(gene_symbol) %>%
    summarise(hom_ref_cont=sum(major_hom_cont), het_cont=sum(het_cont), hom_var_cont=sum(minor_hom_cont),
              hom_ref_case=sum(major_hom_case), het_case=sum(het_case), hom_var_case=sum(minor_hom_case)) %>% 
    mutate(carriers_in_cont=het_cont+hom_var_cont,
           carriers_in_case=het_case+hom_var_case) %>% 
    group_by(gene_symbol) %>%
    summarise(n=sum(carriers_in_cont), s=sum(carriers_in_case)) %>%
    filter(!(n > n_cont | s > n_case)) %>%
    rowwise() %>%
    dplyr::mutate(fisher_carriers.pvalue=fisher(s, n_case-s, n, n_cont-n)$p.value,
                  fisher_carriers.or=fisher(s, n_case-s, n, n_cont-n)$estimate) %>% 
    dplyr::select(gene_symbol, fisher_carriers.pvalue, fisher_carriers.or)
}

fisher_carriers_or <- function(df, n_case, n_cont){
  #` Fisher's exact test on carriers contingency table
  df %>%
    group_by(gene_symbol) %>%
    summarise(hom_ref_cont=sum(major_hom_cont), het_cont=sum(het_cont), hom_var_cont=sum(minor_hom_cont),
              hom_ref_case=sum(major_hom_case), het_case=sum(het_case), hom_var_case=sum(minor_hom_case)) %>% 
    mutate(carriers_in_cont=het_cont+hom_var_cont,
           carriers_in_case=het_case+hom_var_case) %>% 
    group_by(gene_symbol) %>%
    summarise(n=sum(carriers_in_cont), s=sum(carriers_in_case)) %>%
    filter(!(n > n_cont | s > n_case)) %>%
    rowwise() %>%
    dplyr::mutate(fisher_carriers_or=fisher(s, n_case-s, n, n_cont-n)$estimate) %>% 
    dplyr::select(gene_symbol, fisher_carriers_or)
}

agg_df <- function(df){
  #` Summarize genotypes counts for Fisher's alleles test and linear regression 
  df %>%
    dplyr::group_by(gene_symbol) %>%
    summarise(hom_ref_cont=sum(hom_ref_cont), het_cont=sum(het_cont), hom_var_cont=sum(hom_var_cont),
              hom_ref_case=sum(hom_ref_case), het_case=sum(het_case), hom_var_case=sum(hom_var_case)) %>% 
    rowwise()
}

lin_reg_rvas <- function(df){
  df %>% agg_df() %>% 
    rowwise() %>%
    do(data.frame(gene_symbol=.$gene_symbol, linear_regression=l_reg(c(.$hom_ref_case, .$het_case, .$hom_var_case),
                                                                     c(.$hom_ref_cont, .$het_cont, .$hom_var_cont))$pval))
}

fisher_alleles <- function(df){
  #` Fisher's exact test on allele contingency table
  df %>% agg_df() %>% 
    mutate(fisher_alleles=fisher(het_case+hom_var_case*2, hom_ref_case*2+het_case,
                                 het_cont+hom_var_cont*2, hom_ref_cont*2+het_cont)$p.value) %>%
    dplyr::select(gene_symbol, fisher_alleles)
}

c_alpha <- function(df, n_case, n_cont, sim=100){
  #` Wrapper function based on AssotesteR::CALPHA()
  y <- c(rep(1, n_case), rep(0, n_cont))
  df <- df  %>% dplyr::select(hom_var_case, het_case, hom_var_cont, het_cont)
  if(nrow(df)>1){
    X <- apply(df, 1, function(x){c(sample(c(rep(2, x[1]),
                                             rep(1, x[2]),
                                             rep(0, n_case-(x[1]+x[2]))
    ), n_case),
    sample(c(rep(2, x[3]),
             rep(1, x[4]),
             rep(0, n_cont-(x[3]+x[4]))
    ), n_cont))
    })
    table_query <- cbind(y,X)
    # result=calpha(table_query, perm=sim)$perm.pval
    result=CALPHA(y, X, perm=sim)$perm.pval
  } else {
    result=1
  }
  return(data.frame(c_alpha=result))
}

kbac_p <- function(df, n_case, n_cont, alpha=NULL, sim=100){
  #` Wrapper function based on vartools::kbac()
  y <- c(rep(1, n_case), rep(0, n_cont))
  df <- df  %>% dplyr::select(hom_var_case, het_case, hom_var_cont, het_cont)
  if(nrow(df)>1){
    X <- apply(df, 1, function(x){c(sample(c(rep(2, x[1]),
                                             rep(1, x[2]),
                                             rep(0, n_case-(x[1]+x[2]))
    ), n_case),
    sample(c(rep(2, x[3]),
             rep(1, x[4]),
             rep(0, n_cont-(x[3]+x[4]))
    ), n_cont))
    })
    table_query <- cbind(y,X)
    result=kbac(table_query, alpha=alpha, num.permutation=sim)$perm.pval
  } else {
    result=1
  }
  return(data.frame(kbac_p=result))
}

wss <- function(df, n_case, n_cont, sim=100){
  #` Wrapper function based on AssotesteR::WSS()
  y <- c(rep(1, n_case), rep(0, n_cont))
  df <- df %>% dplyr::select(hom_var_case, het_case,
                                                                         hom_var_cont, het_cont)
  if(nrow(df)>1){
    X <- apply(df, 1, function(x){c(sample(c(rep(2, x[1]),
                                             rep(1, x[2]),
                                             rep(0, n_case-(x[1]+x[2]))
    ), n_case),
    sample(c(rep(2, x[3]),
             rep(1, x[4]),
             rep(0, n_cont-(x[3]+x[4]))
    ), n_cont))
    })
    result=WSS(y, X, perm=sim)$perm.pval
  } else {
    result=1
    
  }
  return(data.frame(wss=result))
}

asum <- function(df, n_case, n_cont, sim=100){
  #` Wrapper function based on AssotesteR::ASUM()
  y <- c(rep(1, n_case), rep(0, n_cont))
  df <- df %>% dplyr::select(hom_var_case, het_case,hom_var_cont, het_cont)
  
  if(nrow(df)>1){
    X <- apply(df, 1, function(x){c(sample(c(rep(2, x[1]),
                                             rep(1, x[2]),
                                             rep(0, n_case-(x[1]+x[2]))
    ), n_case),
    sample(c(rep(2, x[3]),
             rep(1, x[4]),
             rep(0, n_cont-(x[3]+x[4]))
    ), n_cont))
    })
    result=ASUM(y, X, perm=sim)$perm.pval
  } else {
    result=1
  }
  return(data.frame(asum=result))
}

RVAS <- function(df, n_case.=n_case, n_cont.=n_cont,
                 gnomad_nfe_af=10, gnomad_afr_af=10, gnomad_amr_af=10, gnomad_eas_af=10, gnomad_sas_af=10,
                 conseq_filter=c('missense_variant')) {
  #` Apply filters
  df <- df %>% dplyr::filter(gnomAD_NFE_AF < gnomad_nfe_af,
                             gnomAD_AFR_AF < gnomad_afr_af,
                             gnomAD_AMR_AF < gnomad_amr_af,
                             gnomAD_EAS_AF < gnomad_eas_af,
                             gnomAD_SAS_AF < gnomad_sas_af,Consequence %in% conseq_filter)
  #` Apply RVAS-tests and outer_join the results
  result <- df %>% nest_by(gene_symbol, .keep=T) %>% summarise(mac_case_sum=sum(data$MAC_case),
                                                               mac_cont_sum=sum(data$MAC_cont),
                                                               maf_case=mac_case_sum/(n_case_eur*2),
                                                               maf_cont=mac_cont_sum/(n_cont_eur*2),
                                                               lin_reg_rvas(data),
                                                               fisher_carriers(data, n_case., n_cont.),
                                                              #  fisher_carriers_or(data, n_case., n_cont.),
                                                               fisher_alleles(data),
                                                               c_alpha(data, n_case., n_cont., sim=100),
                                                               kbac_p(data, n_case., n_cont., sim=100),
                                                               asum(data, n_case., n_cont., sim=100),
                                                               wss(data, n_case., n_cont., sim=100))

  # CALPHA
  genes_zero <- result[which(result$c_alpha==0),]$gene_symbol
  if (length(genes_zero)>0){
    print('recalculate p.values (test: C-alpha, sim: 1000) for: ')
    print(genes_zero)
    r_upd <- df %>% filter(gene_symbol %in% genes_zero) %>% nest_by(gene_symbol) %>%
      summarise(c_alpha(data, sim=1000, n_case., n_cont.))
    result[which(result$c_alpha==0),]$c_alpha <- r_upd$c_alpha
    genes_zero <- result[which(result$c_alpha==0),]$gene_symbol
  }
  if (length(genes_zero)>0){
    print('recalculate p.values (test: C-alpha, sim: 10000) for: ')
    print(genes_zero)
    r_upd <- df %>% filter(gene_symbol %in% genes_zero) %>% nest_by(gene_symbol) %>%
      summarise(c_alpha(data, sim=10000, n_case., n_cont.))
    result[which(result$c_alpha==0),]$c_alpha <- r_upd$c_alpha
  }
  
  # WSS
  genes_zero <- result[which(result$wss==0),]$gene_symbol
  if (length(genes_zero) > 0){
    print('recalculate p.values (test: WSS, sim: 1000) for: ')
    print(genes_zero)
    r_upd <- df %>% filter(gene_symbol %in% genes_zero) %>% nest_by(gene_symbol) %>%
      summarise(wss(data, sim=1000, n_case., n_cont.))
    result[which(result$wss==0),]$wss <- r_upd$wss
    genes_zero <- result[which(result$wss==0),]$gene_symbol
  }
  if (length(genes_zero) > 0){
    print('recalculate p.values (test: WSS, sim: 10000) for: ')
    print(genes_zero)
    r_upd <- df %>% filter(gene_symbol %in% genes_zero) %>% nest_by(gene_symbol) %>%
      summarise(wss(data, sim=10000, n_case., n_cont.))
    result[which(result$wss==0),]$wss <- r_upd$wss 
  }
  
  # ASUM
  genes_zero <- result[which(result$asum==0),]$gene_symbol
  if (length(genes_zero) > 0){
    print('recalculate p.values (test: ASUM, sim: 1000) for: ')
    print(genes_zero)
    r_upd <- df %>% filter(gene_symbol %in% genes_zero) %>% nest_by(gene_symbol) %>%
      summarise(asum(data, n_case., n_cont., sim=1000))
    result[which(result$asum==0),]$asum <- r_upd$asum
    genes_zero <- result[which(result$asum==0),]$gene_symbol
  }
  if (length(genes_zero) > 0){
    print('recalculate p.values (test: ASUM, sim: 10000) for: ')
    print(genes_zero)
    r_upd <- df %>% filter(gene_symbol %in% genes_zero) %>% nest_by(gene_symbol) %>%
      summarise(asum(data, sim=10000, n_case., n_cont.))
    result[which(result$asum==0),]$asum <- r_upd$asum
  }
  
  print('Combining p.values by Simes method')
  result <- result %>% mutate(c_alpha = replace(c_alpha, c_alpha==0, 0.00001),
                              asum = replace(asum, asum==0, 0.00001),
                              kbac_p = replace(kbac_p, kbac_p==0, 0.00001),
                              wss = replace(wss, wss==0, 0.00001)) %>%
    mutate(simes=combinePValues(fisher_carriers.pvalue, c_alpha, asum,
                                kbac_p,
                                wss, method = 'simes')) %>%
    arrange(simes) %>% 
    dplyr::rename("C-alpha"=c_alpha, ASUM=asum,
                  KBAC=kbac_p,
                  WSS=wss, Simes=simes)
  return(result)
}

 mantelhaen <- function(df){
  genes <- df$gene_symbol %>% unique()
  pvalues <- list()
  for (g in genes){
    tt <- dplyr::filter(df, gene_symbol==g)
    cc <- c()
    for (i in 1:nrow(tt)){
      cc <- c(cc, tt[[i, 'ref_case']], tt[[i, 'alt_case']], tt[[i, 'ref_cont']], tt[[i, 'alt_cont']])
    }
    print(g)
    query <- array(cc, dim = c(2, 2, 3), dimnames = list(gt=c('ref', 'alt'), fsgs=c('case', 'cont'), pop=c('AMR', 'EAS', 'SAS')))
    pvalues[g] <- mantelhaen.test(query)$p.value
  }
  pvalues <- pvalues %>% unlist() %>% data.frame(mantelhaen.pvalues=.)
  pvalues$gene_symbol <- rownames(pvalues)
  pvalues <- pvalues %>% as_data_frame() %>% dplyr::select(gene_symbol, mantelhaen.pvalues)
  return(pvalues)
}

#### QQplot function ####
lambdaGC <- function (data, plot = FALSE, hline=-log10(0.05/2500),
                      proportion = 1, method = "regression", filter = TRUE, df = 1, exp=NA, ...){
  #` Reworked GenABEL::estlambda(). Only plot function was changing.
  pvals <- copy(data)
  len_pvals <- length(pvals)
  print(len_pvals)
  data <- data[which(!is.na(data))]
  if (proportion > 1 || proportion <= 0) 
    stop("proportion argument should be greater then zero and less than or equal to one")
  ntp <- round(proportion * length(data))
  if (ntp < 1) 
    stop("no valid measurements")
  if (ntp == 1) {
    warning(paste("One measurement, lambda = 1 returned"))
    return(list(estimate = 1, se = 999.99))
  }
  if (ntp < 10)
    warning(paste("number of points is too small:", ntp))
  if (min(data) < 0) 
    stop("data argument has values <0")
  if (max(data) <= 1) {
    data <- qchisq(data, 1, lower.tail = FALSE)
  }
  if (filter) {
    data[which(abs(data) < 1e-08)] <- NA
  }
  data <- sort(data)
  ppoi <- ppoints(data)
  ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
  data <- data[1:ntp]
  ppoi <- ppoi[1:ntp]
  out <- list()
  
  if(is.na(exp)) exp=-log10(1:len_pvals/len_pvals)

  dd <- data_frame(exp=exp, obs=-log10(sort(pvals)))
  if (method == "regression") {
    s <- summary(lm(data ~ 0 + ppoi))$coeff
    out$estimate <- s[1, 1]
    out$se <- s[1, 2]
  }
  else if (method == "median") {
    out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, df)
    out$se <- NA
  }
  else if (method == "KS") {
    limits <- c(0.5, 100)
    out$estimate <- estLambdaKS(data, limits = limits, df = df)
    if (abs(out$estimate - limits[1]) < 1e-04 || abs(out$estimate - 
                                                     limits[2]) < 1e-04) 
      warning("using method='KS' lambda too close to limits, use other method")
    out$se <- NA
  }
  else {
    stop("'method' should be either 'regression' or 'median'!")
  }
  if (plot) {
    dd$gene_symbol <- ''
    p <- ggplot(dd, aes(x=exp, y=obs))+
      geom_point(shape=19, size=4, alpha=0.8, color="#3C5488FF")+
      geom_hline(yintercept=hline, linetype='dashed', color='#DC0000FF')+
      geom_abline(intercept = 0, slope = 1, alpha = 0.8, color='black')+
      scale_x_continuous(expression(paste("Expected -log"[10], plain(P))))+
      scale_y_continuous(expression(paste("Observed -log"[10], plain(P))))+
      theme_classic()
    show(p)
  }
  print(out)
  return(p)
}

#### Manhattan plot ####
manhattan_plot <- function(gwas.dat) {
  # Manhattan plot function.
  gwas.dat <- gwas.dat %>% rename(CHR=chr, BP=pos, P=pval) %>%
    filter(!CHR %in% c(23, 24))

  sig.dat <- gwas.dat %>% 
    subset(P < 0.05)
  notsig.dat <- gwas.dat %>% 
    subset(P >= 0.05) %>%
    slice(sample(nrow(.), nrow(.) / 5))
  gwas.dat <- rbind(sig.dat,notsig.dat)
  nCHR <- length(unique(gwas.dat$CHR))
  gwas.dat$BPcum <- 0
  s <- 0
  nbp <- c()
  for (i in unique(gwas.dat$CHR)){
    nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$BP)
    gwas.dat[gwas.dat$CHR == i,"BPcum"] <- gwas.dat[gwas.dat$CHR == i,"BP"] + s
    s <- s + nbp[i]
  }

  axis.set <- gwas.dat %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ylim <- abs(floor(log10(min(gwas.dat$P)))) + 2 
  sig <- nrow(gwas.dat)
  manhplot <- ggplot(gwas.dat, aes(x = BPcum, y = -log10(P), 
                                   color = as.factor(CHR), size = -log10(P))) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = -log10(0.05/sig), color = "grey40", linetype = "dashed") + 
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#3C5488FF", "#DC0000FF"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log10(p)") + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
    )
  return(manhplot)
}

manhattan_mirror <- function (top, bottom, tline, bline, log10 = TRUE, yaxis, opacity = 1, 
           annotate_snp, annotate_p, toptitle = NULL, bottomtitle = NULL, 
           highlight_snp, highlight_p, highlighter = "red", chrcolor1 = "#AAAAAA", 
           chrcolor2 = "#4D4D4D", freey = FALSE, background = "variegated", 
           chrblocks = FALSE, file = "gmirror", hgt = 7, hgtratio = 0.5, 
           wi = 12, res = 300) 
 {
   topn <- names(top)
   bottomn <- names(bottom)
   top$Location <- "Top"
   bottom$Location <- "Bottom"
   d <- rbind(top, bottom)
   d$POS <- as.numeric(as.character(d$POS))
   d$CHR <- factor(d$CHR, levels = c("1", "2", "3", "4", "5", 
                                     "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 
                                     "16", "17", "18", "19", "20", "21", "22", "X", "Y"))
   d_order <- d[order(d$CHR, d$POS), ]
   d_order$pos_index <- seq.int(nrow(d_order))
   d_order_sub <- d_order[, c("SNP", "CHR", "POS", "pvalue", 
                              "pos_index")]
   maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$pos_index), 
                                                             ])
   minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$pos_index), 
                                                             ])
   milimits <- do.call(rbind, minRows)
   malimits <- do.call(rbind, maxRows)
   lims <- merge(milimits, malimits, by = "CHR")
   names(lims) <- c("Color", "snpx", "px", "posx", "posmin", 
                    "snpy", "py", "posy", "posmax")
   lims$av <- (lims$posmin + lims$posmax)/2
   lims <- lims[order(lims$Color), ]
   lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), length.out = nrow(lims), 
                        each = 1)
   nchrcolors <- nlevels(factor(lims$Color))
   colnames(d_order)[2] <- "Color"
   newcols <- c(rep(x = c(chrcolor1, chrcolor2), length.out = nchrcolors, 
                    each = 1), "#FFFFFF", "#EBEBEB")
   names(newcols) <- c(levels(factor(lims$Color)), "shade_ffffff", 
                       "shade_ebebeb")
   if (log10 == TRUE) {
     d_order$pval <- -log10(d_order$pvalue)
     yaxislab1 <- expression(paste("-log"[10], "(p-value)", 
                                   sep = ""))
     yaxislab2 <- expression(paste("-log"[10], "(p-value)", 
                                   sep = ""))
     if (!missing(tline)) {
       tredline <- -log10(tline)
     }
     if (!missing(bline)) {
       bredline <- -log10(bline)
     }
   }
   else {
     d_order$pval <- d_order$pvalue
     yaxislab1 <- yaxis[1]
     yaxislab2 <- yaxis[2]
     if (!missing(tline)) {
       tredline <- tline
     }
     if (!missing(bline)) {
       bredline <- bline
     }
   }
   yaxismax1 <- ifelse(freey == FALSE, max(d_order$pval[which(d_order$pval < 
                                                                Inf)]), max(d_order$pval[which(d_order$pval < Inf) & 
                                                                                           d_order$Location == "Top"]))
   yaxismax2 <- ifelse(freey == FALSE, max(d_order$pval[which(d_order$pval < 
                                                                Inf)]), max(d_order$pval[which(d_order$pval < Inf) & 
                                                                                           d_order$Location == "Bottom"]))
   yaxismin1 <- ifelse(freey == FALSE, 0, min(d_order$pval[d_order$Location == 
                                                             "Top"]))
   yaxismin2 <- ifelse(freey == FALSE, 0, min(d_order$pval[d_order$Location == 
                                                             "Bottom"]))
   backpanel1 <- ifelse(background == "white", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)")
   backpanel2 <- ifelse(background == "white", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin2, ymax = Inf, fill=factor(shademap)), alpha = 0.5)")
   p1 <- ggplot() + eval(parse(text = backpanel1))
   if ("Shape" %in% topn) {
     p1 <- p1 + geom_point(data = d_order[d_order$Location == 
                                            "Top", ], aes(x = pos_index, y = pval, color = factor(Color), 
                                                          shape = factor(Shape)), alpha = opacity)
   }
   else {
     p1 <- p1 + geom_point(data = d_order[d_order$Location == 
                                            "Top", ], aes(x = pos_index, y = pval, color = factor(Color)), 
                           alpha = opacity)
   }
   p1 <- p1 + scale_x_continuous(breaks = lims$av, labels = lims$Color, 
                                 expand = c(0, 0))
   if (chrblocks == TRUE) {
     p1 <- p1 + geom_rect(data = lims, aes(xmin = posmin - 
                                             0.5, xmax = posmax + 0.5, ymin = -Inf, ymax = min(d_order$pval), 
                                           fill = as.factor(Color)), alpha = 1)
   }
   p1 <- p1 + scale_colour_manual(name = "Color", values = newcols) + 
     scale_fill_manual(name = "Color", values = newcols)
   p1 <- p1 + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
                    axis.title.x = element_blank(), legend.position = "top", 
                    legend.title = element_blank())
   p2 <- ggplot() + eval(parse(text = backpanel2))
   if ("Shape" %in% bottomn) {
     p2 <- p2 + geom_point(data = d_order[d_order$Location == 
                                            "Bottom", ], aes(x = pos_index, y = pval, color = factor(Color), 
                                                             shape = factor(Shape)), alpha = opacity)
   }
   else {
     p2 <- p2 + geom_point(data = d_order[d_order$Location == 
                                            "Bottom", ], aes(x = pos_index, y = pval, color = factor(Color)), 
                           alpha = opacity)
   }
   p2 <- p2 + scale_x_continuous(breaks = lims$av, labels = lims$Color, 
                                 expand = c(0, 0))
   if (chrblocks == TRUE) {
     p2 <- p2 + geom_rect(data = lims, aes(xmin = posmin - 
                                             0.5, xmax = posmax + 0.5, ymin = -Inf, ymax = min(d_order$pval), 
                                           fill = as.factor(Color)), alpha = 1)
   }
   p2 <- p2 + scale_colour_manual(name = "Color", values = newcols) + 
     scale_fill_manual(name = "Color", values = newcols)
   p2 <- p2 + theme(axis.text.x = element_text(angle = 90), 
                    panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
                    axis.title.x = element_blank(), legend.position = "bottom", 
                    legend.title = element_blank())
   if (!missing(highlight_snp)) {
     if ("Shape" %in% topn) {
       p1 <- p1 + geom_point(data = d_order[d_order$SNP %in% 
                                              highlight_snp & d_order$Location == "Top", ], 
                             aes(x = pos_index, y = pval, shape = Shape), 
                             colour = highlighter)
       p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
     }
     else {
       p1 <- p1 + geom_point(data = d_order[d_order$SNP %in% 
                                              highlight_snp & d_order$Location == "Top", ], 
                             aes(x = pos_index, y = pval), colour = highlighter)
     }
     if ("Shape" %in% bottomn) {
       p2 <- p2 + geom_point(data = d_order[d_order$SNP %in% 
                                              highlight_snp & d_order$Location == "Bottom", 
                                            ], aes(x = pos_index, y = pval, shape = Shape), 
                             colour = highlighter)
       p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
     }
     else {
       p2 <- p2 + geom_point(data = d_order[d_order$SNP %in% 
                                              highlight_snp & d_order$Location == "Bottom", 
                                            ], aes(x = pos_index, y = pval), colour = highlighter)
     }
   }
   if (!missing(highlight_p)) {
     if ("Shape" %in% topn) {
       p1 <- p1 + geom_point(data = d_order[d_order$pvalue < 
                                              highlight_p[1] & d_order$Location == "Top", ], 
                             aes(x = pos_index, y = pval, shape = Shape), 
                             colour = highlighter)
       p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
     }
     else {
       p1 <- p1 + geom_point(data = d_order[d_order$pvalue < 
                                              highlight_p[1] & d_order$Location == "Top", ], 
                             aes(x = pos_index, y = pval), colour = highlighter)
     }
     if ("Shape" %in% bottomn) {
       p2 <- p2 + geom_point(data = d_order[d_order$pvalue < 
                                              highlight_p[2] & d_order$Location == "Bottom", 
                                            ], aes(x = pos_index, y = pval, shape = Shape), 
                             colour = highlighter)
       p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
     }
     else {
       p2 <- p2 + geom_point(data = d_order[d_order$pvalue < 
                                              highlight_p[2] & d_order$Location == "Bottom", 
                                            ], aes(x = pos_index, y = pval), colour = highlighter)
     }
   }
   if (!missing(tline)) {
     for (i in 1:length(tline)) {
       p1 <- p1 + geom_hline(yintercept = tredline[i], colour = "red", linetype='dashed')
     }
   }
   if (!missing(bline)) {
     for (i in 1:length(bline)) {
       p2 <- p2 + geom_hline(yintercept = bredline[i], colour = "red", linetype='dashed')
     }
   }
   if (!missing(annotate_p)) {
     if (!requireNamespace(c("ggrepel"), quietly = TRUE) == 
         TRUE) {
       print("Consider installing 'ggrepel' for improved text annotation")
       p1 <- p1 + geom_text(data = d_order[d_order$pvalue < 
                                             annotate_p[1] & d_order$Location == "Top", ], 
                            aes(pos_index, pval, label = SNP))
       p2 <- p2 + geom_text(data = d_order[d_order$pvalue < 
                                             annotate_p[2] & d_order$Location == "Bottom", 
                                           ], aes(pos_index, pval, label = SNP))
     }
     else {
       p1 <- p1 + ggrepel::geom_text_repel(data = d_order[d_order$pvalue < 
                                                            annotate_p[1] & d_order$Location == "Top", ], 
                                           aes(pos_index, pval, label = SNP))
       p2 <- p2 + ggrepel::geom_text_repel(data = d_order[d_order$pvalue < 
                                                            annotate_p[2] & d_order$Location == "Bottom", 
                                                          ], aes(pos_index, pval, label = SNP))
     }
   }
   if (!missing(annotate_snp)) {
     if (!requireNamespace(c("ggrepel"), quietly = TRUE) == 
         TRUE) {
       print("Consider installing 'ggrepel' for improved text annotation")
       p1 <- p1 + geom_text(data = d_order[d_order$SNP %in% 
                                             annotate_snp & d_order$Location == "Top", ], 
                            aes(pos_index, pval, label = SNP))
       p2 <- p2 + geom_text(data = d_order[d_order$SNP %in% 
                                             annotate_snp & d_order$Location == "Bottom", 
                                           ], aes(pos_index, pval, label = SNP))
     }
     else {
       p1 <- p1 + ggrepel::geom_text_repel(data = d_order[d_order$SNP %in% 
                                                            annotate_snp & d_order$Location == "Top", ], 
                                           aes(pos_index, pval, label = SNP))
       p2 <- p2

     }
   }
   p1 <- p1 + ylab(yaxislab1)
   p2 <- p2 + ylab(yaxislab2)
   if (chrblocks == TRUE) {
     if (freey == TRUE) {
       print("Sorry, drawing chrblocks with freey=TRUE is currently unsupported and will be ignored.")
     }
     else {
       p1 <- p1 + theme(axis.text.x = element_text(vjust = 1), 
                        axis.ticks.x = element_blank()) + ylim(c(yaxismin1, 
                                                                 yaxismax1))
       p2 <- p2 + scale_y_reverse(limits = c(yaxismax2, 
                                             yaxismin2)) + theme(axis.text.x = element_blank(), 
                                                                 axis.ticks.x = element_blank())
     }
   }
   else {
     p1 <- p1 + theme(axis.text.x = element_text(vjust = 1), 
                      axis.ticks.x = element_blank()) + scale_y_continuous(limits = c(yaxismin1, 
                                                                                      yaxismax1), expand = expansion(mult = c(0, 0.1)))
     p2 <- p2 + scale_y_reverse(limits = c(yaxismax2, yaxismin2), 
                                expand = expansion(mult = c(0.1, 0))) + theme(axis.text.x = element_blank(), 
                                                                              axis.ticks.x = element_blank())
   }
   if (background == "white") {
     p1 <- p1 + theme(panel.background = element_rect(fill = "white"))
     p2 <- p2 + theme(panel.background = element_rect(fill = "white"))
   }
   p1 <- p1 + guides(fill = "none", color = "none")
   p2 <- p2 + guides(fill = "none", color = "none")

   return(p1/p2)
}

#### Plot functions ####
plot_pca1 <- function(df) {
  p <- df %>%
    mutate(Clusters=as.factor(cluster)) %>% 
    # filter(!Clusters %in% c(5, 7, 6, 4, 3)) %>%
    ggplot(aes(x=PC1, y=PC2, col=Clusters, fill=Clusters))+
    geom_point()+
    # geom_point(pch=21, col='black', alpha=0.7)+
    scale_alpha(guide = 'none')+
    guides(color = guide_legend(override.aes = list(alpha = 1)))+
    theme_classic()+
    scale_color_npg()+
    scale_fill_npg()+
    guides(fill = guide_legend(override.aes = list(size = 5, alpha=1)))
  show(p)
}

plot_pca2 <- function(df) {
  p <- df %>%
    mutate(alpha=case_when(SuperPopulation=='not in 1kg' ~ 0,
                           TRUE~1),
           SuperPopulation=fct_relevel(SuperPopulation, 
                                       'EUR', 'AFR', 'AMR', 'EAS', 'SAS', 'not in 1kg'),
           shape=factor(case_when(SuperPopulation=='not in 1kg' ~ 'not in 1kg',
                                  TRUE~'in 1kg'))) %>% 
    arrange(desc(alpha)) %>% filter(SuperPopulation!='not in 1kg') %>% 
    ggplot(aes(x=PC1, y=PC2, col=SuperPopulation, alpha=alpha))+
    geom_point()+
    scale_color_npg(name='1kg samples')+
    theme_classic()+
    scale_alpha(guide = 'none')+
    guides(col = guide_legend(override.aes = list(size = 5, alpha=1)),
           shape = guide_legend(override.aes = list(size = 5, alpha=1)))
  show(p)
}

plot_pca3 <- function(df) {
  p <- df %>%  mutate(alpha=case_when(fsgs=='case' ~ 1,
                                      fsgs=='control' ~ 0.9)) %>% 
    ggplot(aes(x=PC1, y=PC2, col=fsgs, alpha=alpha))+
    geom_point()+
    scale_alpha(guide = 'none')+
    guides(color = guide_legend(override.aes = list(size=5)))+
    theme_classic()+
    scale_color_npg(name='')+
    scale_fill_npg()
  show(p)
}

plot_bar1 <- function(df) {
  df %>% group_by(fsgs) %>% summarise(Counts=n()) %>%
    mutate(lbl=case_when(fsgs=='control' ~ paste(Counts, "exomes"),
                         fsgs=='case' ~ paste(Counts, "target panels"))) %>% 
    ggplot(aes(x=fsgs, y=Counts, fill=fsgs))+
    geom_bar(stat='identity', position='dodge', col='black', show.legend = F)+
    geom_text(aes(label=lbl, x=fsgs, y=Counts), vjust = -0.6)+
    theme_classic()+
    scale_x_discrete('')+
    scale_color_npg(name='')+
    scale_fill_npg(name='FSGS status')
}

plot_bar2 <- function(df, ord){
  dd <- df %>% 
    mutate(cluster=factor(cluster, rev(ord))) %>% 
    group_by(fsgs, cluster) %>%
    summarise(Occurence = n()) %>% arrange(desc(Occurence), cluster, fsgs)
  
  p <- dd %>% ungroup() %>% 
    ggplot(aes(x=cluster, y=Occurence, fill=fsgs))+
    geom_bar(stat='identity', position='dodge', width=.6, col='black')+
    geom_text(aes(label=Occurence, x=cluster, y=Occurence), size=3,
              position = position_dodge(width=0.5), hjust = -0.1, vjust=0.6)+
    theme_classic()+
    scale_x_discrete('Cluster')+
    scale_y_continuous('Occurence', limits = c(0, 5000))+
    scale_fill_npg(name='')+
    coord_flip()
  show(p)
}

plot_match <- function(pca_before, pca_after, pca_subset=c('PC1', 'PC2', 'PC3'), title='') {
  p <- rbind(pca_before, pca_after) %>%
    dplyr::select(pca_subset, fsgs, matching) %>% gather(PCX, PC_value, -fsgs, -matching) %>% 
    mutate(matching=case_when(matching=='before' ~ 'Before matching',
                              matching=='after' ~ 'After matching'),
           matching=factor(matching, levels=c('Before matching', 'After matching'))) %>% 
    ggplot(aes(x=PC_value, fill=fsgs))+
    geom_density(alpha=0.8)+
    ggtitle(title,)+
    facet_rep_grid(matching~PCX,
                   switch='x', scales='free', repeat.tick.labels = T)+
    scale_fill_npg(name='', labels=c('Case', 'Control'))+
    scale_x_continuous('', guide = guide_axis(check.overlap = TRUE))+
    scale_y_continuous('Density')+
    theme_classic()+
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
  show(p)
} 
