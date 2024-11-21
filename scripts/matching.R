# Matching
source('scripts/functions.R')

library(dplyr)
library(ggplot2)
library(data.table)
library(ggsci)
library(lemon)
library(tidyverse)
library(MatchIt)

set.seed(42)

# Download global population PCA-dataset with predefined AutoGMM-clustering (Thomas L. Athey, 2020).
df <- as_data_frame(fread('data/all_pca.annotated.tsv'))
df <- as_tibble(fread('data/plink.all.eigenvec')) %>% 
  select(-V2) %>% select(V1, V3:V12)
colnames(df) <- c("s","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
fsgs <- as_tibble(fread('data/fsgs.phe')) %>% rename(s=ID) %>% select(s, FSGS)
df <- left_join(df, fsgs)
df <- df %>% mutate(fsgs=case_when(FSGS==2 ~ 'case',
                                   FSGS==1 ~ 'control'))
######################
## Data preparation ##
######################

# This section provides a basic understanding of case-controls distribution after AutoGMM-clustering.

# Global PCA-scatter
plot_pca1(df)

# Obtain 1000 genome annotation
kg1 <- as_data_frame(fread('data/1kg_annotations.txt')) %>% dplyr::select(Sample, SuperPopulation, Population) %>%
  rename(s=Sample)
df <- left_join(df, kg1, by='s')
df <- df %>% replace_na(list(SuperPopulation='not in 1kg'))

# Presence of 1kg-population per cluster
df %>% group_by(cluster, SuperPopulation) %>% summarise(n=n())
plot_pca2(df)

# Presence of case/control per cluster
df %>% group_by(fsgs, cluster) %>% summarise(n = n()) %>% arrange(desc(n), cluster, fsgs)
plot_pca3(df)
plot_bar1(df)
plot_bar2(df, c(3, 5, 0, 2, 4, 6, 7, 1))

# Prepare data for further analysis
df <- df %>% mutate(FSGS=case_when(fsgs=='control' ~ 0,
                                   fsgs=='case' ~ 1),
                    matching='before')

# Split PCA-dataset into population specific set of clusters
## major clusters
eur <- df %>% filter(cluster %in% c(0, 3, 5))  # European
afr <- df %>% filter(cluster==2)  # African

## minor clusters
eas <- df %>% filter(cluster==7)  # East Asian
sas <- df %>% filter(cluster==1)  # South Asian
amr <- df %>% filter(cluster %in% c(4, 6))  # Ad Mixed American

############################
## Whole dataset matching ##
############################

# Here we are reproducing authors global meta matching (based on weighted euclidean metrics) (Minxian Wang, 2019)
# with slightly different approach in same case-controls proportions (1:1). For applying exactly same approach
# can be used another settings of matchit function (see "?matchit").

meta_result <- matchit(FSGS ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,
                       method = 'genetic', discard = 'both', reestimate = T, replace = T, ratio = 2)

# Visual estimation of the amount of filtered data
plot(meta_result, type = "jitter", interactive = FALSE)
# A visual representation of the overall assessment of the quality of the match.
plot(meta_result, type = 'histogram')

# Determine which PC contribute more to data variance.
plot(summary(meta_result))
# Let's focus on the three most significant PC.
plot(meta_result, type = 'qq', which.xs = c('PC2', 'PC3', 'PC5'))
plot(meta_result, type = 'ecdf', which.xs = c('PC2', 'PC3', 'PC5'))

# Obtain matched dataset.
meta_matched <- match.data(meta_result)
meta_matched$matching <- 'after'
# Quantifying how much data is left
meta_matched %>% group_by(fsgs) %>% summarise(n=n())

# Let's take a look at the difference in the distribution of cases and controls in the most significant PC
plot_match(df, dplyr::select(meta_matched, -distance, -weights), c('PC2', 'PC3', 'PC5'), 'Whole dataset matching')


rbind(df %>% filter(PC2>-0.02, PC3<0.25), dplyr::select(meta_matched, -distance, -weights)) %>%
  dplyr::select(PC2,PC3, fsgs, matching) %>%
  mutate(matching=case_when(matching=='before' ~ 'Before matching',
                            matching=='after' ~ 'After matching'),
         matching=factor(matching, levels=c('Before matching', 'After matching'))) %>% 
  ggplot(aes(x=PC2, PC3, col=fsgs))+
  geom_point()+
  facet_wrap(~matching,
                 switch='x')+
  scale_fill_npg(name='', labels=c('Case', 'Control'))+
  scale_x_continuous('PС2', guide = guide_axis(check.overlap = TRUE))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        plot.title = element_text(hjust = 0.5))

# Write the matched samples
write(meta_matched$s , 'data/new.x2.matched_samples.txt')

##################################
## Population specific matching ##
##################################

# In this section, we propose a local population-based matching approach. The intermediate quality of case-control
# matching is estimated based on built in metrics of Matchit package (Ho DE, 2011). The resulting quality of 
# case-control matching is assessed by the QQ-plot on synonymous rare variants (see assoc.R). All functions used
# below are the same as in the section above. Only the composition of the selected clusters changes.

## European matching ##
eur_result <- matchit(FSGS ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = eur,
                      method = 'genetic', discard = 'both', reestimate = T, replace = T, ratio = 5)

# Visual estimation of the amount of filtered data
plot(eur_result, type = "jitter", interactive = FALSE)
# A visual representation of the overall assessment of the quality of the match.
plot(eur_result, type = 'histogram')

# Determine which PC contribute more to data variance.
plot(summary(eur_result))
# Let's focus on the three most significant PC.
plot(eur_result, type = 'qq', which.xs = c('PC1', 'PC3', 'PC5'))
plot(eur_result, type = 'ecdf', which.xs = c('PC1', 'PC3', 'PC5'))

# Obtain matched dataset.
eur_matched <- match.data(eur_result)
eur_matched$matching <- 'after'

# Let's take a look at the difference in the distribution of cases and controls in the most significant PC
plot_match(eur, dplyr::select(eur_matched, -distance, -weights), c('PC1', 'PC3', 'PC5'), 'European matching')

# Write the matched samples
write(eur_matched$s , 'data/eur.matched_samples.txt')

## African matching ##
afr_result <- matchit(FSGS ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = afr,
                       method = 'genetic', discard = 'both', reestimate = T, replace = T, ratio = 5)

plot(afr_result, type = "jitter", interactive = FALSE)
plot(afr_result, type = "histogram")
plot(summary(afr_result))
plot(afr_result, type = 'qq', which.xs = c('PC1', 'PC8', 'PC9'))
plot(afr_result, type = 'ecdf', which.xs = c('PC1', 'PC8', 'PC9'))
afr_matched <- match.data(afr_result)
afr_matched$matching <- 'after'
plot_match(afr, dplyr::select(afr_matched, -distance, -weights), c('PC1', 'PC8', 'PC9'), 'African matching')
write(afr_matched$s , 'data/afr.matched_samples.txt')

## Ad Mixed American matching ##
amr_result <- matchit(FSGS ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = amr,
                      method = "genetic", replace = TRUE, mahvars = ~ PC1+PC2+PC3+PC4+PC6+PC7,
                      ratio = 3,
                      pop.size = 10)

plot(amr_result, type = "jitter", interactive = FALSE)
plot(amr_result, type = "histogram")
plot(summary(amr_result))
plot(amr_result, type = 'qq', which.xs = c('PC1', 'PC5', 'PC6'))
plot(amr_result, type = 'ecdf', which.xs = c('PC1', 'PC5', 'PC6'))
amr_matched <- match.data(amr_result)
amr_matched$matching <- 'after'
plot_match(amr, dplyr::select(amr_matched, -distance, -weights), c('PC1', 'PC5', 'PC6'), 'Ad Mixed American matching')
write(amr_matched$s , 'data/amr.matched_samples.txt')
