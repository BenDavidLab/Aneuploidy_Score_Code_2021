#Recreating the estimated ploidy + Aneuploidy score for the cell line using Jame's code:
library(tidyverse)
library(matrixStats)
library(modi)
library(plyr)
library(dplyr)
library(matrixStats)
library(BiocManager)
library(ggplot2)

### spatstat.geom v2.3-0
# "http://cran.r-project.org/src/contrib/Archive/spatstat.geom/spatstat.geom_2.3-0.tar.gz"
install.packages("http://cran.r-project.org/src/contrib/Archive/spatstat.geom/spatstat.geom_2.3-0.tar.gz", repos=NULL, type="source")
library(spatstat.geom)

# Code: -------------------------------------------------------------------
# James used public-21q2-110d, which is DepMap Public 21Q2
CCLE.segment.cn <- read.csv('./CCLE_segment_cn_21Q2.csv')

#Load ploidy table

est_ploidy_simple <- function(CN_seg_df, fig_dir = NULL, return_fig = FALSE, CCLE_ABS_TABLE = NULL, min_ploidy = 1.25, max_ploidy = 10, n_ploidy_gridvals = 500,
                              diploid_prior_weight = 0.01) {
  poss_sc_facs <- seq(from = min_ploidy, to = max_ploidy, length.out = n_ploidy_gridvals)
  
  #create simple prior to slightly bias towards diploid 
  tent_fn <- function(x, cp = 2, slope = diploid_prior_weight, min_prob = 0.1) {
    y <- 1+(x-cp)*slope
    gt <- x > cp
    y[gt] <- 1-(x[gt]-cp)*slope
    y[y < min_prob] <- min_prob
    return(y)
  }
  prior <- tent_fn(poss_sc_facs)
  
  cls <- unique(CN_seg_df$DepMap_ID)
  simp_ploidy_res <- ldply(cls, function(targ) {
    df <- CN_seg_df %>% filter(DepMap_ID == targ) %>% mutate(width = End - Start)
    
    scaled_vals <- df$Segment_Mean %o% poss_sc_facs
    err <- colSums(df$width * (scaled_vals - round(scaled_vals))^2)/sum(df$width)
    res <- tibble(sc = poss_sc_facs, err = err, p = prior) %>% 
      mutate(post = p*(1-err))
    
    #find scaling factor that maximizes the 'posterior'
    best_sc <- res %>% arrange(dplyr::desc(post)) %>% head(1)
    
    #compute the residual var of the integer CN fit 
    df = df %>% mutate(resid = Segment_Mean*best_sc$sc - round(Segment_Mean*best_sc$sc)) #calculate residual from integer model
    rvar <- weightedVar(df$resid, df$width)
    
    #make fig if neeeded
    if (!is.null(fig_dir) | return_fig) {
      print(targ)
      abs_ploidy <- CCLE_ABS_TABLE %>% filter(DepMap_ID == targ) %>% pull(ploidy)
      abs_sc <- CCLE_ABS_TABLE %>% filter(DepMap_ID == targ) %>% pull(`Subclonal genome fraction`)
      
      orig_thresh <- 0.5*best_sc$sc
      new_thresh <- round(best_sc$sc)/2
      
      g1 <- ggplot(res, aes(sc, 1-err)) +
        geom_point() +
        geom_point(aes(sc, post), color = 'blue') +
        geom_point(data = best_sc, color = 'red', size = 2) +
        xlim(0,10)
      g2 <- ggplot(df, aes(x = Segment_Mean*best_sc$sc, y = ..density.., weight = width)) +
        geom_histogram(bins = 200) +
        geom_vline(xintercept = seq(0, 10, by = 1), linetype = 'dashed') +
        xlim(0, 10) +
        geom_vline(xintercept = orig_thresh, color = 'red') +
        geom_vline(xintercept = new_thresh, color = 'blue') +
        ggtitle(sprintf('Inferred ploidy, %s, %.2f, resid-var: %.4f', unique(df$Source), best_sc$sc, rvar))
      g3 <- ggplot(df, aes(x = Segment_Mean*abs_ploidy, y = ..density.., weight = width)) +
        geom_histogram(bins = 200) +
        geom_vline(xintercept = seq(0, 10, by = 1), linetype = 'dashed') +
        xlim(0, 10) +
        ggtitle(sprintf('ABS ploidy, Sc frac: %.2f, ploidy: %.2f', abs_sc, abs_ploidy))
      p <- cowplot::plot_grid(g1, g2, g3, ncol = 1)
      
      if (return_fig) {
        print(p)
      }
      if (!is.null(fig_dir)) {
        ggsave(file.path(fig_dir, paste0(targ, '_simple_ploidy.png')), plot = p, height = 9, width = 5)
      }
    }
    resTable = tibble(DepMap_ID = targ,
                      min_err = best_sc$err,
                      est_ploidy = best_sc$sc,
                      resid_var = rvar,
                      Source = unique(df$Source))
  })
}

# James used diploid_prior_weight = 0.0075 to create simp_ploidy_res 
simp_ploidy_res <- est_ploidy_simple(CCLE.segment.cn, diploid_prior_weight = 0.0075)


#Aneuploidy calculation (as in Cohen & Sarir)
#*** This code is what was used to estimate aneuploidy scores in Cohen-Sharir et al.

#assigns chromosome arm to each segment using centromere range data
get_which_arm <- function(seg_start, seg_end, cent_start, cent_end) {
  output <- rep(NA, length(seg_start))
  seg_cent = 0.5*(seg_start + seg_end)
  output[seg_cent < cent_start] <- 'p'
  output[seg_cent > cent_end] <- 'q'
  return(output)
}


#call arm events based on comparing arm-level CN to ploidy values
arm_call <- function(CN_vals, ploidy_vals){
  arm_CN <- round(CN_vals)
  ploidy_CN <- round(ploidy_vals)
  arm_calls <- rep(0, times = length(CN_vals))
  arm_calls[arm_CN < ploidy_CN] <- -1
  arm_calls[arm_CN > ploidy_CN] <- 1
  arm_calls
}



#handle cases where a segment crosses over the centromere by splitting it there
split_cent_crosses <- function(seg_df) {
  cross_segs <- seg_df$Start < seg_df$centStart & seg_df$End > seg_df$centEnd
  left_side <- seg_df[cross_segs, ] %>% 
    mutate(End = 0.5*(centStart + centEnd))
  right_side <- seg_df[cross_segs, ] %>% 
    mutate(Start = 0.5*(centStart + centEnd))
  seg_df[cross_segs, ] <- left_side
  seg_df <- rbind(seg_df, right_side)
  return(seg_df)
}



# Load data ---------------------------------------------------------------
#load hg38 centromere info
cent_df <- rCGH::hg38 %>% 
  dplyr::mutate(chrom = as.character(chrom)) %>% 
  dplyr::select(chrom,
                centStart = centromerStart,
                centEnd = centromerEnd)


seg_df <- CCLE.segment.cn %>%    
  left_join(simp_ploidy_res, by = 'DepMap_ID') %>%
  mutate(CN = round(Segment_Mean * est_ploidy)) %>% 
  filter(!(Chromosome %in% c('X', 'Y'))) %>% 
  dplyr::select(chrom = Chromosome,
                Start,
                End,
                CN,
                DepMap_ID) %>% 
  dplyr::mutate(chrom = as.character(chrom)) %>% 
  left_join(cent_df, by = 'chrom') %>% 
  split_cent_crosses() %>% #split segments at centromeres
  dplyr::mutate(width = End-Start) #compute segment widths



# Calculate arm-calls -----------------------------------------------------
#assign arms to each segment and compute total arm lengths and relative weights
# !!! weighted.median func was updated in spatstat.geom v2.3-1. make sure to run it with spatstat.geom v2.3-0 to get the same results.
arm_calls <- seg_df %>% 
  dplyr::mutate(arm = get_which_arm(Start, End, centStart, centEnd)) %>% 
  dplyr::filter(!is.na(arm)) %>% 
  dplyr::group_by(DepMap_ID, chrom, arm) %>% 
  dplyr::summarize(wmed_CN = weighted.median(CN, width, na.rm=T)) %>%
  ungroup() %>% 
  left_join(
    simp_ploidy_res %>% dplyr::select(DepMap_ID, ploidy = est_ploidy),
    by = 'DepMap_ID'
  ) %>% dplyr::mutate(arm_call = arm_call(wmed_CN, ploidy)) 


#save matrix of arm-calls
arm_ord <- c('1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q', '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q', '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q', '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q')
arm_calls <- arm_calls %>% 
  dplyr::mutate(chrom_arm = paste0(chrom, arm),
                chrom_arm = factor(chrom_arm, levels = arm_ord)) %>% 
  dplyr::select(DepMap_ID, chrom_arm, arm_call) %>% 
  filter(!is.na(chrom_arm)) 


# Compute aneuploidy features ---------------------------------------------
arm_CNVs_df <- arm_calls %>% 
  dplyr::mutate(arm_event = abs(arm_call)) %>% 
  dplyr::group_by(DepMap_ID) %>% 
  dplyr::summarise(num_arm_events = sum(arm_event, na.rm=T)) %>% 
  ungroup() 

df <- merge(arm_CNVs_df, simp_ploidy_res, by = c("DepMap_ID"))
