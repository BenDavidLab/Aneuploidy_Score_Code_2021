# Simple ploidy estimate
# James McFarland
# 7/9/2021
#Load data:
CCLE.segment.cn <- load.from.taiga(data.name='public-21q2-110d', data.version=13, data.file='CCLE_segment_cn')
## Fetching https://cds.team/taiga/api/dataset/public-21q2-110d/13 
## Status 200 
## loading cached data version from  /Users/jmmcfarl/.taiga/public-21q2-110d_13.toc
CCLE.ABSOLUTE.combined.table <- load.from.taiga(data.name='ccle-absolute-cn', data.version=5, data.file='CCLE_ABSOLUTE_combined_table')
## Fetching https://cds.team/taiga/api/dataset/ccle-absolute-cn/5 
## Status 200 
## loading cached data version from  /Users/jmmcfarl/.taiga/ccle-absolute-cn_5.toc

#Defining the function
#This function takes as input the segmented copy number data (with relative CN estimates) and returns a table with ploidy estimates for each line (and some other stats).

# It���s based on a very simple approach of trying to find a scaling factor for the relative copy number data (i.e.��ploidy) that produces a CN distribution that is most centered around integer values. It does this by minimizing a squared error loss of the integer-CN model, with a slight (heuristically applied) bias to encourage near-diploid ploidy estimates in more ambiguous cases. Importantly, it assumes that there is perfect cancer cell purity, and low subclonal fraction (i.e.��assumes the true distribution should be peaked at integer values after scaling by ploidy).

#You can use these est_ploidy values to convert relative CN to absolute in the segmented profile like this:
  
  #CN_seg_df %<>% mutate(abs_CN = Segment_Mean*est_ploidy)

# and for the gene level CN (which is in log2(1+x) units) you can do this:
  #gene_abs_CN = (2^gene_rel_CN-1) * est_ploidy

#' Title
#'
#' @param CN_seg_df: segmented CN file formatted as in DepMap release files 
#' @param fig_dir: path to directory to print out figs
#' @param return_fig: Bool to indicate whether to print out fig to console (for testing) 
#' @param CCLE_ABS_TABLE: table of ABSOLUTE stats (used for making diagnostic plots only) 
#' @param min_ploidy: minimum ploidy values tested 
#' @param max_ploidy: max ploidy value tested 
#' @param n_ploidy_gridvals number of grid points used to test ploidy vals
#' @param diploid_prior_weight strength of prior biasing towards diploidy ploidy
#'
#' @return
#' @export
#'
#' @examples
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
    df %<>% mutate(resid = Segment_Mean*best_sc$sc - round(Segment_Mean*best_sc$sc)) #calculate residual from integer model
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
    
    tibble(DepMap_ID = targ,
           min_err = best_sc$err,
           est_ploidy = best_sc$sc,
           resid_var = rvar,
           Source = unique(df$Source))
  })
}

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
library(spatstat)
#assign arms to each segment and compute total arm lengths and relative weights
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
Compare with original aneuploidy scores (using ABSOLUTE data)
old_aneuploidy_data <- read_csv('~/CPDS/data/CCLE/aneuploidy_data.csv')

comb_aneuploidy <- full_join(arm_CNVs_df, old_aneuploidy_data %>% dplyr::select(DepMap_ID, old_arm_events = num_arm_events), by = 'DepMap_ID')
ggplot(comb_aneuploidy, aes(old_arm_events, num_arm_events)) + 
  geom_point(pch = 21, fill = 'black', color = 'white', stroke = 0.1) + 
  geom_abline() +
  ggpubr::stat_cor()


comb_aneuploidy %>% 
  left_join(CCLE.ABSOLUTE.combined.table) %>% 
  ggplot(aes(old_arm_events, num_arm_events, fill = `Subclonal genome fraction` > 0.4)) + 
  geom_jitter(pch = 21, color = 'white', stroke = 0.1) + 
  geom_abline() +
  ggpubr::stat_cor()


Create high and low-aneuploidy groups
(As in Cohen-Sharir)

#return boolean vector which is true if data is greater than upper threshold, and FALSE if lower than lower threshold
#use as_quantiles = T to specify whether the thresholds should be interpreted as quantiles
two_threshold_split <- function(data, thresh_vals = c(0.25, 0.75), as_quantiles = T) {
  group <- rep(NA, length(data))
  if (as_quantiles) {
    thresh_vals = quantile(data, thresh_vals, na.rm=T)
  }
  group[data >= thresh_vals[2]] <- TRUE
  group[data <= thresh_vals[1]] <- FALSE
  group
}
comb_aneuploidy %<>% 
  mutate(many_arm_events_old = two_threshold_split(old_arm_events),
         many_arm_events = two_threshold_split(num_arm_events))
