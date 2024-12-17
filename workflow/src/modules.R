

#' Prepare data in the right format to make a manhattan plot
#' @param dt A dataframe with the columns: `id`, `pvalue`, `chrom`, `locus`
#' @return a list
prepare_manhattan_dt <- function(dt){
  # https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
  dt <- dt %>%
    dplyr::mutate(chrom = as.numeric(chrom))
  data_cum <- dt %>% 
    dplyr::mutate(chrom=as.numeric(chrom)) %>%
    dplyr::group_by(chrom) %>% 
    dplyr::summarise(max_bp = max(locus)) %>% 
    dplyr::mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
    dplyr::select(chrom, bp_add)

  new_dt <- dt %>% 
    dplyr::inner_join(data_cum, by = "chrom") %>% 
    dplyr::mutate(bp_cum = locus + bp_add) %>%
    dplyr::arrange(chrom, bp_cum)

  axis_set <- new_dt %>% 
    dplyr::group_by(chrom) %>% 
    dplyr::summarize(center = mean(bp_cum)) %>%
    dplyr::arrange(chrom)

  ylim <- new_dt %>% 
    dplyr::filter(pvalue == min(pvalue)) %>% 
    dplyr::mutate(ylim = abs(floor(log10(pvalue))) + 2) %>% 
    dplyr::pull(ylim)

  return(list(x=new_dt, axis_set=axis_set, ylim=ylim))
}

plot_manhattan_dt <- function(manhattan_list, signif=0.05/nrow(manhattan_list$x), plot_id=TRUE){
  # # https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

  tfwas_data <- manhattan_list$x
  axis_set <- manhattan_list$axis_set
  ylim <- manhattan_list$ylim

  ppt <- ggplot(tfwas_data, aes(x = bp_cum, y = -log10(pvalue), color = as_factor(chrom), size = -log10(pvalue))) +
    geom_hline(yintercept = -log10(signif), color = "grey40", linetype = "dashed") + 
    geom_point(alpha = 0.75, size=0.3) +
    scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#000000", "#999999"), unique(length(axis_set$chrom)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "-log(p)") + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )
  
  if(plot_id == TRUE){
    ppt + geom_text_repel(aes(label = ifelse(-log10(pvalue) > -log10(signif), as.character(id),'')), 
    box.padding = unit(0.2, "lines"), point.padding = 0.5, force = 100, segment.size  = 0.1, segment.color = "grey50")
   } else {
    ppt 
    }
}


## pvalue vs uniform

qqunif = 
  function(p,BH=T,CI=T,mlog10_p_thres=30,...)
  {
    ## thresholded by default at 1e-30
    p=na.omit(p)
    nn = length(p)
    xx =  -log10((1:nn)/(nn+1))
    
    p_thres = 10^{-mlog10_p_thres}
    if( sum( p < p_thres) )
    {
      warning(paste("thresholding p to ",p_thres) )
      p = pmax(p, p_thres)
    }
    plot( xx,  -sort(log10(p)),
          xlab=expression(Expected~~-log[10](italic(p))),
          ylab=expression(Observed~~-log[10](italic(p))),
          cex.lab=1.4,mgp=c(2,1,0),
          ... )
    abline(0,1,col='gray')
    if(BH)
    {
      abline(-log10(0.05),1, col='orange',lty=1)
      abline(-log10(0.10),1, col='orange',lty=2)
      abline(-log10(0.25),1, col='orange',lty=3)
      legend('topleft', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
             col=c('orange','orange','orange'),lty=1:3, cex=1)
      abline(h=-log10(0.05/nn)) ## bonferroni
    }
    if(CI)
    {
      ## create the confidence intervals
      c95 <- rep(0,nn)
      c05 <- rep(0,nn)
      ## the jth order statistic from a
      ## uniform(0,1) sample
      ## has a beta(j,n-j+1) distribution
      ## (Casella & Berger, 2002,
      ## 2nd edition, pg 230, Duxbury)
      ## this portion was posted by anonymous on
      ## http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
      
      for(i in 1:nn)
      {
        c95[i] <- qbeta(0.95,i,nn-i+1)
        c05[i] <- qbeta(0.05,i,nn-i+1)
      }
      
      lines(xx,-log10(c95),col='gray')
      lines(xx,-log10(c05),col='gray')
    }
  }