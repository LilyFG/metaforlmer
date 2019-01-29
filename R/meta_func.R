#' Meta-analysis of multiple predictors in a list of models with the same structure.
#'
#' \code{meta_models} returns a list containing metafor outputs for each predictor and a table of the combined results.
#'
#' This code currently pulls the beta and standard error from the each predictor
#' in the model and uses the \code{\link{rma.uni()}} function from
#' \code{\link{metafor}} to run meta-analyses over each predictor.
#'
#' @param model_list A list of lmer or glmer models with the same structure.
#'
#' @return Returns a dataframe for quick reference and a list of the full
#'  meta-analysis outputs from \code{\link{metafor}}.


meta_models <- function(model_list){
  studies <- names(model_list)
  coefs <- data.frame(lapply(model_list, function(x){
    summary(x)$coef[1:nrow(summary(x)$coef),1:2, drop = F]}))
  coefs <- coefs[,grepl(paste(gsub(" ", ".", studies),collapse="|"), colnames(coefs))]
  metas <- lapply(as.list(row.names(coefs)), function(predictor){
    data <- data.frame(yi = as.numeric(coefs[predictor,seq(1,ncol(coefs)-1,2)]),
                       sei = as.numeric(coefs[predictor,seq(2,ncol(coefs),2)]))
    rma.uni(yi=yi, sei=sei, data = data, slab = studies, method = "REML")
  })
  names(metas) <- row.names(coefs)
  df <- data.frame(
                   beta = unlist(lapply(metas, function(x)x$beta)),
                   se = unlist(lapply(metas, function(x)x$se)),
                   ci.lb = unlist(lapply(metas, function(x)x$ci.lb)),
                   ci.ub = unlist(lapply(metas, function(x)x$ci.ub)),
                   z = unlist(lapply(metas, function(x)x$zval)),
                   p = unlist(lapply(metas, function(x)x$pval)),
                   psig = ifelse(unlist(lapply(metas, function(x)x$pval))<0.05,"*",""),
                   I2 = unlist(lapply(metas, function(x)x$I2)),
                   Q = unlist(lapply(metas, function(x)x$QE)),
                   Qp = unlist(lapply(metas, function(x)x$QEp)),
                   Qpsig = ifelse(unlist(lapply(metas, function(x)x$QEp))<0.05,"*","")
                   )
  list(metas=metas, df=df)
}

#
# meta_func_sem <- function(med_list, studies){
#   pe1 <- med_list[[1]]$pe
#   paths <- paste(med_list[[1]]$pe$lhs, med_list[[1]]$pe$op, med_list[[1]]$pe$rhs)
#   labels <- med_list[[1]]$pe$label
#   pe <- lapply(med_list, function(x)x$pe)
#   pe <- as.data.frame(data.table::rbindlist(pe))
#   pe$path <- paste(pe$lhs, pe$op, pe$rhs)
#   metas <- lapply(as.list(paths), function(path){
#     data <- data.frame(yi = as.numeric(pe[pe$path==path,"est"]),
#                        sei = as.numeric(pe[pe$path==path,"se"]))
#     suppressWarnings(rma.uni(yi=yi, sei=sei, data = data, slab = studies, method = "REML"))
#   })
#   names(metas) <- paths
#   df <- data.frame(
#     lhs = pe1$lhs,
#     op = pe1$op,
#     rhs = pe1$rhs,
#     label = labels,
#     est = unlist(lapply(metas, function(x)x$beta)),
#     se = unlist(lapply(metas, function(x)x$se)),
#     ci.lower = unlist(lapply(metas, function(x)x$ci.lb)),
#     ci.upper = unlist(lapply(metas, function(x)x$ci.ub)),
#     z = unlist(lapply(metas, function(x)x$zval)),
#     pvalue = unlist(lapply(metas, function(x)x$pval)),
#     psig = ifelse(unlist(lapply(metas, function(x)x$pval))<0.05,"*",""),
#     I2 = unlist(lapply(metas, function(x)x$I2)),
#     Q = unlist(lapply(metas, function(x)x$QE)),
#     Qp = unlist(lapply(metas, function(x)x$QEp)),
#     Qpsig = ifelse(unlist(lapply(metas, function(x)x$QEp))<0.05,"*","")
#   )
#   list(metas=metas, df=df)
# }
#


#' Meta-analysis of multiple predictors in a list of models with the same structure.
#'
#' \code{ggforest} returns a list containing metafor outputs for each predictor and a talbe of the combined results.
#'
#' This function takes the list returned by .
#'
#' @param x The list returned by the meta_models function.
#' @param intercept A Boolean. Do you want to plot the model intercept? Defaults to \code{FALSE}.
#' @param labels A vector same length as number of predictors to contain pretty names for the predictors.
#' @param hetero A character string with the desired heterogeneity measure to be displayed.
#'  Options: "none" {default}; "I2"; "Q"
#' @param palette A vector containing hex codes for the desired colour palette. Defaults to Darjeeling1 from
#'  \code{\link{wes_palettes}}.
#'
#' @return Returns a dataframe for quick reference and a list of the full
#'  meta-analysis outputs from \code{\link{metafor}}.

ggforest <- function(x, intercept=F, labels = NULL, hetero = "none", palette = wesanderson::wes_palettes$Darjeeling1){
  require("ggplot2")
  # Function to convert REM results in `rma`-format into a data.frame
  rma2df = function(x){
    df <- as.data.frame(data.table::rbindlist(lapply(x, function(x){
    rbind(
      data.frame(Study = "RE Model", LogFC = x$b, CILB=x$ci.lb, CIUB=x$ci.ub,
                 p = x$pval, group = "Meta-\nanalysis",
                 I2 = x$I2,
                 Qp = x$QEp,
                 stringsAsFactors = FALSE),
      data.frame(Study = x$slab, LogFC = x$yi,
                 CILB=x$yi - 2*sqrt(x$vi),
                 CIUB=x$yi + 2*sqrt(x$vi),
                 p = x$pval, group = "Experiments",
                 I2 = NA,
                 Qp = NA,
                 stringsAsFactors = FALSE)
    )})))
    df$predictor <- rep(names(x), each = nrow(df)/length(x))
    df
  }
  remresdf = rma2df(x)
  remresdf <- transform(remresdf, interval = CIUB - CILB)
  remresdf <- transform(remresdf, RelConf = 1/interval)
  if(!intercept){
    remresdf <- remresdf[remresdf$predictor!="(Intercept)",]
    remresdf$predictor <- factor(remresdf$predictor, levels=unique(remresdf$predictor))
    if(!is.null(labels)){
      remresdf$predictor <- factor(remresdf$predictor, labels=labels)
    }
  }
  #print(remresdf)
  p = ggplot(remresdf,
             aes(LogFC, factor(Study, c(rev(studies), "RE Model")), xmax=CIUB, xmin=CILB, shape = group, colour = predictor)) +
    #coord_cartesian(xlim=c(min(remresdf$CILB), max(remresdf$CIUB))) +
    #scale_alpha_discrete(range = c(0.2, 1)) +
    geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) +
    geom_errorbarh(alpha=0.5, color="black", height = 0.1) +
    geom_point(aes(size = RelConf)) +
    scale_size(range = c(2, 5), guide=FALSE) +
    geom_point(data = subset(remresdf, Study=="RE Model"), size=3) +
    theme_bw() +
    theme(text = element_text(size=16)) +
    facet_grid(group~predictor, scales= 'free', space = 'free_y') +
    labs(x = "Observed Outcome", y = NULL)+
    scale_colour_manual(values = palette)+
    guides(colour = F, shape = F)+
    scale_x_continuous(breaks = scales::pretty_breaks(3), labels = function(x) round(as.numeric(x), digits=2))+
    expand_limits(x = c(-0.02,0.02))
  if(hetero == "I2"){
    p <- p+geom_label(aes(x = Inf, y = -Inf,
                          hjust = 1,
                          vjust = 0,
                          label=ifelse(is.na(I2), I2, paste0("I\u00B2 = ",round(I2,2), "%"))),
                      #position = position_nudge(y = -0.4),
                      colour = 'black',label.size=0)
  }else if(hetero == "Qp"){
    p <- p+geom_label(aes(label=ifelse(is.na(Qp), NA, ifelse(Qp<0.001, "***",ifelse(Qp<0.01, "**", ifelse(Qp<0.05, "*", ""))))), position = position_nudge(y = 0.2), colour = 'black',label.size=0)
  }
  return(p)
}
