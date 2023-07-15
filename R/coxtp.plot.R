#' plotting results from a fitted `coxtp` object
#' 
#' This function creates a plot of the time-varying coefficients from a fitted `coxtp` model. 
#'
#' @param fit model obtained from `coxtp`.
#' @param parm covariate name fitted in the model to be plotted. If `NULL`, all covariates are plotted.
#' @param CI if `TRUE`, confidence intervals are displayed. Default value is `TRUE`.
#' @param level the level of confidence interval. Default value is `0.95`.
#' @param exponentiate if `TRUE`, exponential scale of the fitted coefficients (hazard ratio) for each covariate is plotted. 
#' If `FALSE`, the fitted time-varying coefficients (log hazard ratio) are plotted.
#' @param xlim the limits for the x axis.
#' @param ylim the limits for the y axis.
#' @param xlab the title for the x axis.
#' @param ylab the title for the y axis.
#' @param title the title for the plot.
#' @param linetype the line type for the plot.
#' @param color the aesthetics parameter for the plot.
#' @param fill the aesthetics parameter for the plot.
#' @param times the time points for which the time-varying coefficients to be plotted. 
#' The default value is the unique observed event times in the dataset fitting the time-varying effects model.
#' @param allinone if `TRUE`, the time-varying trajectories for different covariates are combined into a single plot. Default value is `FALSE`.
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw theme element_text element_blank margin labs ggtitle
#' 
#' @exportS3Method plot coxtp
#' 
#' @examples
# data(ExampleData)
# z <- ExampleData$x
# time <- ExampleData$time
# event <- ExampleData$event
# fit <- coxtp(event = event, z = z, time = time)
# plot(fit$lambda1)
#' 
plot.coxtp <- function(fit, parm, CI=TRUE, level=0.95, exponentiate=FALSE, 
                       xlab, ylab, xlim, ylim, allinone=FALSE, 
                       title, linetype, color, fill, times) {
  
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="coxtp") stop("Object fit is not of class 'coxtp'!")
  # if (!is.logical(save)) stop("Invalid save!")
  # if (!is.logical(exponentiate)) stop("Invalid exponentiate!")
  term.event <- attr(fit, "response")
  if (missing(xlab)) xlab <- "time"
  if (missing(ylab)) ylab <- ifelse(exponentiate,"hazard ratio","coefficient")
  missingxlim <- missing(xlim); missingylim <- missing(ylim); 
  missingtitle <- missing(title); missinglty <- missing(linetype)
  missingfill <- missing(fill); missingcolor <- missing(color)
  defaultcols <- c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
  defaultltys <- c("solid", "dashed", "dotted", "dotdash", "longdash")
  expand <- c(1,1)/100
  ls.tvef <- confint(fit, times, parm, level)$tvef
  if (length(ls.tvef)==0) stop("No time-varying effect chosen!")
  # if (missing(labels)) labels <- names(ls.tvef)
  labels <- names(ls.tvef)
  # if (!require(ggplot2)) install.packages('ggplot2')
  library(ggplot2)
  options(stringsAsFactors=F)
  
  if (!allinone) {
    ls.plts <- lapply(names(ls.tvef), function(tv) {
      df.tv <- data.frame(ls.tvef[[tv]], as.numeric(rownames(ls.tvef[[tv]])))
      names(df.tv) <- c("est", "lower", "upper", "time")
      for (col in names(df.tv)) {
        range.tmp <- range(df.tv[!is.infinite(df.tv[,col]),col])
        df.tv[is.infinite(df.tv[,col]) & df.tv[,col] < 0, col] <- range.tmp[1]
        df.tv[is.infinite(df.tv[,col]) & df.tv[,col] > 0, col] <- range.tmp[2]
      }
      row.names(df.tv) <- NULL
      if (exponentiate) df.tv[,-4] <- exp(df.tv[,-4])
      plt <- ggplot(data=df.tv, aes(x=time)) +
        geom_hline(yintercept=ifelse(exponentiate,1,0),
                   color="black", size=0.3, linetype="dashed") +
        geom_line(aes(y=est, linetype="estimate"), size=0.9)
      if (CI) {
        plt <- plt +
          geom_ribbon(aes(ymin=lower, ymax=upper,
                          fill=paste0(round(100*level),"% CI")), alpha=0.4)
      }
      if (missingxlim) {
        plt <- plt + scale_x_continuous(name=xlab, expand=expand)
      } else {
        if (!is.numeric(xlim)) stop("Invalid xlim!")
        plt <- plt + scale_x_continuous(name=xlab, expand=expand, limits=xlim)
      }
      if (missingylim) {
        plt <- plt +
          scale_y_continuous(name=ylab, expand=expand)
      } else {
        if (!is.numeric(ylim)) stop("Invalid ylim!")
        plt <- plt +
          scale_y_continuous(name=ylab, expand=expand, limits=ylim)
      }
      plt +
        scale_linetype_manual("", values="solid") +
        scale_fill_manual("", values="grey") +
        # ggtitle(paste0(tv, " (", term.event, ")")) + theme_bw() +
        # ggtitle(paste0(tv)) + theme_bw() +
        theme(plot.title=element_text(hjust=0),
              panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(), panel.border=element_blank(),
              axis.line=element_line(color="black"),
              axis.title=element_text(size=18, margin=margin(t=0,r=0,b=0,l=0)),
              axis.text=element_text(size=14), text=element_text(size=14),
              legend.title=element_blank(), legend.text=element_text(size=14),
              legend.position=c(0.5, 1), legend.box="horizontal")
    })
    if (missingtitle) {
      return(ggpubr::ggarrange(plotlist=ls.plts, common.legend=T))
    } else {
      final.plt <- ggpubr::ggarrange(plotlist=ls.plts, common.legend=T)
      return(ggpubr::annotate_figure(final.plt, top = ggpubr::text_grob( title, face = "bold", size = 14)))
    }
  } else {
    if (length(names(ls.tvef)) > 5) stop("Number of parameters greater than 5!")
    df <- do.call(rbind, lapply(names(ls.tvef), function(tv) {
      df.tv <- data.frame(ls.tvef[[tv]], as.numeric(rownames(ls.tvef[[tv]])))
      names(df.tv) <- c("est", "lower", "upper", "time")
      for (col in names(df.tv)) {
        range.tmp <- range(df.tv[!is.infinite(df.tv[,col]),col])
        df.tv[is.infinite(df.tv[,col]) & df.tv[,col] < 0, col] <- range.tmp[1]
        df.tv[is.infinite(df.tv[,col]) & df.tv[,col] > 0, col] <- range.tmp[2]
      }
      if (exponentiate) df.tv[,-4] <- exp(df.tv[,-4])
      df.tv[,"parm"] <- tv
      row.names(df.tv) <- NULL
      return(df.tv)}))
    plt <- ggplot(data=df, aes(x=time, group=parm)) +
      geom_hline(yintercept=ifelse(exponentiate,1,0),
                 color="black", size=0.3, linetype="dashed") +
      geom_line(aes(y=est, linetype=parm, color=parm), size=0.9)
    if (CI) {
      plt <- plt +
        geom_ribbon(aes(ymin=lower, ymax=upper, fill=parm), alpha=0.1)
    }
    if (missingxlim) {
      plt <- plt + scale_x_continuous(name=xlab, expand=expand)
    } else {
      if (!is.numeric(xlim)) stop("Invalid xlim!")
      plt <- plt + scale_x_continuous(name=xlab, expand=expand, limits=xlim)
    }
    if (missingylim) {
      plt <- plt +
        scale_y_continuous(name=ylab, expand=expand)
    } else {
      if (!is.numeric(ylim)) stop("Invalid ylim!")
      plt <- plt +
        scale_y_continuous(name=ylab, expand=expand, limits=ylim)
    }
    if (!missinglty) {
      plt <- plt + scale_linetype_manual(NULL, values=linetype, labels=labels)
    } else {
      plt <- plt + scale_linetype_manual(NULL, values=defaultltys[1:length(names(ls.tvef))], 
                                         labels=labels)
    }
    if (!missingcolor) {
      plt <- plt + scale_color_manual(NULL, values=color, labels=labels)
    } else {
      plt <- plt + scale_color_manual(NULL, values=defaultcols[1:length(names(ls.tvef))],
                                      labels=labels)
    }
    if (!missingfill & CI) {
      plt <- plt + scale_fill_manual(NULL, values=fill, labels=labels)
    } else if (CI) {
      plt <- plt + scale_fill_manual(NULL, values=defaultcols[1:length(names(ls.tvef))],
                                     labels=labels)
    }
    if (!missingtitle) plt <- plt + ggtitle(title)
    plt <- plt + guides(linetype=guide_legend(nrow=1), fill=guide_legend(nrow=1), 
                        color=guide_legend(nrow=1)) +
      theme_bw() +
      theme(plot.title=element_text(hjust=0),
            panel.background=element_blank(), panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), panel.border=element_blank(),
            axis.line=element_line(color="black"),
            axis.title=element_text(size=18, margin=margin(t=0,r=0,b=0,l=0)),
            axis.text=element_text(size=14), text=element_text(size=14),
            legend.title=element_blank(), legend.text=element_text(size=14),
            legend.position=c(0.5, 1), legend.box="horizontal")
    # if (save) {
      return(plt)
    # } else {
      # plt
    # }
  }
  
}





#' #' Plotting result from a `ic.coxtp` object
#' #'
#' #' @param fit Fitted "coxtp" model object from the main function "coxtp".
#' #' @param IC The Criteria selected for the plot. Default is "AIC", which uses AIC to select the tunning parameter.
#' #' @param times The time scale for the plot
#' #' @param CI This argument controls the confidence interval on the plot. With 'CI' = TRUE, the confidence interval will be plotted.
#' #' @param level The confidence level.
#' #'
#' #' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw theme element_text element_blank margin labs ggtitle
#' #'
#' #' @exportS3Method plot ic.coxtp
#' #'
#' #' @examples
#' #' data("ExampleData")
#' #' z     <- ExampleData$x
#' #' time  <- ExampleData$time
#' #' event <- ExampleData$event
#' #' model   <- coxtp(event = event, z = z, time = time, lambda_spline = c(1))
#' #' plot(model, IC = "AIC", allinone = TRUE)
#' #'
#' plot.ic.coxtp <- function(model, IC="AIC", times, CI=TRUE, level=0.95, exponentiate=FALSE,
#'                        xlab, ylab, xlim, ylim, save=FALSE, allinone=FALSE,
#'                        title, linetype, fill, color, labels, expand, ...){
#' 
#'   if (missing(model)) stop ("Argument fit is required!")
#'   if (class(model)!="ic.coxtp") stop("Object fit is not of class 'ic.coxtp'!")
#' 
#'   if (!IC %in% c("AIC", "TIC", "GIC")) stop("IC has to be one of AIC, TIC and GIC!")
#'   # #Test uses
#'   # coef="V1"
#'   # IC="AIC"
#'   # fit=fit
#'   # ##
#'   if(IC == 'AIC'){
#'     fit <- model$model.AIC
#'   } else if(IC == 'TIC'){
#'     fit <- model$model.AIC
#'   } else{
#'     fit <- model$model.GIC
#'   }
#' 
#'   if (!is.logical(save)) stop("Invalid save!")
#'   if (!is.logical(exponentiate)) stop("Invalid exponentiate!")
#'   term.event <- attr(fit, "response")
#'   if (missing(xlab)) xlab <- "time"
#'   if (missing(ylab)) ylab <- ifelse(exponentiate,"hazard ratio","coefficient")
#'   missingxlim <- missing(xlim); missingylim <- missing(ylim);
#'   missingtitle <- missing(title); missinglty <- missing(linetype)
#'   missingfill <- missing(fill); missingcolor <- missing(color)
#'   defaultcols <- c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
#'   defaultltys <- c("solid", "dashed", "dotted", "dotdash", "longdash")
#'   if (missing(expand)) expand <- c(1,1)/100
#'   # ls.tvef <- confint.coxtp(fit, times, parm, level)$tvef
#'   ls.tvef <- confint(fit, times, level)$tvef
#' 
#'   if (length(ls.tvef)==0) stop("No time-varying effect chosen!")
#'   if (missing(labels)) labels <- names(ls.tvef)
#'   # if (!require(ggplot2)) install.packages('ggplot2')
#'   library(ggplot2)
#'   options(stringsAsFactors=F)
#' 
#'   if (!allinone) {
#'     ls.plts <- lapply(names(ls.tvef), function(tv) {
#'       df.tv <- data.frame(ls.tvef[[tv]], as.numeric(rownames(ls.tvef[[tv]])))
#'       names(df.tv) <- c("est", "lower", "upper", "time")
#'       for (col in names(df.tv)) {
#'         range.tmp <- range(df.tv[!is.infinite(df.tv[,col]),col])
#'         df.tv[is.infinite(df.tv[,col]) & df.tv[,col] < 0, col] <- range.tmp[1]
#'         df.tv[is.infinite(df.tv[,col]) & df.tv[,col] > 0, col] <- range.tmp[2]
#'       }
#'       row.names(df.tv) <- NULL
#'       if (exponentiate) df.tv[,-4] <- exp(df.tv[,-4])
#'       plt <- ggplot(data=df.tv, aes(x=time)) +
#'         geom_hline(yintercept=ifelse(exponentiate,1,0),
#'                    color="black", size=0.3, linetype="dashed") +
#'         geom_line(aes(y=est, linetype="estimate"), size=0.9)
#'       if (CI) {
#'         plt <- plt +
#'           geom_ribbon(aes(ymin=lower, ymax=upper,
#'                           fill=paste0(round(100*level),"% CI")), alpha=0.4)
#'       }
#'       if (missingxlim) {
#'         plt <- plt + scale_x_continuous(name=xlab, expand=expand)
#'       } else {
#'         if (!is.numeric(xlim)) stop("Invalid xlim!")
#'         plt <- plt + scale_x_continuous(name=xlab, expand=expand, limits=xlim)
#'       }
#'       if (missingylim) {
#'         plt <- plt +
#'           scale_y_continuous(name=ylab, expand=expand)
#'       } else {
#'         if (!is.numeric(ylim)) stop("Invalid ylim!")
#'         plt <- plt +
#'           scale_y_continuous(name=ylab, expand=expand, limits=ylim)
#'       }
#'       plt +
#'         scale_linetype_manual("", values="solid") +
#'         scale_fill_manual("", values="grey") +
#'         ggtitle(paste0(tv, " (", term.event, ")", " ", IC)) + theme_bw() +
#'         theme(plot.title=element_text(hjust=0),
#'               panel.background=element_blank(), panel.grid.major=element_blank(),
#'               panel.grid.minor=element_blank(), panel.border=element_blank(),
#'               axis.line=element_line(color="black"),
#'               axis.title=element_text(size=18, margin=margin(t=0,r=0,b=0,l=0)),
#'               axis.text=element_text(size=14), text=element_text(size=14),
#'               legend.title=element_blank(), legend.text=element_text(size=14),
#'               legend.position=c(0.5, 1), legend.box="horizontal")
#'     })
#'     if (save) {
#'       return(ls.plts)
#'     } else {
#'       return(ggpubr::ggarrange(plotlist=ls.plts, common.legend=T, ...))
#'     }
#'   } else {
#'     if (length(names(ls.tvef)) > 5) stop("Number of parameters greater than 5!")
#'     df <- do.call(rbind, lapply(names(ls.tvef), function(tv) {
#'       df.tv <- data.frame(ls.tvef[[tv]], as.numeric(rownames(ls.tvef[[tv]])))
#'       names(df.tv) <- c("est", "lower", "upper", "time")
#'       for (col in names(df.tv)) {
#'         range.tmp <- range(df.tv[!is.infinite(df.tv[,col]),col])
#'         df.tv[is.infinite(df.tv[,col]) & df.tv[,col] < 0, col] <- range.tmp[1]
#'         df.tv[is.infinite(df.tv[,col]) & df.tv[,col] > 0, col] <- range.tmp[2]
#'       }
#'       if (exponentiate) df.tv[,-4] <- exp(df.tv[,-4])
#'       df.tv[,"parm"] <- tv
#'       row.names(df.tv) <- NULL
#'       return(df.tv)}))
#'     plt <- ggplot(data=df, aes(x=time, group=parm)) +
#'       geom_hline(yintercept=ifelse(exponentiate,1,0),
#'                  color="black", size=0.3, linetype="dashed") +
#'       geom_line(aes(y=est, linetype=parm, color=parm), size=0.9)
#'     if (CI) {
#'       plt <- plt +
#'         geom_ribbon(aes(ymin=lower, ymax=upper, fill=parm), alpha=0.1)
#'     }
#'     if (missingxlim) {
#'       plt <- plt + scale_x_continuous(name=xlab, expand=expand)
#'     } else {
#'       if (!is.numeric(xlim)) stop("Invalid xlim!")
#'       plt <- plt + scale_x_continuous(name=xlab, expand=expand, limits=xlim)
#'     }
#'     if (missingylim) {
#'       plt <- plt +
#'         scale_y_continuous(name=ylab, expand=expand)
#'     } else {
#'       if (!is.numeric(ylim)) stop("Invalid ylim!")
#'       plt <- plt +
#'         scale_y_continuous(name=ylab, expand=expand, limits=ylim)
#'     }
#'     if (!missinglty) {
#'       plt <- plt + scale_linetype_manual(NULL, values=linetype, labels=labels)
#'     } else {
#'       plt <- plt + scale_linetype_manual(NULL, values=defaultltys[1:length(names(ls.tvef))],
#'                                          labels=labels)
#'     }
#'     if (!missingcolor) {
#'       plt <- plt + scale_color_manual(NULL, values=color, labels=labels)
#'     } else {
#'       plt <- plt + scale_color_manual(NULL, values=defaultcols[1:length(names(ls.tvef))],
#'                                       labels=labels)
#'     }
#'     if (!missingfill & CI) {
#'       plt <- plt + scale_fill_manual(NULL, values=fill, labels=labels)
#'     } else if (CI) {
#'       plt <- plt + scale_fill_manual(NULL, values=defaultcols[1:length(names(ls.tvef))],
#'                                      labels=labels)
#'     }
#'     if (!missingtitle) plt <- plt + ggtitle(title)
#'     plt <- plt + guides(linetype=guide_legend(nrow=1), fill=guide_legend(nrow=1),
#'                         color=guide_legend(nrow=1)) +
#'       theme_bw() +
#'       theme(plot.title=element_text(hjust=0),
#'             panel.background=element_blank(), panel.grid.major=element_blank(),
#'             panel.grid.minor=element_blank(), panel.border=element_blank(),
#'             axis.line=element_line(color="black"),
#'             axis.title=element_text(size=18, margin=margin(t=0,r=0,b=0,l=0)),
#'             axis.text=element_text(size=14), text=element_text(size=14),
#'             legend.title=element_blank(), legend.text=element_text(size=14),
#'             legend.position=c(0.5, 1), legend.box="horizontal")
#'     if (save) {
#'       return(plt)
#'     } else {
#'       plt
#'     }
#'   }
#' 
#' 
#' }


