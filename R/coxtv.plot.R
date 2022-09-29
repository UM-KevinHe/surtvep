

#' Plotting result from coxtp() function
#'
#' @param fit Model get from coxtp
#' @param IC The Creteria selected for the plotting
#' @param coef The variable that needed to be plotted
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw theme element_text element_blank margin labs ggtitle
#' 
#' @return
#' @export
#' 
#' @examples
plot.surtiver <- function(fit, times, parm, CI=TRUE, level=0.95, exponentiate=FALSE, 
                          xlab, ylab, xlim, ylim, save=FALSE, allinone=FALSE, 
                          title, linetype, fill, color, labels, expand, ...) {
  
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="coxtv") stop("Object fit is not of class 'coxtv'!")
  # if (!is.logical(save)) stop("Invalid save!")
  # if (!is.logical(exponentiate)) stop("Invalid exponentiate!")
  term.event <- attr(fit, "response")
  if (missing(xlab)) xlab <- "time"
  missingxlim <- missing(xlim); missingylim <- missing(ylim); 
  missingtitle <- missing(title); missinglty <- missing(linetype)
  missingfill <- missing(fill); missingcolor <- missing(color)
  defaultcols <- c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
  defaultltys <- c("solid", "dashed", "dotted", "dotdash", "longdash")
  if (missing(expand)) expand <- c(1,1)/100
  ls.tvef <- confint.surtiver(fit, times, parm, level)$tvef
  if (length(ls.tvef)==0) stop("No time-varying effect chosen!")
  if (missing(labels)) labels <- names(ls.tvef)
  if (!require(ggplot2)) install.packages('ggplot2')
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
        ggtitle(paste0(tv, " (", term.event, ")")) + theme_bw() +
        theme(plot.title=element_text(hjust=0),
              panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(), panel.border=element_blank(),
              axis.line=element_line(color="black"),
              axis.title=element_text(size=18, margin=margin(t=0,r=0,b=0,l=0)),
              axis.text=element_text(size=14), text=element_text(size=14),
              legend.title=element_blank(), legend.text=element_text(size=14),
              legend.position=c(0.5, 1), legend.box="horizontal")
    })
    if (save) {
      return(ls.plts)
    } else {
      return(ggpubr::ggarrange(plotlist=ls.plts, common.legend=T, ...))
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
    if (save) {
      return(plt)
    } else {
      plt
    }
  }
  
}