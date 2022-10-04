

#' Plotting result from coxtp() function
#'
#' @param fit Model get from coxtp
#' @param IC The Creteria selected for the plotting
#' @param coef The variable that needed to be plotted
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw theme element_text element_blank margin labs ggtitle
#' 
#' @return
#' @exportS3Method plot coxtp
#' 
#' @examples
# plot.coxtp <- function(model, IC="AIC", coef){
IC="AIC"
  if (!IC %in% c("AIC", "TIC", "GIC")) stop("IC has to be one of AIC, TIC and GIC!")
  # #Test uses
  # coef="V1"
  # IC="AIC"
  # fit=fit
  # ##
  if(IC == 'AIC'){
    fit <- model$model.AIC
  } else if(IC == 'TIC'){
    fit <- model$model.AIC
  } else{
    fit <- model$model.GIC
  }
  
  if (missing(model)) stop ("Argument model is required!")
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
  if (missing(expand)) expand <- c(1,1)/100
  ls.tvef <- confint.surtiver(fit, times, parm, level)$tvef
  if (length(ls.tvef)==0) stop("No time-varying effect chosen!")
  if (missing(labels)) labels <- names(ls.tvef)
  # if (!require(ggplot2)) install.packages('ggplot2')
  library(ggplot2)
  options(stringsAsFactors=F)
  
  if(is.null(fit$model.AIC)){
    model_plot<-fit$model_result
  } else {
    if(IC == 'AIC'){
      model_plot <- fit$model.AIC
    } else if(IC == 'TIC'){
      model_plot <- fit$model.TIC
    } else{
      model_plot <- fit$model.GIC
    }
  }


  B.spline <- as.matrix(model_plot$bases)
  theta_plot  <- model_plot$theta_list[[length(model_plot$theta_list)]]

  beta     <- B.spline%*%t(theta_plot)
  var      <-  model_plot$VarianceMatrix

  colnames(beta)=fit$z_names
  p        <- fit$p
  knot     <- dim(B.spline)[2]
  list     <- 1:dim(B.spline)[1]
  beta_low <- matrix(0,dim(B.spline)[1],p)
  beta_up  <- matrix(0,dim(B.spline)[1],p)
  for(i in 1:p){
    beta_t_1   <- beta[,i]
    var2       <- var[((i-1)*knot+1):((i-1)*knot+knot),((i-1)*knot+1):((i-1)*knot+knot)]
    temp       <- 1.96*sqrt(vapply(list, function(x) matrix(B.spline[x,],1,knot)%*%var2%*%t(matrix(B.spline[x,],1,knot)),FUN.VALUE=numeric(1)))
    low        <- beta_t_1-temp
    up         <- beta_t_1+temp

    beta_low[,i] <- low
    beta_up[,i]  <- up
  }

  colnames(beta_low) <- paste0(fit$z_names,"_low")
  colnames(beta_up) <- paste0(fit$z_names,"_up")

  beta              <- as.data.frame(beta)
  beta              <- cbind(beta, beta_low, beta_up)
  y=beta[[coef]]
  ymin=beta[[paste0(coef,"_low")]]
  ymax=beta[[paste0(coef,"_up")]]
  time=model_plot$uniqfailtimes


  plot<-ggplot(data=beta, aes(x=time)) +
    geom_line(aes(y= y),size = 0.9,color = 'red') +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill="red", alpha = 0.2) +
    scale_y_continuous(name=ylab) +
    theme_bw() +  theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(text= element_text(size=14)) + theme(axis.text= element_text(size=14)) +
    theme(axis.title.y = element_text(margin= margin(t=0, r=10, b=0, l=0))) +labs(x=xlab) +
    ggtitle(paste0("Effect of ",coef, " When holding other covariates constant"))
  plot
  return(plot)
}


