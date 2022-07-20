

#' Plotting result from coxtp() function
#'
#' @param fit Model get from coxtp
#' @param IC The Creteria selected for the plotting
#' @param coef The variable that needed to be plotted
#'
#' @return
#' @export
#'
#' @examples
coxtp.plot <- function(fit, IC="AIC", coef,xlab="Time",ylab="Hazard Ratio (log-scale)"){
  # if (missing(fit)) stop ("Argument fit is required!")
  # if (class(fit)!="surtvep") stop("Object fit is not of class 'surtvep'!")

  # #Test use
  # coef="V1"
  # IC="AIC"
  # fit=fit
  # ##
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


