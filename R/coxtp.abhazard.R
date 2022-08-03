coxtp.abhazard=function(fit,var){
  
  
  #Test use, delete later
  sim_data=sim_data_p5
  #Fit the model first:
  event=sim_data[,"event"]
  time=sim_data[,"time"]
  data=sim_data[,!colnames(sim_data) %in% c("event","time")]
  
  # lambda_spline_all=c(0.001,0.01,0.1,1,10,100,1000)
  # fit_penalized <- coxtp(event = event, z = data, time = time,lambda_spline=lambda_spline_all)
  lambda_spline=1000
  fit_penalized <- coxtp(event = event, z = data, time = time,lambda_spline=lambda_spline)
  
  model1  = fit_penalized$model_result
  
  baseline=coxtp.baseline(fit=model1, delta=event,z=data,time=time)
  data_predict=c(1,0,0,0,0)
  predict=coxtp.predict(model1,baseline,newdata=data_predict)
  
  
}