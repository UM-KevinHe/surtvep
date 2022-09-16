

data(support)

head(support)

event = support$death
time  = support$d.time

z = as.matrix(as.numeric(as.factor(support$sex)))
colnames(z) = "gender"
fit = coxtp(event = event, time = time, z = z)

coxtp.plot(fit, coef = "gender")


table(support$income)
sum(table(support$income))


data <- support
data$income_cat=array(NA,length(support[,1]))
data[is.na(data$income),]$income = "NA"
data[data$income=="NA",]$income_cat = "NA"
data[data$income=="under $11k",]$income_cat = "income_under11k"
data[data$income=="$11-$25k",]$income_cat = "income_11_25k"
data[data$income=="$25-$50k",]$income_cat = "income_25_50k"
data[data$income==">$50k",]$income_cat = "income_above50k"
data$income_cat = factor(data$income_cat, levels = c("NA", "income_under11k", "income_11_25k", "income_25_50k", "income_above50k"))

income_NA  = (data$income_cat == "NA")
income_under11k = (data$income_cat == "income_under11k")
income_11_25k = (data$income_cat == "income_11_25k")
income_25_50k = (data$income_cat == "income_25_50k")
income_above50k = (data$income_cat == "income_above50k")

z  = cbind(income_NA,income_under11k,income_11_25k,income_25_50k,income_above50k)

fit = coxtp(event = event, time = time, z = z, ICLastOnly = TRUE, lambda_spline = 10, btr = "dynamic")
coxtp.plot(fit, coef = ""income_11_25k"")


Rcpp::sourceCpp("src/PenalizeStopCpp.cpp")







