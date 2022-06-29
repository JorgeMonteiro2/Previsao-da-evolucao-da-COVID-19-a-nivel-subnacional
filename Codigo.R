set.seed(1)
library(keras)
library (tensorflow)
library(forecast)
library (MLmetrics)
library (Boruta)
library(ranger)
library(rBayesianOptimization)
library(ggplot2)

set_random_seed (1)

#Select working directory
setwd(choose.dir())

#Select target country
country = "Portugal"

# Load Data ------------------------------------------------------------

#load cases
df <- read.csv("Cases.csv")
regions <- unique(df[df$ï..country == country ,2])

ginfo <-  c(1:4) 
univariate <- c(5)
regional <- c()
national <- c()

# differencing
df$differencing_1 <- 0
df$differencing_2 <- 0

for(i in unique(df$nuts_code)) {
	df_tmp <- df[df$nuts_code == i,]
	df$differencing_1[df$nuts_code == i][2:(nrow(df_tmp))] <-  df_tmp$rate_14_day_per_100k[2:nrow(df_tmp)] / df_tmp$rate_14_day_per_100k[1:(nrow(df_tmp)-1)] *100
}
for(i in unique(df$nuts_code)) {
	df_tmp <- df[df$nuts_code == i,]
	df$differencing_2[df$nuts_code == i][2:(nrow(df_tmp))] <-  df_tmp$differencing_1[2:nrow(df_tmp)] / df_tmp$differencing_1[1:(nrow(df_tmp)-1)] *100
}
univariate <- c(univariate, (ncol(df)-1):ncol(df) )

#One-hot encoding of target regions
for(i in 1:length(regions)){
new_column <- df[,2] == regions[i]
new_column[new_column == TRUE] <- 1
new_column[new_column == FALSE] <- 0
df <- data.frame(df, new_column)
colnames(df)[ncol(df)] <- regions[i]
regional <- c(regional,ncol(df)) }

df_tmp <- read.csv("Tests.csv")

df_tmp$year_week <- gsub("W","",as.character(df_tmp$year_week))
df$n_testing = 0
df$n_positivity = 0
national <- c(national, ncol(df)-1, ncol(df) )

df$r_testing = 0
df$r_positivity = 0
regional <- c(regional , ncol(df)-1, ncol(df) )

df$country_code = 0
ginfo <- c(ginfo, ncol(df))

for(i in 1:nrow(df_tmp)){
df$n_testing[(df$year_week == df_tmp$year_week[i]) & (df$ï..country == df_tmp$ï..country[i])] <- df_tmp$testing_rate[i]
df$n_positivity [(df$year_week == df_tmp$year_week[i]) & (df$ï..country == df_tmp$ï..country[i])] <- df_tmp$positivity_rate[i]
df$r_testing[(df$year_week == df_tmp$year_week[i]) & (df$nuts_code == df_tmp$region[i])] <- df_tmp$positivity_rate[i]
df$r_positivity [(df$year_week == df_tmp$year_week[i]) & (df$nuts_code == df_tmp$region[i])] <- df_tmp$positivity_rate[i]
}

df_tmp <- read.csv("Vaccine.csv")
df_tmp$YearWeek <- gsub("W","",as.character(df_tmp$YearWeek))
df_tmp <- df_tmp[df_tmp$TargetGroup == "ALL",]
df_tmp <- df_tmp[order(df_tmp$Region),]

Region <- df_tmp[1,11]
week <- df_tmp[1,15]
df_tmp2 <- df_tmp[0,]
for (i in 2:nrow(df_tmp) ){
	if(df_tmp$Vaccine[i] == "JANSS"){ df_tmp[i,8] <- df_tmp[i,8]+ df_tmp[i,6]}
	if (Region == df_tmp[i,11]){ 
	df_tmp[i,6] <- df_tmp[i,6] + df_tmp[i-1,6] 
	df_tmp[i,8] <- df_tmp[i,8] + df_tmp[i-1,8] 
	} 

if (week != df_tmp[i,15]) {df_tmp2 <- rbind(df_tmp2,df_tmp[i-1,])}

week <- df_tmp[i,15]
Region <- df_tmp[i,11]
}

df$r_first_dose <- 0
df$r_full_dose <- 0
regional <- c(regional , ncol(df)-1, ncol(df) )

for (i in 1:nrow(df_tmp2)) {
df$r_first_dose[(df$year_week == df_tmp2$YearWeek[i]) & (df$nuts_code == df_tmp2$Region[i])] <- df_tmp2$FirstDose[i]
df$r_full_dose[(df$year_week == df_tmp2$YearWeek[i]) & (df$nuts_code == df_tmp2$Region[i])] <- df_tmp2$SecondDose[i]
}

df$n_first_dose <- 0
df$n_full_dose <- 0
national <- c(national, ncol(df)-1, ncol(df) )
df$country_code <- substr(df$nuts_code,1,2)

for (i in 1:nrow(df_tmp2)) {
df$n_first_dose[(df$year_week == df_tmp2$YearWeek[i]) & (df$country_code == df_tmp2$Region[i])] <- df_tmp2$FirstDose[i] / df_tmp2$Population[i]
df$n_full_dose[(df$year_week == df_tmp2$YearWeek[i]) & (df$country_code == df_tmp2$Region[i])] <- df_tmp2$SecondDose[i] / df_tmp2$Population[i]
}


df_tmp <- read.csv("Variants.csv")
variants <- unique(df_tmp$variant)
for (i in 1:length(variants)) {
df <- data.frame(df, new = 0)
national <- c(national, ncol(df)) }
colnames(df)[(ncol(df)-length(variants)+1) : ncol(df)] <- variants

for (i in 1:nrow(df_tmp)) {
df[(df$year_week == df_tmp$year_week[i]) & (df$ï..country == df_tmp$ï..country[i]), colnames(df) == df_tmp$variant[i]] <- df_tmp$percent_variant[i]}
			
df[is.na(df)] <- 0

#Add regional mean, median and max
df$r_mean <- 0
df$r_median <- 0
df$r_max <- 0
regional <- c(regional, (ncol(df)-2) : ncol(df))

for (i in unique(df$ï..country)){
	df_tmp <- df[df$ï..country == i,]
		for (j in unique(df_tmp$year_week)){
			df_tmp2 <- df_tmp[df_tmp$year_week == j,]
				df$r_mean[df$ï..country == i & df$year_week == j] <- mean(df_tmp2$rate_14_day_per_100k)
				df$r_median[df$ï..country == i & df$year_week == j] <- median(df_tmp2$rate_14_day_per_100k)
				df$r_max[df$ï..country == i & df$year_week == j] <- max(df_tmp2$rate_14_day_per_100k)
}}

df_tmp <- read.csv("Global_Mobility_Report.csv")

df_tmp <- df_tmp[df_tmp$sub_region_1=="",]
df_tmp$date<- strftime(df_tmp$date, format = "%Y-W%V")

df_tmp2 <- df_tmp[0,]
for( i in 2:nrow(df_tmp)) {
if(df_tmp$date[i] != df_tmp$date[i+1]){df_tmp2 <- rbind(df_tmp2,df_tmp[i,])}
}
df_tmp2$date <- gsub("W","",as.character(df_tmp2$date))

gmr <- colnames(df_tmp2)[10:15]
for ( i in 1:length(gmr)) { df <- data.frame(df, new =0)}
colnames(df)[(ncol(df)-length(gmr)+1) : ncol(df)] <- gmr
df_tmp2 <- df_tmp2[df_tmp2$country_region_code %in% df$country_code,]

for (i in 1:nrow(df_tmp2)) {
df[(df$year_week == df_tmp2$date[i]) & (df$country_code == df_tmp2$country_region_code[i]), (ncol(df)-length(gmr)+1) : ncol(df)] <- df_tmp2[i, (ncol(df_tmp2)-length(gmr)+1) : ncol(df_tmp2)]
}
national <- c(national, (ncol(df)-5) : ncol(df) )


df$week  <- 0
ginfo <- c(ginfo, ncol(df))
for(i in 1:nrow(df)) {df$week[i] <- which(df$year_week[i] == unique(df$year_week))}
response_measures <- c()

df_tmp <- read.csv("response_graphs_data.csv")
measures <- unique(df_tmp$Response_measure)
for ( i in 1:length(measures)) { 
	df <- data.frame(df, new =0)
	response_measures <- c(response_measures, ncol(df))}
colnames(df)[(ncol(df)-length(measures)+1) : ncol(df)] <- measures
df_tmp <- df_tmp[df_tmp$Country %in% df$ï..country,]

df_tmp$date_start<- strftime(df_tmp$date_start, format = "%Y-%V")
df_tmp$date_end<- strftime(df_tmp$date_end, format = "%Y-%V")

df_tmp$date_end [is.na(df_tmp$date_end)] <- strftime(Sys.Date(), format = "%Y-%V")

for(i in 1:nrow(df_tmp)){ 	
	df[(df$ï..country == df_tmp$Country[i]) & (df$year_week >= df_tmp$date_start[i]), which(colnames(df) == df_tmp$Response_measure[i] )] <- 1
	df[(df$ï..country == df_tmp$Country[i]) & (df$year_week >= df_tmp$date_end[i]), which(colnames(df) == df_tmp$Response_measure[i] )] <- 0
}
national <- c(national, c( (ncol(df) - length(measures)+1) : ncol(df)))

# End of data loading -----------------------------------------------------

wforecast <- 2
timesteps <- 4
wval <- 3
wtrain <- 5

df <- df[,c(ginfo,univariate,regional,national)]
ginfo <- 1:length(ginfo)
univariate <- (length(ginfo)+1) :length(c(ginfo, univariate))
regional <- (length(c(ginfo,univariate))+1) :length(c(ginfo, univariate,regional))
national <- (length(c(ginfo,univariate,regional))+1) :length(c(ginfo, univariate,regional,national))
future <- c()

for ( i in c(regional,national)[which(c(regional,national) %in% which(!(colnames(df) %in% c("r_mean", "r_median","r_max"))) )]) {
	df$new = 0
	colnames(df)[ncol(df)] <- paste("f_", colnames(df)[i])
	future <- c(future,ncol(df)) }

f_measures <- 1:length(measures) + ncol(df) - length(measures)
df$target <- 0

for(i in unique(df$nuts_code)) {
	df_tmp <- df[df$nuts_code == i,]
	df$target[df$nuts_code == i][1:(nrow(df_tmp)-wforecast)] <- df_tmp$rate_14_day_per_100k[(1+wforecast):nrow(df_tmp)]
	df[df$nuts_code == i,][1:(nrow(df_tmp)-wforecast),future ] <- df_tmp[(1+wforecast):nrow(df_tmp), c(regional,national)]

}

df_tmp <- list()
for (i in 1:timesteps) {
	df_tmp[[i]] <- df[0,]
	for (j in unique(df$nuts_code) ){
		df_tmp[[i]] <- rbind(df_tmp[[i]], df[df$nuts_code == j,][(timesteps-i+1):(nrow(df[df$nuts_code == j,])-i+1),])
}
}

df_array <- array (0, dim = c(nrow(df_tmp[[1]]), timesteps, length(c(univariate,regional,national, future))  ))
for (i in timesteps:1) { df_array[,i,] <- as.matrix(df_tmp[[i]][,c(univariate,regional,national,future)]  ) }
df_array[is.na(df_array)] <- 0
df_array[!is.finite(df_array)] <-0
univariate <- univariate - length(ginfo)
regional <- regional - length(ginfo)
national <- national - length(ginfo)
future <- future - length(ginfo)
f_measures <- f_measures - length(ginfo)

df_info <- df_tmp[[1]][,c(ginfo, univariate[1])]
df_info$week_to_predict <- df_info$week + wforecast
observed <- df_tmp[[1]][,ncol(df_tmp[[1]])]
df_info$observed <- observed

first_predict <- min(df_info$week_to_predict[df_info$country_code== "PT"]) + wtrain + wval + wforecast -1
weeks_predict <- first_predict: (first_predict + 53)

df_info$forecast_Dumb <- 0
df_info$forecast_ARIMA <- 0
df_info$forecast_TBATS<- 0
df_info$forecast_Holt<- 0
df_info$forecast_CSsplines <- 0

df_info$forecast_GRU_U<- 0
df_info$forecast_GRU_R<- 0
df_info$forecast_GRU_N <- 0
df_info$forecast_GRU_F_Measures <- 0
df_info$forecast_GRU_F <- 0

df_info$forecast_PT_GRU_U<- 0
df_info$forecast_PT_GRU_R<- 0
df_info$forecast_PT_GRU_N <- 0
df_info$forecast_PT_GRU_F_Measures <- 0
df_info$forecast_PT_GRU_F <- 0

df_info$forecast_RF_U<- 0
df_info$forecast_RF_R<- 0
df_info$forecast_RF_N <- 0
df_info$forecast_RF_F_Measures <- 0
df_info$forecast_RF_F <- 0


prediction <- 0

for ( week in weeks_predict) {
	for (region in regions) {	
		df_info$forecast_Dumb[df_info$week_to_predict == week & df_info$region_name == region] <- df$rate_14_day_per_100k[df$week == (week-wforecast) & df$region_name == region]

		model<-auto.arima(df$rate_14_day_per_100k[df$week < (week-wforecast+1) & df$region_name == region])
		prediction <- forecast(model, h=wforecast)[[4]][wforecast] 
		df_info$forecast_ARIMA[df_info$week_to_predict == week & df_info$region_name == region] <- prediction
		df_info$forecast_ARIMA[df_info$forecast_ARIMA < 0] <- 0

		model<-tbats(df$rate_14_day_per_100k[df$week < (week-wforecast+1) & df$region_name == region])
		prediction <- forecast(model, h=wforecast)[[2]][wforecast]  
		df_info$forecast_TBATS[df_info$week_to_predict == week & df_info$region_name == region] <- prediction
		df_info$forecast_TBATS[df_info$forecast_TBATS< 0] <- 0
		
		model<-holt(df$rate_14_day_per_100k[df$week < (week-wforecast+1) & df$region_name == region])
		prediction <- forecast(model, h=wforecast)[[2]][wforecast]  
		df_info$forecast_Holt[df_info$week_to_predict == week & df_info$region_name == region] <- prediction
		df_info$forecast_Holt[df_info$forecast_Holt< 0] <- 0
		
		df_info$forecast_CSsplines [df_info$week_to_predict == week & df_info$region_name == region] <- splinef(df$rate_14_day_per_100k[df$week < (week-wforecast+1) & df$region_name == region], h=wforecast)[[5]][wforecast] 
		df_info$forecast_CSsplines [df_info$forecast_CSsplines < 0] <- 0

}}


#		-----------------------------GRU Predictions-----------------------------------
set.seed(1)
wval <- 4


keras_fit <- function(layer, units, drop, r_drop, var,lookback,wtrain){
set.seed(1)
set_random_seed (1)
units=round(units)
var = round(var)
wtrain =round (wtrain)

train <- which(df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval)))
i <- 0
while(!(   (length(train)-i)/batch_size == round((length(train)-i)/batch_size)   )   ){
i= i+1}
train <- train[1:(length(train)-i)]

variables <- t_variables[1:(var + length(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]]))]
x_train <-array(df_array[train,1:lookback ,variables], dim = c(length(train),lookback ,length(variables))) 
validation_data <- list(array(df_array[val,1:lookback ,variables], c(length(val), lookback , length(variables))) ,y=observed[val])

if(layer ==1) {
model_val<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}
if(layer ==2) {
model_val<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback ,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1) ) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1) ) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}
if(layer ==3) {
model_val <- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback ,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}

model_val%>% compile(
loss = 'mean_squared_error',
  optimizer = optimizer_adam(),
  metrics = c('mean_squared_error')  
)

history<-  model_val %>% fit(x = x_train,y=observed[train], 
	     validation_data = validation_data, 
           batch_size = batch_size,
           epochs     = 400, 
           verbose    = 0, 
           shuffle    = TRUE,
	     callback_early_stopping(monitor="val_loss", mode="min", verbose=0, patience=50, restore_best_weights = TRUE))

result <- list(Score = 1/min(history$metrics$val_loss), 
               Pred = 0)

return (result)
}

new_model <- c(round(quantile(1:length(weeks_predict))[1:4]))
update_model <- seq(from = 1 , to= length(weeks_predict), by= 2 )
update_model <- c(update_model, new_model)

batch_size = 128
t_variables <- c()

parameters <-data.frame()
for (i in 1:length(new_model)){
parameters <- rbind(parameters,data.frame(mode = c(0:5), new = weeks_predict[new_model[i]], b_layer = 0, b_units = 0, b_drop = 0, b_r_drop = 0, b_var =0 ,b_lookback=0,b_wtrain =3, best = 0))}
var_used <- list()
for (i in which(parameters$mode == 0)){ var_used[[i]] <- c(1)}
models <- list()
parameters$best[parameters$mode==0] <- 10000000000000000

for(mode in 1:5){

for (week in weeks_predict) { 

if(week %in% weeks_predict[new_model]){
if(mode == 1){ t_variables<- c(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]],univariate)}
if(mode == 2){ t_variables<- c(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]],regional)}
if(mode == 3){ t_variables<- c(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]],national)}
if(mode == 4){ t_variables<- c(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]],f_measures)}
if(mode == 5){ t_variables<- c(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]],future)}

k_clear_session()
layer <- 1
units <- 100
drop <- 0
var <- 0.5
r_drop <- 0
lookback <- 1
b_layer <- 1
b_units <- 100
b_drop <- 0
b_r_drop <- 0
b_var <- 1
b_lookback <- 1
best <- 0
test <- 1
wtrain <- 5
b_wtrain <- 2.5

train <- which(df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval)))
val <- which(df_info$week_to_predict %in% c(1:wval + week-wforecast-wval))

i <- 0
while(!(   (length(train)-i)/batch_size == round((length(train)-i)/batch_size)   )   ){
i= i+1}	
train <- train[1:(length(train)-i)]

i <- 0
while(!(   (length(val )-i)/batch_size == round((length(val )-i)/batch_size)   )   ){
i= i+1}	
val <- val [1:(length(val )-i)]

set.seed(1)
boruta <- Boruta(df_array[train,1,t_variables], y = observed[train])
t_variables <- t_variables[order(attStats(boruta)[,1], decreasing = TRUE)]
t_variables <- c(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]],t_variables[!t_variables %in% var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]]])

search_grid_keras <- data.frame(layer=0, units=0, drop=0, r_drop=0, var=0,lookback=0,wtrain=0, Value =0)

while(test < 9) {
layer <- b_layer
units <- b_units
drop <- b_drop
r_drop <- b_r_drop
var <- b_var
lookback <- b_lookback
wtrain <- b_wtrain

if(test == 1) { wtrain = wtrain * 2}
if(test == 2) { 
	var = round(var * 1.5)
	if(var + length(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]]) > length(t_variables)) {var = length(t_variables) - length(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]])}
	}
if(test == 3) { lookback = lookback +1 }
if(test == 4) { layer =2 }
if(test == 5) { layer =3 }
if(test == 6) { units = units * 2 }
if(test == 7) { drop = drop + 0.1 }
if(test == 8) { r_drop = r_drop + 0.1 }

run <- keras_fit (layer, units, drop, r_drop, var,lookback,wtrain)
search_grid_keras <- rbind(search_grid_keras,c(layer, units, drop, r_drop, var,lookback,wtrain, run$Score))

if ( run$Score > best){
best <- run$Score
b_layer <- layer
b_units <- units
b_drop <-drop
b_r_drop <- r_drop
b_var <- var
b_lookback <- lookback
b_history <- history
b_wtrain <- wtrain
	if(test == 6 & units == 200) { test = test+1 }
	if(test == 5) { 
		b_units = 25
		best = 0
		test = test + 1}
	if(test == 4){test = test + 1 }
	if( (test == 3) & (lookback == timesteps) ) { test = test + 1 }
	if( (test == 2) & (var + length(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]]) == length(t_variables)) ) { test = test + 1 }
	if( (test == 1) & (wtrain ==20) ){ test = test + 1 }
} else {
	if(test == 8) { test = test + 1 }
	if(test == 7) { test = test + 1 }
	if(test == 6) { test = test + 1 }
	if(test == 5) { 
		b_units = 25
		best = 0
		test = test + 1}
	if(test == 4) {test = test + 1 }
	if(test == 3) {test = test + 1 }
	if(test == 2) {test = test + 1 }
	if(test == 1) {test = test + 1 }
	}
}
search_grid_keras <- search_grid_keras[2:nrow(search_grid_keras),]

search_bound_keras <- list(layer = c(1L,3L), units=c(25,max(search_grid_keras$units)), drop = c(0,max(search_grid_keras$drop)*1.5), r_drop= c(0,max(search_grid_keras$r_drop)*1.5), var=c(1,length(t_variables) - length(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]])) ,lookback = c(1L,4L),wtrain = c(3,max(search_grid_keras$wtrain)))
if( search_bound_keras$drop[2] >1){ search_bound_keras$drop[2] =1}
if( search_bound_keras$r_drop[2] >1){ search_bound_keras$r_drop[2] =1}

bayes_keras <- BayesianOptimization(FUN = keras_fit, bounds = search_bound_keras, 
                     init_points = 0, init_grid_dt = search_grid_keras,
                     n_iter = 3, acq = "ucb")

layer <- bayes_keras$Best_Par[[1]]
units <- bayes_keras$Best_Par[[2]]
drop <- bayes_keras$Best_Par[[3]]
r_drop <- bayes_keras$Best_Par[[4]]     
var <- bayes_keras$Best_Par[[5]] 
lookback <- bayes_keras$Best_Par[[6]]  
wtrain <- bayes_keras$Best_Par[[7]]

set.seed(1)
set_random_seed (1)
units=round(units)
var = round(var)
wtrain =round (wtrain)

train <- which(df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval)))
i <- 0
while(!(   (length(train)-i)/batch_size == round((length(train)-i)/batch_size)   )   ){
i= i+1}
train <- train[1:(length(train)-i)]

variables <- t_variables[1:(var + length(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]]))]
x_train <-array(df_array[train,1:lookback ,variables], dim = c(length(train),lookback ,length(variables))) 
validation_data <- list(array(df_array[val,1:lookback ,variables], c(length(val), lookback , length(variables))) ,y=observed[val])

if(layer ==1) {
model<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}
if(layer ==2) {
model<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback ,length(variables)), dropout = drop , recurrent_dropout = r_drop , kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1) ) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}
if(layer ==3) {
model<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback ,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}

model%>% compile(
loss = 'mean_squared_error',
  optimizer = optimizer_adam(),
  metrics = c('mean_squared_error')  
)

history<-  model%>% fit(x = x_train,y=observed[train], 
	     validation_data = validation_data, 
           batch_size = batch_size,
           epochs     = 400, 
           verbose    = 0, 
           shuffle    = TRUE,
	     callback_early_stopping(monitor="val_loss", mode="min", verbose=0, patience=50, restore_best_weights = TRUE))
print(plot(history))
best <- min(history$metrics$val_loss)

	if(best >= parameters$best[which(parameters$mode==(mode-1) & parameters$new == week)]) {
	models[[which(parameters$mode==(mode) & parameters$new == week)]] <- models[[which(parameters$mode==(mode-1) & parameters$new == week)]]
	model <- models[[which(parameters$mode==(mode) & parameters$new == week)]]
	var_used[[which(parameters$mode==(mode) & parameters$new == week)]] <- var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]]
	variables <- var_used[[which(parameters$mode==(mode) & parameters$new == week)]]
	parameters[which(parameters$mode==(mode) & parameters$new == week), 3:10] <- parameters[which(parameters$mode==(mode-1) & parameters$new == week), 3:10]
	parameters$b_var[which(parameters$mode==(mode) & parameters$new == week)] <- 0
	lookback <- parameters$b_lookback[which(parameters$mode==(mode) & parameters$new == week)]
	wtrain <- parameters$b_wtrain[which(parameters$mode==(mode) & parameters$new == week)]
	} else {
		var_used[[which(parameters$mode==(mode) & parameters$new == week)]] <- variables
		parameters[which(parameters$mode==(mode) & parameters$new == week), 3:10] <- c(layer, units, drop, r_drop, var,lookback,wtrain, best)
		models[[which(parameters$mode==(mode) & parameters$new == week)]] <- model
		}
}

if(week %in% weeks_predict[update_model]) {

train <- which(df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval)))
val <- which(df_info$week_to_predict %in% c(1:wval + week-wforecast-wval))
i <- 0
while(!(   (length(train)-i)/batch_size == round((length(train)-i)/batch_size)   )   ){
i= i+1}	
train <- train[1:(length(train)-i)]
i <- 0
while(!(   (length(val )-i)/batch_size == round((length(val )-i)/batch_size)   )   ){
i= i+1}	
val <- val [1:(length(val )-i)]

x_train <-array(df_array[train,1:lookback ,variables], dim = c(length(train),lookback ,length(variables))) 
validation_data <- list(array(df_array[val,1:lookback ,variables], c(length(val), lookback , length(variables))) ,y=observed[val])
history<-  model%>% fit(x = x_train,y=observed[train], 
	     validation_data = validation_data, 
           batch_size = batch_size,
           epochs     = 400, 
           verbose    = 0, 
           shuffle    = TRUE,
	     callback_early_stopping(monitor="val_loss", mode="min", verbose=0, patience=50, restore_best_weights = TRUE))
print(week)
}

predict <- 1:batch_size
predict[(batch_size-length(regions)+1):batch_size] <- which(df_info$week_to_predict == week & df_info$region_name %in% regions)

x_predict <-array(df_array[predict,1:lookback ,variables], dim = c(length(predict),lookback ,length(variables))) 

df_info[which(df_info$week_to_predict == week & df_info$region_name %in% regions), 14+mode] <- predict(model, x_predict,batch_size = batch_size) [(batch_size-length(regions)+1):batch_size]
}
#Monitor Results
eval <- df_info [df_info$week_to_predict %in% weeks_predict & df_info$region_name %in% regions,]
results_RMSE <- data.frame(region = c(regions,"global"), DUMB=0, ARIMA = 0, TBATS = 0, holt = 0, CSsplines =0, GRU_U=0, GRU_R=0, GRU_N=0, GRU_F_Measures= 0,GRU_F = 0, RF_U=0, RF_R=0, RF_N=0, RF_F_Measures= 0,RF_F = 0, CNN_U=0, CNN_R=0, CNN_N=0, CNN_F_Measures= 0,CNN_F = 0)
for( i in 1:length(regions)){
for (j in 1:20){
results_RMSE[i,j+1] <- RMSE (eval$observed [eval$region_name == regions[i]], eval [eval$region_name == regions[i], j+length(ginfo)+3]  )}}
for (j in 1:20){ results_RMSE[length(regions)+1,j+1] <- RMSE (eval$observed [eval$region_name %in% regions], eval [eval$region_name %in% regions, j+length(ginfo)+3]  )}
print(results_RMSE)
}



#		-----------------------------  GRU_PT  -----------------------------------
set.seed(1)
wval <- 4


keras_fit <- function(layer, units, drop, r_drop, var,lookback,wtrain){
set.seed(1)
set_random_seed (1)
units=round(units)
var = round(var)
wtrain =round (wtrain)

train <- which( (df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval))) & (df_info$region_name %in% regions) )
i <- 0
while(!(   (length(train)-i)/batch_size == round((length(train)-i)/batch_size)   )   ){
i= i+1}
train <- train[1:(length(train)-i)]

variables <- t_variables[1:(var + length(var_used[[which(parameters$mode==(mode-1) & parameters$new == week)]]))]
x_train <-array(df_array[train,1:lookback ,variables], dim = c(length(train),lookback ,length(variables))) 
validation_data <- list(array(df_array[val,1:lookback ,variables], c(length(val), lookback , length(variables))) ,y=observed[val])

if(layer ==1) {
model_val<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}
if(layer ==2) {
model_val<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback ,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1) ) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1) ) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}
if(layer ==3) {
model_val <- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback ,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}

model_val%>% compile(
loss = 'mean_squared_error',
  optimizer = optimizer_adam(),
  metrics = c('mean_squared_error')  
)

history<-  model_val %>% fit(x = x_train,y=observed[train], 
	     validation_data = validation_data, 
           batch_size = batch_size,
           epochs     = 400, 
           verbose    = 0, 
           shuffle    = TRUE,
	     callback_early_stopping(monitor="val_loss", mode="min", verbose=0, patience=50, restore_best_weights = TRUE))

result <- list(Score = 1/min(history$metrics$val_loss), 
               Pred = 0)

return (result)
}

batch_size = length(regions)
t_variables <- c()

parameters_PT <-data.frame()
for (i in 1:length(new_model)){
parameters_PT <- rbind(parameters_PT,data.frame(mode = c(0:5), new = weeks_predict[new_model[i]], b_layer = 0, b_units = 0, b_drop = 0, b_r_drop = 0, b_var =0 ,b_lookback=0,b_wtrain =3, best = 0))}
var_used_PT <- list()
for (i in which(parameters_PT$mode == 0)){ var_used_PT[[i]] <- c(1)}
models_PT <- list()
parameters_PT$best[parameters_PT$mode==0] <- 10000000000000000

for(mode in 1:5){

for (week in weeks_predict) { 

if(week %in% weeks_predict[new_model]){
if(mode == 1){ t_variables<- c(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]],univariate)}
if(mode == 2){ t_variables<- c(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]],regional)}
if(mode == 3){ t_variables<- c(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]],national)}
if(mode == 4){ t_variables<- c(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]],f_measures)}
if(mode == 5){ t_variables<- c(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]],future)}

k_clear_session()
layer <- 1
units <- 100
drop <- 0
var <- 0.5
r_drop <- 0
lookback <- 1
b_layer <- 1
b_units <- 100
b_drop <- 0
b_r_drop <- 0
b_var <- 1
b_lookback <- 1
best <- 0
test <- 1
wtrain <- 5
b_wtrain <- 2.5

train <- which( (df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval))) & (df_info$region_name %in% regions) )
val <- which( (df_info$week_to_predict %in% c(1:wval + week-wforecast-wval)) & (df_info$region_name %in% regions) )

i <- 0
while(!(   (length(train)-i)/batch_size == round((length(train)-i)/batch_size)   )   ){
i= i+1}	
train <- train[1:(length(train)-i)]

i <- 0
while(!(   (length(val )-i)/batch_size == round((length(val )-i)/batch_size)   )   ){
i= i+1}	
val <- val [1:(length(val )-i)]

set.seed(1)
boruta <- Boruta(df_array[train,1,t_variables], y = observed[train])
t_variables <- t_variables[order(attStats(boruta)[,1], decreasing = TRUE)]
t_variables <- c(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]],t_variables[!t_variables %in% var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]]])

search_grid_keras <- data.frame(layer=0, units=0, drop=0, r_drop=0, var=0,lookback=0,wtrain=0, Value =0)

while(test < 9) {
layer <- b_layer
units <- b_units
drop <- b_drop
r_drop <- b_r_drop
var <- b_var
lookback <- b_lookback
wtrain <- b_wtrain

if(test == 1) { wtrain = wtrain * 2}
if(test == 2) { 
	var = round(var * 1.5)
	if(var + length(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]]) > length(t_variables)) {var = length(t_variables) - length(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]])}
	}
if(test == 3) { lookback = lookback +1 }
if(test == 4) { layer =2 }
if(test == 5) { layer =3 }
if(test == 6) { units = units * 2 }
if(test == 7) { drop = drop + 0.1 }
if(test == 8) { r_drop = r_drop + 0.1 }

run <- keras_fit (layer, units, drop, r_drop, var,lookback,wtrain)
search_grid_keras <- rbind(search_grid_keras,c(layer, units, drop, r_drop, var,lookback,wtrain, run$Score))

if ( run$Score > best){
best <- run$Score
b_layer <- layer
b_units <- units
b_drop <-drop
b_r_drop <- r_drop
b_var <- var
b_lookback <- lookback
b_history <- history
b_wtrain <- wtrain
	if(test == 6 & units == 200) { test = test+1 }
	if(test == 5) { 
		b_units = 25
		best = 0
		test = test + 1}
	if(test == 4){test = test + 1 }
	if( (test == 3) & (lookback == timesteps) ) { test = test + 1 }
	if( (test == 2) & (var + length(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]]) == length(t_variables)) ) { test = test + 1 }
	if( (test == 1) & (wtrain ==80) ){ test = test + 1 }
} else {
	if(test == 8) { test = test + 1 }
	if(test == 7) { test = test + 1 }
	if(test == 6) { test = test + 1 }
	if(test == 5) { 
		b_units = 25
		best = 0
		test = test + 1}
	if(test == 4) {test = test + 1 }
	if(test == 3) {test = test + 1 }
	if(test == 2) {test = test + 1 }
	if(test == 1) {test = test + 1 }
	}
}
search_grid_keras <- search_grid_keras[2:nrow(search_grid_keras),]

search_bound_keras <- list(layer = c(1L,3L), units=c(25,max(search_grid_keras$units)), drop = c(0,max(search_grid_keras$drop)*1.5), r_drop= c(0,max(search_grid_keras$r_drop)*1.5), var=c(1,length(t_variables) - length(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]])) ,lookback = c(1L,4L),wtrain = c(3,max(search_grid_keras$wtrain)))
if( search_bound_keras$drop[2] >1){ search_bound_keras$drop[2] =1}
if( search_bound_keras$r_drop[2] >1){ search_bound_keras$r_drop[2] =1}

bayes_keras <- BayesianOptimization(FUN = keras_fit, bounds = search_bound_keras, 
                     init_points = 0, init_grid_dt = search_grid_keras,
                     n_iter = 3, acq = "ucb")

layer <- bayes_keras$Best_Par[[1]]
units <- bayes_keras$Best_Par[[2]]
drop <- bayes_keras$Best_Par[[3]]
r_drop <- bayes_keras$Best_Par[[4]]     
var <- bayes_keras$Best_Par[[5]] 
lookback <- bayes_keras$Best_Par[[6]]  
wtrain <- bayes_keras$Best_Par[[7]]

set.seed(1)
set_random_seed (1)
units=round(units)
var = round(var)
wtrain =round (wtrain)

train <- which( (df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval))) & (df_info$region_name %in% regions) )
i <- 0
while(!(   (length(train)-i)/batch_size == round((length(train)-i)/batch_size)   )   ){
i= i+1}
train <- train[1:(length(train)-i)]

variables <- t_variables[1:(var + length(var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]]))]
x_train <-array(df_array[train,1:lookback ,variables], dim = c(length(train),lookback ,length(variables))) 
validation_data <- list(array(df_array[val,1:lookback ,variables], c(length(val), lookback , length(variables))) ,y=observed[val])

if(layer ==1) {
model<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}
if(layer ==2) {
model<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback ,length(variables)), dropout = drop , recurrent_dropout = r_drop , kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1) ) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}
if(layer ==3) {
model<- keras_model_sequential() %>%
  layer_gru(units = length(variables), activation = "relu", stateful = FALSE,return_sequences = TRUE, batch_input_shape = list(batch_size,lookback ,length(variables)), dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>%
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_gru(units = units, activation = "relu", stateful = FALSE,return_sequences = TRUE, dropout = drop , recurrent_dropout = r_drop, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1)) %>% 
  layer_dense(units = 1, kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 1))
}

model%>% compile(
loss = 'mean_squared_error',
  optimizer = optimizer_adam(),
  metrics = c('mean_squared_error')  
)

history<-  model%>% fit(x = x_train,y=observed[train], 
	     validation_data = validation_data, 
           batch_size = batch_size,
           epochs     = 400, 
           verbose    = 0, 
           shuffle    = TRUE,
	     callback_early_stopping(monitor="val_loss", mode="min", verbose=0, patience=50, restore_best_weights = TRUE))
print(plot(history))
best <- min(history$metrics$val_loss)

	if(best >= parameters_PT$best[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]) {
	models_PT[[which(parameters_PT$mode==(mode) & parameters_PT$new == week)]] <- models_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]]
	model <- models_PT[[which(parameters_PT$mode==(mode) & parameters_PT$new == week)]]
	var_used_PT[[which(parameters_PT$mode==(mode) & parameters_PT$new == week)]] <- var_used_PT[[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week)]]
	variables <- var_used_PT[[which(parameters_PT$mode==(mode) & parameters_PT$new == week)]]
	parameters_PT[which(parameters_PT$mode==(mode) & parameters_PT$new == week), 3:10] <- parameters_PT[which(parameters_PT$mode==(mode-1) & parameters_PT$new == week), 3:10]
	parameters_PT$b_var[which(parameters_PT$mode==(mode) & parameters_PT$new == week)] <- 0
	lookback <- parameters_PT$b_lookback[which(parameters_PT$mode==(mode) & parameters_PT$new == week)]
	wtrain <- parameters_PT$b_wtrain[which(parameters_PT$mode==(mode) & parameters_PT$new == week)]
	} else {
		var_used_PT[[which(parameters_PT$mode==(mode) & parameters_PT$new == week)]] <- variables
		parameters_PT[which(parameters_PT$mode==(mode) & parameters_PT$new == week), 3:10] <- c(layer, units, drop, r_drop, var,lookback,wtrain, best)
		models_PT[[which(parameters_PT$mode==(mode) & parameters_PT$new == week)]] <- model
		}
}

if(week %in% weeks_predict[update_model]) {

train <- which( (df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval))) & (df_info$region_name %in% regions) )
val <- which( (df_info$week_to_predict %in% c(1:wval + week-wforecast-wval)) & (df_info$region_name %in% regions) )
i <- 0
while(!(   (length(train)-i)/batch_size == round((length(train)-i)/batch_size)   )   ){
i= i+1}	
train <- train[1:(length(train)-i)]
i <- 0
while(!(   (length(val )-i)/batch_size == round((length(val )-i)/batch_size)   )   ){
i= i+1}	
val <- val [1:(length(val )-i)]

x_train <-array(df_array[train,1:lookback ,variables], dim = c(length(train),lookback ,length(variables))) 
validation_data <- list(array(df_array[val,1:lookback ,variables], c(length(val), lookback , length(variables))) ,y=observed[val])
history<-  model%>% fit(x = x_train,y=observed[train], 
	     validation_data = validation_data, 
           batch_size = batch_size,
           epochs     = 400, 
           verbose    = 0, 
           shuffle    = TRUE,
	     callback_early_stopping(monitor="val_loss", mode="min", verbose=0, patience=50, restore_best_weights = TRUE))
print(week)
}

predict <- 1:batch_size
predict[(batch_size-length(regions)+1):batch_size] <- which(df_info$week_to_predict == week & df_info$region_name %in% regions)

x_predict <-array(df_array[predict,1:lookback ,variables], dim = c(length(predict),lookback ,length(variables))) 

df_info[which(df_info$week_to_predict == week & df_info$region_name %in% regions), 19+mode] <- predict(model, x_predict,batch_size = batch_size) [(batch_size-length(regions)+1):batch_size]
}
#Monitor Results
eval <- df_info [df_info$week_to_predict %in% weeks_predict & df_info$region_name %in% regions,]
results_RMSE <- data.frame(region = c(regions,"global"), DUMB=0, ARIMA = 0, TBATS = 0, holt = 0, CSsplines =0, GRU_U=0, GRU_R=0, GRU_N=0, GRU_F_Measures= 0,GRU_F = 0, RF_U=0, RF_R=0, RF_N=0, RF_F_Measures= 0,RF_F = 0, CNN_U=0, CNN_R=0, CNN_N=0, CNN_F_Measures= 0,CNN_F = 0)
for( i in 1:length(regions)){
for (j in 1:20){
results_RMSE[i,j+1] <- RMSE (eval$observed [eval$region_name == regions[i]], eval [eval$region_name == regions[i], j+length(ginfo)+3]  )}}
for (j in 1:20){ results_RMSE[length(regions)+1,j+1] <- RMSE (eval$observed [eval$region_name %in% regions], eval [eval$region_name %in% regions, j+length(ginfo)+3]  )}
print(results_RMSE)
}



#		----------------------------------------- Random Forest ---------------------------------------------

set.seed(1)
new_model <- c(round(quantile(1:length(weeks_predict))[1:4]))
update_model <- seq(from = 1 , to= length(weeks_predict), by= 2 )
update_model <- c(update_model, new_model)

t_variables <- c()

parameters_RF <-data.frame()
for (i in 1:length(new_model)){
parameters_RF <- rbind(parameters_RF,data.frame(mode = c(0:5), new = weeks_predict[new_model[i]], b_mtry = 0,b_var = 1,b_wtrain =2, best = 0))}
var_used_RF <- list()
for (i in which(parameters_RF$mode == 0)){ var_used_RF[[i]] <- c(1)}
parameters_RF$best[parameters_RF$mode==0] <- 10000000000000000

for(mode in 1:5){

for (week in weeks_predict) { 

if(week %in% weeks_predict[new_model]){
if(mode == 1){ t_variables<- c(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]],univariate)}
if(mode == 2){ t_variables<- c(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]],regional)}
if(mode == 3){ t_variables<- c(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]],national)}
if(mode == 4){ t_variables<- c(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]],f_measures)}
if(mode == 5){ t_variables<- c(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]],future)}

b_mtry <- 1
b_var <- 0.5
b_wtrain <- 2.5
best <- 10000000000000000
test <- 1

wtrain <- 5
train <- which(df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval)))
val <- which(df_info$week_to_predict %in% c(1:wval + week-wforecast-wval))

boruta <- Boruta(df_array[train,1,t_variables], y = observed[train])
t_variables <- t_variables[order(attStats(boruta)[,1], decreasing = TRUE)]
t_variables <- c(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]],t_variables[!t_variables %in% var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]]])
while(test < 4) {

mtry <- b_mtry
var <- b_var
wtrain <- b_wtrain

if(test == 1) { wtrain = wtrain * 2}
if(test == 2) { 
	var = var * 2
	if(var + length(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]]) > length(t_variables)) {var = length(t_variables) - length(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]])}
	}
if(test == 3) { 
	mtry = mtry +1 
	if(mtry > var + length(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]])){ mtry <- var + length(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]]) } } 

train <- which(df_info$week_to_predict %in% c((week-wforecast-wval-wtrain+1):(week-wforecast-wval)))

variables <- t_variables[1:(var + length(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]]))]
x_train <- as.data.frame(df_array[train,1,variables])
x_validation <- as.data.frame(df_array[val,1,variables])
if(length(variables == 1)) { colnames(x_validation) <- colnames(x_train) }

print(c(mtry, var, wtrain,best))
print(variables)

model_val <- ranger(x = as.data.frame(x_train), y = observed[train], mtry = mtry, num.trees  = 2000)

if( RMSE (predict(model_val, data = x_validation)$prediction, observed[val])<best){
best <- RMSE (predict(model_val, data = x_validation)$prediction, observed[val])
b_var <- var
b_mtry <- mtry
b_wtrain <- wtrain
if( (test == 3) & (mtry == var + length(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]])) ) { test = test + 1 }
if( (test == 2) & (var + length(var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]]) == length(t_variables)) ) { test = test + 1 }
} else {
if(test == 3) {test = test+1}
if(test == 2) {test = test+1}
if(test == 1) {test = test+1}
}
}

if(best >= parameters_RF$best[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]) {
	var_used_RF[[which(parameters_RF$mode==(mode) & parameters_RF$new == week)]] <- var_used_RF[[which(parameters_RF$mode==(mode-1) & parameters_RF$new == week)]]
	variables <- var_used_RF[[which(parameters_RF$mode==(mode) & parameters_RF$new == week)]]
	parameters_RF[which(parameters_RF$mode==(mode) & parameters_RF$new == week), 3:6] <- parameters[which(parameters$mode==(mode-1) & parameters$new == week), 3:6]
	parameters_RF$b_var[which(parameters_RF$mode==(mode) & parameters_RF$new == week)] <- 0
	mtry <- parameters_RF$b_mtry[which(parameters_RF$mode==(mode) & parameters_RF$new == week)]
	wtrain <- parameters_RF$b_wtrain[which(parameters_RF$mode==(mode) & parameters_RF$new == week)]
	} else {
		var_used_RF[[which(parameters_RF$mode==(mode) & parameters_RF$new == week)]] <- variables
		parameters_RF[which(parameters_RF$mode==(mode) & parameters_RF$new == week), 3:6] <- c(b_mtry, b_var,b_wtrain, best)
	}
}

if(week %in% weeks_predict[update_model]) {
train <- which(df_info$week_to_predict %in% c((week-wforecast-wtrain+1):(week-wforecast)))
x_train <-as.data.frame(df_array[train,1,variables]) 
model <- ranger(x = as.data.frame(x_train), y = observed[train], mtry = mtry, num.trees  = 2000)
print(week)
}

predict <- which(df_info$week_to_predict == week & df_info$region_name %in% regions)
x_predict <- as.data.frame(df_array[predict,1,variables]) 
if(length(variables == 1)) { colnames(x_predict) <- colnames(x_train) }

df_info[which(df_info$week_to_predict == week & df_info$region_name %in% regions), 24+mode] <- predict(model, data = x_predict)$prediction
}}


#		------------------------------------ RESULTS ----------------------

rsq <- function (x, y) cor(x, y) ^ 2

eval <- df_info [df_info$week_to_predict %in% weeks_predict & df_info$region_name %in% regions,]

results_RMSE <- data.frame(region = c(regions,"global"), DUMB=0, ARIMA = 0, TBATS = 0, holt = 0, CSsplines =0, GRU_U=0, GRU_R=0, GRU_N=0, GRU_F_Measures= 0,GRU_F = 0, PT_U=0, PT_R=0, PT_N=0, PT_F_Measures= 0,PT_F = 0,RF_U=0, RF_R=0, RF_N=0, RF_F_Measures= 0,RF_F = 0)
for( i in 1:length(regions)){
for (j in 1:20){
results_RMSE[i,j+1] <- RMSE (eval$observed [eval$region_name == regions[i]], eval [eval$region_name == regions[i], j+length(ginfo)+3]  )}}
for (j in 1:20){ results_RMSE[length(regions)+1,j+1] <- RMSE (eval$observed [eval$region_name %in% regions], eval [eval$region_name %in% regions, j+length(ginfo)+3]  )}
results_RMSE

write.csv(results_RMSE,"RMSE.csv")

results_MSE <- data.frame(region = c(regions,"global"), DUMB=0, ARIMA = 0, TBATS = 0, holt = 0, CSsplines =0, GRU_U=0, GRU_R=0, GRU_N=0, GRU_F_Measures= 0,GRU_F = 0, PT_U=0, PT_R=0, PT_N=0, PT_F_Measures= 0,PT_F = 0,RF_U=0, RF_R=0, RF_N=0, RF_F_Measures= 0,RF_F = 0)
for( i in 1:length(regions)){
for (j in 1:20){
results_MSE[i,j+1] <- MSE(eval$observed [eval$region_name == regions[i]], eval [eval$region_name == regions[i], j+length(ginfo)+3]  )}}
for (j in 1:20){ results_MSE[length(regions)+1,j+1] <- MSE(eval$observed [eval$region_name %in% regions], eval [eval$region_name %in% regions, j+length(ginfo)+3]  )}
results_MSE

write.csv(results_MSE,"MSE.csv")

t_test <- data.frame(method = c("DUMB", "ARIMA", "TBATS", "holt", "CSsplines", "GRU_U", "GRU_R", "GRU_N", "GRU_F_Measures","GRU_F", "PT_U", "PT_R", "PT_N", "PT_F_Measures","PT_F", "RF_U", "RF_R", "RF_N", "RF_F_Measures","RF_F"), DUMB=0, ARIMA = 0, TBATS = 0, holt = 0, CSsplines =0, GRU_U=0, GRU_R=0, GRU_N=0, GRU_F_Measures= 0,GRU_F = 0, PT_U=0, PT_R=0, PT_N=0, PT_F_Measures= 0,PT_F = 0, RF_U=0, RF_R=0, RF_N=0, RF_F_Measures= 0,RF_F = 0)
for( i in 1:20){
for (j in 1:20){
t_test[i,j+1] <-t.test((eval$observed - eval[,i+9])^2, (eval$observed - eval[,j+9])^2 )[[3]] }}
t_test

write.csv(t_test,"t_test_MSE.csv")

results_R <- data.frame(region = c(regions,"global"), DUMB=0, ARIMA = 0, TBATS = 0, holt = 0, CSsplines =0, GRU_U=0, GRU_R=0, GRU_N=0, GRU_F_Measures= 0,GRU_F = 0, PT_U=0, PT_R=0, PT_N=0, PT_F_Measures= 0,PT_F = 0,RF_U=0, RF_R=0, RF_N=0, RF_F_Measures= 0,RF_F = 0)
for( i in 1:length(regions)){
for (j in 1:20){
results_R[i,j+1] <- rsq (eval$observed [eval$region_name == regions[i]], eval [eval$region_name == regions[i], j+length(ginfo)+3]  )}}
for (j in 1:20){ results_R[length(regions)+1,j+1] <- rsq (eval$observed [eval$region_name %in% regions], eval [eval$region_name %in% regions, j+length(ginfo)+3]  )}
results_R

write.csv(results_R,"R.csv")

results_MAE <- data.frame(region = c(regions,"global"), DUMB=0, ARIMA = 0, TBATS = 0, holt = 0, CSsplines =0, GRU_U=0, GRU_R=0, GRU_N=0, GRU_F_Measures= 0,GRU_F = 0, PT_U=0, PT_R=0, PT_N=0, PT_F_Measures= 0,PT_F = 0,RF_U=0, RF_R=0, RF_N=0, RF_F_Measures= 0,RF_F = 0)
for( i in 1:length(regions)){
for (j in 1:20){
results_MAE[i,j+1] <- MAE (eval$observed [eval$region_name == regions[i]], eval [eval$region_name == regions[i], j+length(ginfo)+3]  )}}
for (j in 1:20){ results_MAE[length(regions)+1,j+1] <- MAE (eval$observed [eval$region_name %in% regions], eval [eval$region_name %in% regions, j+length(ginfo)+3]  )}
results_MAE

write.csv(results_MAE,"MAE.csv")

t_test <- data.frame(method = c("DUMB", "ARIMA", "TBATS", "holt", "CSsplines", "GRU_U", "GRU_R", "GRU_N", "GRU_F_Measures","GRU_F", "PT_U", "PT_R", "PT_N", "PT_F_Measures","PT_F", "RF_U", "RF_R", "RF_N", "RF_F_Measures","RF_F"), DUMB=0, ARIMA = 0, TBATS = 0, holt = 0, CSsplines =0, GRU_U=0, GRU_R=0, GRU_N=0, GRU_F_Measures= 0,GRU_F = 0, PT_U=0, PT_R=0, PT_N=0, PT_F_Measures= 0,PT_F = 0, RF_U=0, RF_R=0, RF_N=0, RF_F_Measures= 0,RF_F = 0)
for( i in 1:20){
for (j in 1:20){
t_test[i,j+1] <-t.test(abs(eval$observed - eval[,i+9]), abs(eval$observed - eval[,j+9]))[[3]] }}
t_test

write.csv(t_test,"t_test_MAE.csv")

for (i in unique (eval$region_name)){eval$year_week[eval$region_name == i] <- df_info$year_week[14:67]}

for (i in unique (eval$region_name)){
tmp <- eval[eval$region_name == i,]
tmp2 <- data.frame(NUTS2 = i, Ano_Semana = tmp$year_week, Método = "Observado", taxa_de_incidência = tmp$observed )
tmp2 <- rbind(tmp2, data.frame(NUTS2 = i, Ano_Semana = tmp$year_week, Método = "CSsplines", taxa_de_incidência = tmp$forecast_CSsplines))
tmp2 <- rbind(tmp2, data.frame(NUTS2 = i, Ano_Semana = tmp$year_week, Método = "GRU", taxa_de_incidência = tmp$forecast_GRU_N))
tmp2 <- rbind(tmp2, data.frame(NUTS2 = i, Ano_Semana = tmp$year_week, Método = "RF", taxa_de_incidência = tmp$forecast_RF_U))
tmp2 <- rbind(tmp2, data.frame(NUTS2 = i, Ano_Semana = tmp$year_week, Método = "NP", taxa_de_incidência = tmp$forecast_Dumb))

filename <- paste0("results\\",i,".png")
img <- ggplot(data = tmp2,
aes(x= Ano_Semana, y= taxa_de_incidência, size = Método, group = Método,colour = Método)) + geom_line(aes(linetype=Método)) + facet_wrap(vars(NUTS2)) + theme_bw() + 
scale_linetype_manual(values=c("solid", "solid","solid", "solid", "dotted"))+
scale_colour_manual(values=c("lightsteelblue2","coral1", "red","black","blue"))+ 
scale_size_manual(values=c(0.7,0.7,0.7,0.7,1)) + 
ylab("Incidência a 14 dias por 100 Mil Habitantes")

ggplot2::ggsave(filename = filename,img, width = 30, height = 7)
}

var_used[[22]]

v <- c()
for (i in which(parameters$mode ==5)) {
v<- c(v, var_used[[i]])}
v <- unique (v)
used <- colnames(df)[v+max(ginfo)]

used


