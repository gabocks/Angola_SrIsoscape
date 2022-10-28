################################################################################
###                  R-Code for the Manuscript                               ###
###  'A bioavailable strontium isoscape of Angola with implications			 ###
###      for the archaeology of the transatlantic slave trade'               ###
###		Submitted to JAS in October 2022									 ###
###                 Authors: Gaelle Bocksberger and Xueye Wang               ###
###                       Date: 08 August 2022                               ###
###                        R version 4.1.2                                   ###
###                                                                          ###
################################################################################

# install the necessary packages if not already installed
need_packages <- c("sf", "terra", "rnaturalearth", "tidyverse","readxl","VSURF","viridis","ranger","pdp","tuneRanger","ggpubr","parallel","doParallel","spatialRF","caret","ggcorrplot","ggspatial","classInt","ggExtra")

not_installed <- need_packages[!(need_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)                               # Install not installed packages

# load the necessary packages
invisible(lapply(need_packages, library, character.only = TRUE))

# set working directory
setwd()

# create a set of folders
dir.create("data") # to store the data, save the supplementary Table 2 here
dir.create("data/rasters") # save the raster files here downloaded from  https://drive.google.com/drive/folders/1v4zYwE33_E5DdeIQpnhBnQh_0IbDUItm?usp=sharing
dir.create("figures") # to store the plots
dir.create("results") # to store the tiff


# set seed for reproducibility
set.seed(1992)

#######################################################################################
## LOAD AND PREPARE THE PREDICTORS
#######################################################################################

# list the names of the predictors of interest
pred_cont_names<-c("agemax","agemean","agemin","ai_reproj","basement_age_reproj","bouger_reproj","dust_reproj","elevation_reproj","geol_agemax","geol_agemean","geol_agemin","MAP","mat_reproj","nfert","orc","pet_reproj","pfert","PhKCL","rcec_reproj","rclay_reproj","rm1_reproj","rph_reproj","sand","silt","srsrq1","srsrq3","ssa")            
pred_cat_names<-c("geol","glim_gum","glim_litho","glim_xx","lc")

# load angola raster data
pred_raster<-rast(list.files("data/rasters/",full.names=T,pattern=".tif$"))
pred_raster<-pred_raster[[c(pred_cont_names,pred_cat_names)]] 

#######################################################################################
## LOAD AND PREPARE THE SAMPLES
#######################################################################################
# load the samples data
samples<-read_excel("data/Wang et al. Supplementary Tables.xlsx",skip=2)
colnames(samples)[3:5]<-c("lat","long","Sriso")

# transform sample data into spatial object and project in Eck IV
samples_proj<-st_as_sf(x = samples,
                  coords = c("long", "lat"),
                  crs = st_crs(4326))%>%
				  st_transform(crs="+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
				  

# extract the predictor data for each sample location
# only continuous data
samples_cont<-terra::extract(pred_raster[[pred_cont_names]],vect(samples_proj),method='bilinear',na.rm=TRUE)
# only categorical data
samples_cat<-terra::extract(pred_raster[[pred_cat_names]],vect(samples_proj),method='simple',na.rm=TRUE)
samples_cat<-bind_cols(ID=samples_cat$ID,data.frame(lapply(samples_cat[,-1],as.factor)))

# join cat and continous data
samples_pred<-full_join(samples_cont,samples_cat)

# response name
response<-"Sriso"

############
## we have 1 NA in the lc data, so we replace it with the closest category
samples_pred[4,"lc"]<-"shrubland"

#######################################################################################
## VSURF: Variable Selection Using Random Forests 
## VSURF is a three steps variable selection procedure based on random forests
## Genuer, R. Poggi, J.-M. and Tuleau-Malot, C. (2015) https://journal.r-project.org/archive/2015-2/genuer-poggi-tuleaumalot.pdf
#######################################################################################

#dim(angola_samples_pred)
results.vsurf<-VSURF(samples_pred[,c(pred_cont_names,pred_cat_names)],pull(samples[,response]), RFimplem = "ranger", parallel = TRUE, ncores = detectCores() - 1, clusterType = "PSOCK") 

# variables chosen after threshold step
results.vsurf$varselect.thres
#[1]  7  8 13 24 15 27  6 20 23 22 16 19  5 11  9 10 18 14 17 12  4  2  1 26  3 28 32 25 21 31 29 30

# variables chosen after interpretation step
results.vsurf$varselect.interp
#[1] 7 8

# variables chosen after prediction step
results.vsurf$varselect.pred
# [1] 7 8

# best predictors names
colnames(samples_pred)[results.vsurf$varselect.pred+1]
#[1] "dust_reproj"      "elevation_reproj"


summary(results.vsurf)
plot(results.vsurf)

#######################################################################################
## CHECK FOR CORRELATION AND COLINEARITY 
#######################################################################################

# using auto.cor and auto.vif from spatialRF
# set the preference order of what we think are the most important variables
preference.order<-c("glim_gum","geol","agemax","geol_agemin","srsrq3","elevation_proj","basement_age_reproj",
					"PhKCL","mat_reproj","bouger_reproj","dust_reproj","glim_litho","MAP",
					"rcec_reproj","rclay_reproj","orc_reproj")

results.spatialRF <- spatialRF::auto_cor(
  x = samples_pred[, 2:33],
  cor.threshold = 0.75,
  #preference.order = preference.order
) %>% 
  spatialRF::auto_vif(
    vif.threshold = 5,
    #preference.order = preference.order
)

#[auto_cor()]: Removed variables: ssa, srsrq1, silt, sand, rph_reproj, rm1_reproj, rcec_reproj, pfert, pet_reproj, mat_reproj, MAP, geol_agemin, geol_agemean, elevation_reproj, agemin, agemean
#[auto_vif()]: Variables are not collinear.
#Warning messages:
#1: In spatialRF::auto_cor(x = samples_pred[, 2:32], cor.threshold = 0.75,  :
#  These columns are non-numeric and will be removed: geol, glim_gum, glim_litho, glim_xx
#2: In spatialRF::auto_cor(x = samples_pred[, 2:32], cor.threshold = 0.75,  :
#  These columns have zero variance and might cause issues: dust_reproj, pfert
#3: In spatialRF::auto_vif(., vif.threshold = 5, ) :
#  These columns have zero variance and might cause issues: dust_reproj
 

results.spatialRF$selected.variables
#[1] "agemax"              "ai_reproj"           "basement_age_reproj"
# [4] "bouger_reproj"       "dust_reproj"         "geol_agemax"        
# [7] "nfert"               "orc"                 "PhKCL"              
#[10] "rclay_reproj"        "srsrq3" 

#######################################################################################
# CREATE CORRELATION PLOT TO ADVISE ON PRED CHOICE
#######################################################################################
## make a correlation heatmap between all the predictors

# extract 10000 random points
ang_5000_cat<-spatSample(pred_raster,5000,na.rm=T,replace=F)

a<-model.matrix(~0+., data=ang_5000_cat) %>% 
  cor(use="pairwise.complete.obs") 

# keep only correlations higher than ±0.5  
a[a>(-0.5)& a<0.5]<-NA
  
ggcorrplot(a,show.diag = F, type="lower", lab=TRUE, lab_size=1)
  
ggsave("figures/correlation_predictors.pdf")

#######################################################################################
# VARIABLE SELECTION AND TRAINING/TESTING
#######################################################################################
# choose most adequate variables taking into account expert kknowledge as well as the results from VSURF, auto_cor and auto_vif
# As well as a visual check of the correlation matrix "correlation_predictors.pdf"
var.selected<-c("elevation_reproj","basement_age_reproj","srsrq3", "dust_reproj",
"geol_agemax","ai_reproj","rclay_reproj","rcec_reproj","orc","bouger_reproj"
,"glim_gum","mat_reproj")

# names of the variables for ploting
names.var.selected<-data.frame(r_names=var.selected,plot_names=c("elevation","terrane age","srsq3","dust","maximum age geology","aridity index","clay content",
					"cation exchange capacity","organic carbon content","bouguer anomaly","lithology","mean annual temperature")) 

# keep only the SELECTED variables and shuffle
training <- cbind(samples_pred[,var.selected],samples[,response])%>%droplevels()

# for reproducibility load the row.names
#training<-training[sample(1:nrow(training),replace=F), ]
#write.csv(data.frame(row_ID=rownames(training)),"data/rownames_training.csv",row.names=F)
row_id<-read.csv("data/rownames_training.csv")
training<-training[row_id$row_ID, ]

##############################################################################
## TUNE mtry AND RUN THE MODEL WITH caret::train(...,method="ranger") ON COMPLETE SET WITH K-FOLD CROSS VALIDATION
##############################################################################

# tune the mtry and min.node.size hyperparameters
samples.task = makeRegrTask(data = training, target = "Sriso")
res = tuneRanger(samples.task, measure = list(rmse),tune.parameters=c("mtry","min.node.size"))
res
#Recommended parameter settings: 
#  mtry min.node.size
#1    12             2
#Results: 
#         rmse exec.time
#1 0.006091133     0.254

tunegrid <- cbind(res$recommended.pars,data.frame(splitrule = "variance"))%>%select( "mtry", "splitrule", "min.node.size")

# for full reproducibility load the seeds for the each iteration
seeds_4model<-read.csv("data/seeds_4model.csv")

# set the cross validation controls
fitControl <- trainControl(
  method = "repeatedcv",# repeated cross validation
  number = 5, # 5-folds
  repeats=10,  # 10 times
  verboseIter=FALSE ,
  returnResamp="final",
  savePredictions="all",
  #allowParallel=TRUE,# With parallel backend
  seeds=as.list(seeds_4model$x))
  

# train the model with caret LOOCV cross-validation and ranger quantile random forest regression 
caret_model<-caret::train(y=training$Sriso, x =select(training,all_of(var.selected)),
						method = "ranger",
						importance="permutation",
						trControl= fitControl,
						metric="RMSE", 
						tuneGrid=tunegrid,
						quantreg = TRUE,
						keep.inbag=T,
						num.trees=500,
						seed=1992,
						respect.unordered.factors=T)

# if you get an invalid connection error, run the next 5 lines
#unregister <- function() {
# env <- foreach:::.foreachGlobals
#  rm(list=ls(name=env), pos=env)
#}
#unregister()

###############################################
## EXTRACT THE RESAMPLING (LOOCV) Rsquare and RMSE 
###############################################

# model R squared out of Bag
caret_model$results$Rsquared
# [1] 0.6681381
# model RMSE on OOB prediction error
caret_model$results$RMSE
#[1] 0.006379622

# standard deviation on the RMSE through the 10 repetitions
sd(caret_model$resample$RMSE)
#  0.001398254

##############################################################################
## GET THE QUANTILES PREDICTIONS FROM ALL THE TREES
##############################################################################

# extract the 2.5 and 97.5 quantiles of the predictions of all the repetitions, calculate the se as the difference between quantile 97.5 and 2.5
predict_se<-caret_model$pred%>%group_by(rowIndex)%>%summarise(predSr=mean(pred),q_5=quantile(pred, probs = c(0.5)),q_25=quantile(pred, probs = c(0.025)),q_975=quantile(pred, probs = c(0.975)),se=q_975-q_25)

# add to the predict data frame
predict_se$obsSr<-training[,response]

#######################
# create a data.frame with the R2 and RMSE results
r2_rmse<-caret_model$resample%>%summarize(RMSE_sd=sd(RMSE),RMSE=mean(RMSE),R2=mean(Rsquared))

##############################################################################	
# PLOT THE RESIDUALS DISTRIBUTION
##############################################################################

# plot obsSr against predSr with SE and correlation coefficient
					
lm_predobs<-ggplot(predict_se,aes(x=obsSr,y=predSr))+
geom_abline(slope=1,intercept=0,col="red",size=1)+
geom_point(alpha=0.5)+
#geom_point(alpha=0.5,aes(x=obsSr,y=predSr),col="orange")+
#geom_errorbar(aes(ymin=predSr-predict_se2$se/2,ymax=predSr+predict_se2$se/2),alpha=0.5,color="blue")+
geom_errorbar(aes(ymin=q_25,ymax=q_975),alpha=0.5)+ # to use the quantiles as "errorbar"
xlim(0.705,0.77)+ylim(0.705,0.77)+
xlab(expression(Observed~''^87*'Sr/'^86*'Sr'))+ylab(expression(Predicted~''^87*'Sr/'^86*'Sr'))+
geom_text(data = r2_rmse,size=10/.pt, aes(x = 0.725, y = 0.765, 
                 label = paste("R2 = ",round(R2,3),", ","RMSE = ",round(RMSE,4),"±",round(RMSE_sd,4), sep = " ")),show.legend = F)+
theme(legend.position="none")+theme_bw()

ggsave("figures/perVSobs_SE.pdf",lm_predobs)


# plot residuals against predicted Strontium
resid_plot<-ggplot(predict_se,aes(y=obsSr-predSr,x=predSr))+
geom_abline(slope=0,intercept=0,col="red",size=1)+
geom_point(alpha=0.5,size=2)+
ylab("Residuals")+xlab(expression(Predicted~''^87*'Sr/'^86*'Sr'))+
#ylim(c(-0.02,0.03))+
theme_bw()

resid_plot<-ggMarginal(resid_plot, type = "density",margins="y",size=7)

# get the percentage within ±0.005
1-(sum(abs(predict_se$obsSr-predict_se$predSr)>0.005)/101)
#[1] 0.6633663


ggsave("figures/residVSpred.pdf",resid_plot)

##############################################################################	
# VARIABLE IMPORTANCE PLOTS
##############################################################################

# plot the permutation score version

imp_4plot <- data.frame(var = row.names(data.frame(caret_model$finalModel$variable.importance)),
                   var_imp = caret_model$finalModel$variable.importance)%>%
				   mutate(var_names=names.var.selected$plot_names[match(var,names.var.selected$r_names)])%>%
				   arrange(var_imp)

							
var_imp<-ggplot(imp_4plot,aes(var_imp, x = reorder(var, var_imp))) +
  geom_point(size = 3, colour = "#ff6767") +
  coord_flip() +
  labs(x = "Predictors", y = "Permutation importance") +
  scale_x_discrete(labels=imp_4plot$var_names)+
  theme_bw()

ggsave("figures/variableImportance.pdf",var_imp)

ggsave("figures/Fig_5.pdf",ggarrange(var_imp,ggarrange(lm_predobs,resid_plot,nrow=1,labels=c("B","C")),nrow=2,labels=c("A","")),width=15,units="cm")

##############################################################################	
## PARTIAL DEPENDENCE PLOTS
##############################################################################
# one variable partial plot
# library(pdp) chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://journal.r-project.org/archive/2017/RJ-2017-016/RJ-2017-016.pdf
# https://christophm.github.io/interpretable-ml-book/pdp.html
# line is marginal effect of covariate on model prediction
# y axis is marginal effect

# partial plot of continuous predictors
par.df_cont<-NULL
pred.continuous<-var.selected[var.selected%in%pred_cont_names]

for (i in pred.continuous){
par.df <- pdp::partial(caret_model,train=training, pred.var = i, prob = TRUE)%>%# ice=T,center=T
			data.frame()%>%
			pivot_longer(cols=all_of(i),names_to="pred",values_to="value")

par.df_cont<-bind_rows(par.df_cont,par.df)
}		

par.df_cont$plot_names<-names.var.selected$plot_names[match(par.df_cont$pred,names.var.selected$r_names)]

rug_4plot<-select(training,any_of(pred.continuous))%>%
			pivot_longer(cols=everything(),values_to="value",names_to="pred")%>%
			mutate(plot_names=names.var.selected$plot_names[match(pred,names.var.selected$r_names)])
			
			
ggplot()+
geom_line(data=par.df_cont,aes(x=value,y=yhat),color="black")+#color="darkgrey",alpha=0.5 group=yhat.id
geom_rug(data=rug_4plot,aes(x=value),alpha=0.5)+
facet_wrap(~plot_names,scales="free")+
ylab(expression(Predicted~''^87*'Sr/'^86*'Sr'))+xlab("")

ggsave("figures/partial_dependence_continuous.pdf")


# partial plot of categorical predictors
par.df_cat<-NULL
pred.categorical<-var.selected[var.selected%in%pred_cat_names]

for (i in pred.categorical){
par.df <- pdp::partial(caret_model,train=training, pred.var = i, prob = TRUE)

par.df_cat<-bind_rows(par.df_cat,par.df)
}		

rug_4plot<-select(training,any_of(pred.categorical))%>%
			pivot_longer(cols=everything(),values_to="value",names_to="pred")%>%
			mutate(plot_names=names.var.selected$plot_names[match(pred,names.var.selected$r_names)])
			

ggplot()+
geom_boxplot(data=par.df_cat,aes(x=glim_gum,y=yhat))+
geom_rug(data=training,aes(x=glim_gum,y=0.71),alpha=0.5,position="jitter",sides="b")+
ylab("Marginal effect")+xlab("")+
ylim(c(0.72795,0.72865))


ggsave("figures/partial_dependence_categorical.pdf")


##############################################################################	
## PREDICT OVER GEOGRAPHIC AREA
##############################################################################

## extract and reproject the boundaries of Angola
angola<- ne_countries(country = 'angola', returnclass = "sf")
angola_proj<-st_transform(angola, "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## extract a raster and rename (only serves as a recipient for the predicted values)
prediction_raster<-pred_raster[[1]]
values(prediction_raster)<-NA # replace values with NA
names(prediction_raster)<-"predicted_Sr"

## predict the ranger model over a data.frame of the predictors
# create the data.frame of the values in the rasters
pred_df<-terra::as.data.frame(pred_raster[[names(pred_raster)%in%var.selected]],cells=T,na.rm=T)

# predict over the whole spatial area
predict_spatial<-predict.train(caret_model,newdata=pred_df[,-1],type="raw")

# place the predicted values back into a raster
values(prediction_raster)[pred_df$cell]<-predict_spatial
# aggregate the tiles for a lighter plot
prediction_raster_agg <- terra::aggregate(prediction_raster, fact = 2,na.rm=T)
# transform to wgs84
prediction_raster_wgs<-project(prediction_raster_agg,y="+proj=longlat +datum=WGS84")

# save the predicted map as tiff
writeRaster(prediction_raster,"results/angola_model_predictions.tif",overwrite=TRUE)

###
# extract 10 natural breaks of the predicted data
nat_breaks<-classIntervals(predict_spatial, n = 10, style = 'fisher')$brks

predict_plot_natbreaks<-ggplot()+
geom_sf(data=angola,fill=NA,color=NA)+
geom_tile(data=terra::as.data.frame(prediction_raster_wgs,xy=T),aes(x=x,y=y,fill=cut(predicted_Sr,nat_breaks)))+
scale_fill_viridis(expression(''^87*'Sr/'^86*'Sr'),direction=-1,option="mako",discrete=T,labels=round(nat_breaks,3))+
theme_bw()+
annotation_scale(width_hint=0.3)+
annotation_north_arrow(style = north_arrow_fancy_orienteering,location="tr",height=unit(0.7,"cm"),width=unit(0.4,"cm"))+
labs(x="",y="")+
theme(legend.key.width = unit(1, "line"))

ggsave("figures/angola_model_predictions_natbreaks.pdf",predict_plot_natbreaks,width=30,units="cm")

##############################################################################	
## compute the SE for assignR
##############################################################################

## extract a raster and rename (only serves as a recipient)
prediction_se_raster<-pred_raster[[1]]
names(prediction_se_raster)<-"predicted_Sr_SE"
values(prediction_se_raster)<-NA 

# predict SE
# because the vector becomes too heavy, it crashes, so I'm looping it and saving every 1000
pred_df_split<-split(pred_df, sample( nrow(pred_df),25, replace=F))

for (i in 1:25){

predict_se<-predict(caret_model$finalModel,pred_df_split[[i]],type="se",num.threads=3)$se

predict_se<-cbind(data.frame(se=predict_se),id=pred_df_split[[i]]$cell)

write.csv(predict_se,file=paste("results/temp_predict/predict_se_",i,".csv",sep=""), row.names=F)

}

all_predict_se<- list.files(path="results/temp_predict",pattern="predict_se", full.names = TRUE) %>% 
			lapply(read.csv) %>% 
			bind_rows%>%
			arrange(id)
 
# place the predicted values back into a raster
values(prediction_se_raster)[pred_df$cell]<-all_predict_se$se[match(pred_df$cell,all_predict_se$id)]
# transform to wgs84
prediction_se_raster_wgs<-project(prediction_se_raster,y="+proj=longlat +datum=WGS84")

# get the range of values for the colour scale
col_breaks_se=round(seq(from=round(min(all_predict_se$se),3)+0.001,to=round(max(all_predict_se$se),3)-0.001,length.out=5),3)

se_plot<-ggplot()+
geom_sf(data=angola,,col=NA,fill=NA)+
geom_tile(data=terra::as.data.frame(prediction_se_raster_wgs,xy=T),aes(x=x,y=y,fill=predicted_Sr_SE))+
scale_fill_viridis("predicted\nstandard\nerror",direction=-1,breaks=col_breaks_se,labels=col_breaks_se)+
theme_bw()+
annotation_scale(width_hint=0.3)+
annotation_north_arrow(location="tr",height=unit(0.8,"cm"),width=unit(0.5,"cm"),style = north_arrow_fancy_orienteering,)+
labs(x="",y="")+
theme(legend.position="bottom",
		legend.key.width = unit(2.5, "line"))

ggsave("figures/angola_model_SE.pdf",se_plot,width=15,units="cm")

# save the SE raster
writeRaster(prediction_se_raster,"results/angola_model_SE.tif",overwrite=TRUE)

