# Master degree project  
## Title: Explore the biogeographic characteristics of microbiome on diverse datasets and predict the geographical source of the microbiome  

#### Course: Binp51  
#### Credits: 45  
#### Student: Yali Zhang  
#### Supervisor: Eran Elhaik  

This readme file contains work flow and all the scripts used in research. The required publicly available input files are provided in *Data_file* folder on Github (https://github.com/YaliZhang98/Binp51_master_degree_project).

## 1 Construction of mGPS interface
All the scripts and functions in this part were performed in R (version 4.0.3).  
The mGPS interface was built based on the “shiny” package (version 1.7.1).  
The mGPS implementation consists of three main components: feature elimination, prediction model construction, and sample prediction using the built model. All these procedures and visualization of results were included in the following script.  
In addition, the code file (mGPS_interface.r) for launching the mGPS interface is provided separately in *mGPS_interface* folder on Github, along with the corresponding readme files and tutorial. It also has a separate Github storage address: https://github.com/YaliZhang98/mGPS_interface

```r
# Loading the required packages
library(geosphere)
library(sp)
library(rworldmap)
library(caret)
library(maps)
library(geosphere)
library(caret)
library(plyr)
library(rgeos)
library(mapplots)
library(shiny)
library(dplyr)
library(xgboost)
library(e1071)
library(rgdal)
library(maptools)
library(sp)

setwd(rprojroot::find_rstudio_root_file())

# Part 1 Function ----------------------------

## Feature selection algorithm that used in feature elimination step
species_select <-
  function(x,
           y,
           remove_correlated = T,
           subsets = NULL,
           cores = 1) {
    doParallel::registerDoParallel(cores)
    
    y <- factor(y)
    
    print("feature elimination start --------------------")
    
    if (remove_correlated == T) {
      correlated <- findCorrelation(
        cor(x, use = "complete.obs"),
        cutoff = 0.98,
        verbose = FALSE,
        names = FALSE,
        exact = FALSE
      )
      
      x <- x[,-c(correlated)]
      print(paste0("correlated features removed:",length(correlated)))
      
    }
    
    len <- ncol(x)
    if(is.null(subsets)){
      subsets <-
        c(floor(len/2), floor(len / 4), floor(len / 8), floor(len / 16), floor(len / 32),floor(len / 64))
    }
    
    rfe_ctrl <- rfeControl(
      functions = rfFuncs,
      method = "cv",
      number =  5,
      verbose = FALSE,
      allowParallel = TRUE
    )
    set.seed(123)
    featureElimination <- rfe(
      x = x,
      y = y,
      sizes = subsets,
      rfeControl = rfe_ctrl,
      tuneLength = 2
    )
    doParallel::registerDoParallel(1) 
    
    print("feature elimination has been done ...")
    
    return(featureElimination)
  }


## Main mGPS algorithm: prediction model construction
mGPS <- 
  function(training = NULL,
           testing = NULL,
           classTarget,
           hierarchy = c('continent','city','latitude','longitude'),
           variables,
           nthread = 1,
           coast = NULL) {
    
    #Check for training set
    if(is.null(training)){
      return(message("No training set given"))
    }else{
      
      training <- droplevels(training)
      
      if (length(hierarchy)==3){
        class_target_orig_name <- training[,hierarchy[1]]
        class_target_make_name <- make.names(training[,hierarchy[1]])
                                           
        training[,hierarchy[1]] <- make.names(training[,hierarchy[1]])
        training[,hierarchy[1]] <- factor(training[,hierarchy[1]])
      }else if(length(hierarchy)==4){
        
        if (classTarget == hierarchy[1]){
          class_target_orig_name <- training[,hierarchy[1]]
          class_target_make_name <- make.names(training[,hierarchy[1]])
        }else if (classTarget == hierarchy[2]){
          class_target_orig_name <- training[,hierarchy[2]]
          class_target_make_name <- make.names(training[,hierarchy[2]])
        }
        
        training[,hierarchy[1]] <- make.names(training[,hierarchy[1]])
        training[,hierarchy[1]] <- factor(training[,hierarchy[1]])
        training[,hierarchy[2]] <- make.names(training[,hierarchy[2]])
        training[,hierarchy[2]] <- factor(training[,hierarchy[2]])
      }
      
      # Train mGPS with 5-fold cross validation of training set for hyperparameter tuning. 
      message("Training mGPS...")
      
      set.seed(1234)
      
      folds <- createFolds(training[,classTarget], k = 5, returnTrain = T)
      
      trControlClass <-  trainControl(
        method = "cv",
        number = 5,  
        verboseIter = FALSE,
        returnData = FALSE,
        search = "grid",
        savePredictions = "final",
        classProbs = T,
        allowParallel = T,
        index = folds )
      
      trControl <-  trainControl(
        method = "cv",
        number = 5,  
        verboseIter = FALSE,
        returnData = FALSE,
        search = "grid",
        savePredictions = "final",
        allowParallel = T,
        index = folds)
      
      tune_grid <- expand.grid(
        nrounds = c(300,600),
        eta = c( 0.05, 0.1),
        max_depth = c(3,6,9),
        gamma = 0,
        colsample_bytree = c(0.6,0.8),
        min_child_weight = c(1),
        subsample = (0.7)
      )
      
      if(length(hierarchy) == 4){
        
        print("--------- hierarchy 1 training ---------")

        Xgb_region <- train(x = training[,variables],y = training[,hierarchy[1]],
                            method = "xgbTree",
                            trControl = trControlClass,
                            tuneGrid = tune_grid,
                            nthread = nthread)
        
        print("hierarchy 1 training has been done ...")
        
        l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,hierarchy[1]]) ])
      }else{
        l1_train <- training[,c(variables)]
      }

      print("--------- hierarchy 2 training --------- ")
      
      if (length(hierarchy)==3){
        class_train <- hierarchy[1]
      }else if(length(hierarchy)==4){
        class_train <- hierarchy[2]
      }
      
      Xgb_class <- train(x = l1_train,y = training[,class_train],
                         method = "xgbTree",
                         trControl = trControlClass,
                         tuneGrid = tune_grid,
                         nthread = nthread)
      
      print("hierarchy 2 training has been done ...")
      
      l2_train <- data.frame(l1_train,Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])
      
      print("--------- latitude training ---------")
      
      if (length(hierarchy)==3){
        lat <- hierarchy[2]
      }else if(length(hierarchy)==4){
        lat <- hierarchy[3]
      }
      
      Xgb_latitude <- train(x = l2_train ,y = training[,lat],
                            method = "xgbTree",
                            trControl = trControl,
                            tuneGrid = tune_grid,
                            nthread = nthread)
      print("Latitude training has been done ...")
      
      l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
      
      print("--------- longitude training ---------")
      
      if (length(hierarchy)==3){
        long <- hierarchy[3]
      }else if(length(hierarchy)==4){
        long <- hierarchy[4]
      }
      
      Xgb_longitude <- train(x = l3_train ,y = training[,long],
                             method = "xgbTree",
                             trControl = trControl,
                             tuneGrid = tune_grid,
                             nthread = nthread)
      
      print("Longitude training has been done ...")
    }
    
    #check for test set, return trained model if no test set. 
      model <- function(test,variables){
        
        if(length(hierarchy) == 4){
          regPred <- predict(Xgb_region, newdata = test[,variables])
          regProbs <- predict(Xgb_region, newdata = test[,variables],type ="prob")
          
          l1_test <- data.frame(test[,variables], regProbs)
        }else{
          l1_test <- test[,variables]
        }
        classPred <- predict(Xgb_class, newdata = l1_test)
        classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
        
        l2_test <-  data.frame(l1_test, classProbs) 
        latPred <- predict(Xgb_latitude, newdata = l2_test)
        
        l3_test <- data.frame(l2_test, latPred)
        longPred <- predict(Xgb_longitude, newdata = l3_test)
        
        if (length(hierarchy)==3){
          
          classPred_new <- c()
          for (p in classPred){
            posi <- which(class_target_make_name %in% p)
            classPred_new <- c(classPred_new,class_target_orig_name[posi])
          }
            
          return(list(classPred_new, latPred, longPred))
        }else if(length(hierarchy)==4){
          
          if (classTarget == hierarchy[1]){
            
            regPred_new <- c()
            for (p in regPred){
              posi <- which(class_target_make_name %in% p)
              regPred_new <- c(regPred_new,class_target_orig_name[posi])
            }
            
            return(list(regPred_new, latPred, longPred))
          }else if (classTarget == hierarchy[2]){
            
            classPred_new <- c()
            for (p in classPred){
              posi <- which(class_target_make_name %in% p)
              classPred_new <- c(classPred_new,class_target_orig_name[posi])
            }
            
            return(list(classPred_new, latPred, longPred))
          }
        }
      }

## Sample prediction using the built model: generate mGPS predictions for test set
      
      message("Generating predictions...")

      features <- variables
      newdata_in <- as.data.frame(setNames(replicate(length(features),numeric(0), simplify = F), features))
      for (f in features){
        if (f %in% colnames(testing)){
          newdata_in[c(1:dim(testing)[1]),f] <- testing[,f]
        }else{
          newdata_in[c(1:dim(testing)[1]),f] <- data.frame(0)
        }
      }
      
      if(length(hierarchy) == 4){
        regPred <- predict(Xgb_region, newdata = newdata_in)
        regProbs <- predict(Xgb_region, newdata = newdata_in,type ="prob")
        
        l1_test <- data.frame(newdata_in, regProbs)
      }else{
        l1_test <- newdata_in
      }
      
      classPred <- predict(Xgb_class, newdata = l1_test)
      classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
      
      l2_test <-  data.frame(l1_test, classProbs) 
      latPred <- predict(Xgb_latitude, newdata = l2_test)
      
      l3_test <- data.frame(l2_test, latPred)
      longPred <- predict(Xgb_longitude, newdata = l3_test)
      
      print("Prediction has been done ...")
      
      #adjust out of bounds predictions
      longPred[longPred > 180] <- 180
      longPred[longPred < -180] <- -180
      latPred[latPred > 90] <- 90
      latPred[latPred < -90] <- -90
      
      if(length(hierarchy) == 3){
        model_store <- list(Xgb_class,Xgb_latitude,Xgb_longitude,"model" = model)
      }else  if(length(hierarchy) == 4){
        model_store <- list(Xgb_region,Xgb_class,Xgb_latitude,Xgb_longitude,"model" = model)
      }
      
      if (length(hierarchy)==3){
        prediction_store <- list(classPred, latPred, longPred)
      }else if(length(hierarchy)==4){
        if (classTarget == hierarchy[1]){
          prediction_store <- list(regPred, latPred, longPred)
        }else if (classTarget == hierarchy[2]){
          prediction_store <- list(classPred, latPred, longPred)
        }
      }
      return(list(model_store, prediction_store))
  }


## Pre-process of input file and mGPS implementation

### Pre-process of input file: sample size cut off
data_preprocess_f <- function(train_f,target_in,hierarchy,remove_small){
  
  data_in <- droplevels(train_f)
  
  if(length(hierarchy) == 4){
    data_in[,hierarchy[1]] <- factor(data_in[,hierarchy[1]])
    data_in[,hierarchy[2]] <- factor(data_in[,hierarchy[2]])
  }else{
    data_in[,hierarchy[1]] <- factor(data_in[,hierarchy[1]])
  }
  
  # Remove sparse samples locations and dubiously labeled samples. 
  remove_small <- as.numeric(remove_small)
  
  if (remove_small > 0){
    small_cities <-  names(which(summary(data_in[,target_in]) < remove_small))
    remove_samples <- which(data_in[,target_in] %in% small_cities)
    
    if (length(remove_samples) != 0 ){
      data_in <- droplevels(data_in[-c(remove_samples), ])
    }
  }
  
  # Remove samples that are mislabeling  
  for (i in 1:length(hierarchy)){
    empty <- which(data_in[,hierarchy[i]] == "" | is.na(data_in[,hierarchy[i]]))
    
    if (length(empty) != 0  ){
      data_in <- data_in[-c(empty),]
    }
  }
  
  data_in[,hierarchy[1]] <- droplevels(data_in[,hierarchy[1]])
  
  if(length(hierarchy) == 4){
    data_in[,hierarchy[2]] <- droplevels(data_in[,hierarchy[2]])
  }
  
  print("metadata has been pre-processed ...")
  return(data_in)
}


### Feature elimination of input data using feature selection algorithm
featureElimination_f <- function(metasub_data,classTarget_in,range_1,range_2,subsets_in){
  
  # Find GITs 
  range_1 <- as.numeric(range_1)
  range_2 <- as.numeric(range_2)
  
  featureElim <- species_select(x = metasub_data[, c(range_1:range_2)],
                                y = metasub_data[,classTarget_in],
                                remove_correlated = F, 
                                subsets = subsets_in,
                                cores = 8)
  return(featureElim)
}


### Output the optimal features found in feature elimination step
v_f <- function(featureElim){
  optVars <- featureElim$optVariables
  v <- varImp(featureElim$fit, type = 1, scale = F)
  v[,"taxa"] <- row.names(v)
  v <- v[order(v$Overall,decreasing = T),]
  write.csv(v, file = "Outputs/Optimal_features.csv")
  return(v)
}


### Output the accuracy of each feature subset size
git_subset_f <- function(featureElim){
  git_subset <- data.frame("n_vars" = featureElim$results$Variables, 
                           "accuracy" = featureElim$results$Accuracy)
  write.csv(git_subset,file = "Outputs/Features_subset_accuracy.csv")
  return(git_subset)
}


### mGPS implementation in Module one: Build a new prediction model using mGPS
model_accuracy_f <- function(metasub_data,optVars,classTarget_in,hierarchy_in){
  set.seed(18)
 
  # Generate 5 stratified folds for test predictions
  trainFolds <-  createFolds(metasub_data[,classTarget_in], k = 5, returnTrain = T)
  
  GeoPreds <- list()
  
  # Iteratively train the model on each of the 5 training folds and 
  # generate predictions using the corresponding test fold.
  for (i in 1:5){
    
    print(i)
    print("-------------------------------------------")
    
    train <- metasub_data[trainFolds[[i]],]
    test <- metasub_data[-trainFolds[[i]],]
    
    testPreds <-mGPS(training = train, testing = test, classTarget = classTarget_in,variables = optVars,nthread = 8,hierarchy = hierarchy_in, coast=NULL)

    GeoPreds[[i]] <- testPreds[[2]]

    print(i)
    print("has been done .....")
  }
  
  model_last_one <- testPreds[[1]]
  
  # Combine these test predictions into one data set 
  add_preds <- list()
  for (i in 1:5){
    add_preds[[i]] <- cbind(metasub_data[-trainFolds[[i]],] , 
                            "cityPred"= GeoPreds[[i]][[1]], 
                            "latPred" = GeoPreds[[i]][[2]], 
                            "longPred" = GeoPreds[[i]][[3]] )
  }
  
  MetasubDataPreds <- rbind.fill(add_preds)
  hierarchy <- hierarchy_in
  
  # calculate the Haversine distance between the predicted site and the sampling site
  for (i in 1:nrow(MetasubDataPreds)){
    
    if(length(hierarchy) == 3){
    
      MetasubDataPreds[i,"Distance_from_origin"] <-
      geosphere::distm(c(MetasubDataPreds[i,"longPred"],MetasubDataPreds[i,"latPred"]), 
                         c(MetasubDataPreds[i,hierarchy[3]],MetasubDataPreds[i,hierarchy[2]]), 
                         fun = geosphere::distHaversine)/1000
  
    }else  if(length(hierarchy) == 4){
      MetasubDataPreds[i,"Distance_from_origin"] <- 
      
geosphere::distm(c(MetasubDataPreds[i,"longPred"],MetasubDataPreds[i,"latPred"]), 
                         c(MetasubDataPreds[i,hierarchy[4]],MetasubDataPreds[i,hierarchy[3]]), 
                         fun = geosphere::distHaversine)/1000
    }
  }
  
  write.csv(MetasubDataPreds,"Outputs/Prediction_results.csv")
  save(model_last_one,file="Outputs/Prediction_model.Rda")
  
  return(list(MetasubDataPreds,model_last_one)) 
}


### mGPS implementation in Module two: 
#   Build a new prediction model using whole input data as training set
New_model_prediction <- function(train,test,classTarget_in,optVars,hierarchy_in){
  
  testPreds <-mGPS(training = train, testing = test, classTarget = classTarget_in,
                   variables = optVars,nthread = 8,hierarchy = hierarchy_in, coast=NULL)
  
  model <- testPreds[[1]]
  GeoPreds <- testPreds[[2]]
  
  MetasubDataPreds <- cbind(test, 
                            "cityPred"= GeoPreds[[1]], 
                            "latPred" = GeoPreds[[2]], 
                            "longPred" = GeoPreds[[3]] )
    
  write.csv(MetasubDataPreds,"Outputs/Prediction_result.csv")
  save(model,file="Outputs/Prediction_model.Rda")
  
  return(list(MetasubDataPreds,model))
}


### Function: Extract features that used in built prediction model when do the prediction
feature_filter <- function(test_dataset,model){
  features <- (model$finalModel$feature_names)
  newdata <- as.data.frame(setNames(replicate(length(features),numeric(0), simplify = F), features))
  for (f in features){
    if (f %in% colnames(test_dataset)){
      newdata[c(1:dim(test_dataset)[1]),f] <- test_dataset[,f]
    }else{
      newdata[c(1:dim(test_dataset)[1]),f] <- data.frame(0)
    }
  }
  return(newdata)
}


### mGPS implementation in Module three: Use constructed model to predict new samples
MetaSub_prediction <- function(model_store,test_dataset){
  Xgb_region <- model_store[[1]]
  Xgb_class <- model_store[[2]]
  Xgb_latitude <- model_store[[3]]
  Xgb_longitude <- model_store[[4]]
  # prediction_model <- model_store[[5]]
  
  print("---------- Prediction ----------")
  
  model <- function(test_dataset){
    
    newdata <- feature_filter(test_dataset, Xgb_region) 
    
    regProbs <- predict(Xgb_region, newdata, type ="prob")
    
    l1_test <- data.frame(newdata, regProbs)
    
    classPred <- predict(Xgb_class, newdata = l1_test)
    classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
    
    l2_test <-  data.frame(l1_test, classProbs) 
    latPred <- predict(Xgb_latitude, newdata = l2_test)
    
    l3_test <- data.frame(l2_test, latPred)
    longPred <- predict(Xgb_longitude, newdata = l3_test)
    
    longPred[longPred > 180] <- 180
    longPred[longPred < -180] <- -180
    latPred[latPred > 90] <- 90
    latPred[latPred < -90] <- -90
    
    return(list(classPred, latPred, longPred))
  }
  
  testPreds <- model(test_dataset)
  
  DataPreds <- cbind(test_dataset , 
                     "cityPred"= testPreds[[1]], 
                     "latPred" = testPreds[[2]], 
                     "longPred" = testPreds[[3]] )
  
  write.csv(DataPreds,"Outputs/Prediction_results.csv")
  return(DataPreds)
}


### Function: Pull predicted points to marine and land

# Pull to nearest coastline
pull_land <- function(land_preds,hierarchy){
  
  coastlines <- cbind("x"  = maps::SpatialLines2map(rworldmap::coastsCoarse)$x,
                      "y" =maps::SpatialLines2map(rworldmap::coastsCoarse)$y)
  coastlines <- coastlines[complete.cases(coastlines),]
  coastlines <- coastlines[coastlines[,1] < 180 ,]
  
  find_coast <- function(long, lat) {
    
    distances_from_coastline <-
      sp::spDistsN1(coastlines, c(long, lat), longlat = TRUE)
    
    closest_point <-  which.min(distances_from_coastline)
    new_coords <- coastlines[closest_point,]
    
    return(new_coords)
  }
  
  toAdjust <-
    which(is.na(maps::map.where(database = "world", land_preds$longPred, land_preds$latPred)))
  
  adjusted <-
    mapply(find_coast, long = land_preds$longPred[toAdjust], lat = land_preds$latPred[toAdjust])
  
  land_preds$longPred[toAdjust] <- adjusted[1,]
  land_preds$latPred[toAdjust] <- adjusted[2,]
  
  # calculate the Haversine distance between the predicted site and the sampling site
  for (i in 1:nrow(land_preds)){
    
    if(length(hierarchy) == 3){
      land_preds[i,"Distance_from_origin"] <- 
        geosphere::distm(c(land_preds[i,"longPred"],land_preds[i,"latPred"]), 
                         c(land_preds[i,hierarchy[3]],land_preds[i,hierarchy[2]]), 
                         fun = geosphere::distHaversine)/1000
      
    }else  if(length(hierarchy) == 4){
      land_preds[i,"Distance_from_origin"] <- 
        geosphere::distm(c(land_preds[i,"longPred"],land_preds[i,"latPred"]), 
                         c(land_preds[i,hierarchy[4]],land_preds[i,hierarchy[3]]), 
                         fun = geosphere::distHaversine)/1000
    }
  }
  return(land_preds)
}

# Pull to nearest ocean
pull_marine <- function(marine_preds,hierarchy){
  
  # Here need geography annotation files that provided in GitHub.
  seas <- rgdal::readOGR(dsn = "Data/Geo/ne_10m_geography_marine_polys", 
                         layer = "ne_10m_geography_marine_polys")
  coastlines <- cbind("x" =maps::SpatialPolygons2map(seas)$x,
                      "y" =maps::SpatialPolygons2map(seas)$y)
  coastlines <- coastlines[complete.cases(coastlines),]
  coastlines <- coastlines[coastlines[,1] < 180,]
  
  # Function: Find the closest site in ocean
  find_coast2 <- function(long,lat){
    distances_from_coastline <-  sp::spDistsN1(coastlines , c(long,lat), longlat = TRUE)
    closest_point <-  which.min(distances_from_coastline)
    new_coords <- coastlines[closest_point,]
    return(new_coords)
  }
  
  data(wrld_simpl)
  
  # Create a SpatialPoints object
  set.seed(0)
  points <- data.frame(marine_preds$longPred, marine_preds$latPred) 
  pts <- SpatialPoints(points, proj4string=CRS(proj4string(wrld_simpl)))
  # Change CRS if program has error for different CRS in over() step
  proj4string(wrld_simpl) <- pts@proj4string
  # Find which points fall over land
  ii <- !is.na(over(pts, wrld_simpl)$FIPS)
  toAdjust <- marine_preds[which(ii == TRUE),]
  adjusted <- mapply(find_coast2, long = toAdjust$longPred, lat = toAdjust$latPred )
  
  marine_preds[which(ii == TRUE), "latPred"] <- adjusted[2,]
  marine_preds[which(ii == TRUE), "longPred"] <- adjusted[1,]
  
  # calculate the Haversine distance between the predicted site and the sampling site
  for (i in 1:nrow(marine_preds)){
    
    if(length(hierarchy) == 3){
      marine_preds[i,"Distance_from_origin"] <- 
        geosphere::distm(c(marine_preds[i,"longPred"],marine_preds[i,"latPred"]), 
                         c(marine_preds[i,hierarchy[3]],marine_preds[i,hierarchy[2]]), 
                         fun = geosphere::distHaversine)/1000
      
    }else  if(length(hierarchy) == 4){
      marine_preds[i,"Distance_from_origin"] <- 
        geosphere::distm(c(marine_preds[i,"longPred"],marine_preds[i,"latPred"]), 
                         c(marine_preds[i,hierarchy[4]],marine_preds[i,hierarchy[3]]), 
                         fun = geosphere::distHaversine)/1000
    }
  }
  return(marine_preds)
}


## visualization of result

### Plot predicted site of samples in original dataset (Module one). 
plot_map <- function(MetasubDataPreds,hierarchy_in,classTarget_in,x_ran,y_ran){
  
  par(mai=c(2,1,0.5,0.5), mar=par()$mar+c(3,0,0,0))
  
  map <- rworldmap::getMap(resolution = "coarse")
  palette <-c( "darkorchid4","gold2","dodgerblue3","brown","orangered2",
               "mediumspringgreen","deeppink2")
  plot(map, xlim = x_ran,ylim = y_ran,col = "grey",border = "darkgrey", 
       xlab = "", ylab = '', bg = "lightskyblue1")
  title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
  
  for (i in 1:(length(hierarchy_in)-2)){
    MetasubDataPreds[,hierarchy_in[i]] <- factor(MetasubDataPreds[,hierarchy_in[i]])
  }
  
  MetasubDataPreds$cityPred <- factor(MetasubDataPreds$cityPred)
  
  # Get the coordination of continent
  continent_list <- c("east_asia","europe","middle_east","north_america","oceania",
                      "south_america", "sub_saharan_africa",'Australia','North_America',
                      'Europe','Africa','South_america','Asia')
  continent_lats <- c(55,69,8,40,-40,-10,-5,-40,40,69,-5,-10,55)
  continent_longs <- c(125,0,60,-130,140,-80,5,140,-130,0,5,-80,125)
  continent_position <- data.frame(cbind(continent_list,continent_lats,continent_longs))
  
  label_continent <- c()
  flag <- 0
  for ( i in 1:length(levels(MetasubDataPreds[,hierarchy_in[1]]))){
    this_continent <- levels(MetasubDataPreds[,hierarchy_in[1]])[i]
    label_continent <- c(label_continent,this_continent)
    find_lats <- MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,][,"latPred"]
    find_longs <- MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,][,"longPred"]
    
    # Plot predicted site with coordinates
    points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.2,xlim = c(-165,168))
    
    # Plot regional prediction accuracy by continent as pies
    correctly_pred <-  mean(as.numeric(as.character(MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,"cityPred"])== 
                                                      as.character(MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,classTarget_in]))) 
    incorrectly_pred <- (1 - correctly_pred)
    
    if (this_continent %in% continent_position$continent_list){
      add.pie(z = c(correctly_pred, incorrectly_pred), 
              x = as.numeric(continent_position[continent_position$continent_list==this_continent,3]), 
              y = as.numeric(continent_position[continent_list==this_continent,2]),
              edges=200,
              radius=10,
              col=c(palette[i],"black"), labels = "")
    }else{
      add.pie(z = c(correctly_pred, incorrectly_pred), x = as.numeric(-150+flag), y = 105
              ,edges=200,
              radius=10,
              col=c(palette[i],"black") , labels = ""
      )
      flag <- flag + 25
    }
  }
  map.axes(cex.axis = 1.1)
  legend(xpd = T,'bottom',inset=c(0,-0.35), label_continent, pch = "+", 
         col = palette[1:length(label_continent)], cex = 0.8,n=4)
  box( col = 'black')
  
  par(mar=c(5, 4, 4, 2) + 0.1)
}


### Plot the prediction source site of new sample on world map (Module two and three)
plot_prediction <- function(MetasubDataPreds,x_ran,y_ran){
  
  par(mai=c(1,1,0.5,0.5))

  map <- rworldmap::getMap(resolution = "coarse")
  plot(map, xlim = x_ran,ylim = y_ran, col = "grey",border = "darkgrey", 
       xlab = "", ylab = '', bg = "lightskyblue1")
  title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)

  find_lats <- MetasubDataPreds$latPred
  find_longs <- MetasubDataPreds$longPred
  
  points(find_longs,find_lats, type = "p",col = "purple", pch = "+", cex = 1.3)
  map.axes(cex.axis = 1.1)
}


### Plot accuracy bar of constructed model (Module one)
plot_accuracy <- function(MetasubDataPreds,classTarget_in){
  MetasubDataPreds[,classTarget_in] <- factor(MetasubDataPreds[,classTarget_in])
  
  bar_df1 <- data.frame(row.names = c("Overall",levels(MetasubDataPreds[,classTarget_in])))
  
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 100)
    bar_df1[i+1,"0 - 100km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] < 100)
    bar_df1[ 1,"0 - 100km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 100 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 500)
    bar_df1[i+1,"100 - 500km"] <- prop
    
    overall_prop <-mean(MetasubDataPreds[,"Distance_from_origin"] > 100 & MetasubDataPreds[,"Distance_from_origin"] < 500)
    bar_df1[ 1,"100 - 500km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 500 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 1000)
    bar_df1[i+1,"500 - 1000km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 500 & MetasubDataPreds[,"Distance_from_origin"] < 1000)
    bar_df1[ 1,"500 - 1000km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 1000 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 2000)
    bar_df1[i+1,"1000 - 2000km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 1000 & MetasubDataPreds[,"Distance_from_origin"] < 2000)
    bar_df1[ 1,"1000 - 2000km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 2000 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 3000)
    bar_df1[i+1,"2000 - 3000km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 2000 & MetasubDataPreds[,"Distance_from_origin"] < 3000)
    bar_df1[1,"2000 - 3000km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 3000 )
    bar_df1[i+1,"> 3000km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 3000)
    bar_df1[ 1,"> 3000km"] <- overall_prop
  }
  
  # Find number of samples per region
  size <- c()
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    size[i] <- length(which(MetasubDataPreds[,classTarget_in] == this_city))
  }
  
  par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
  bp <- barplot(t(bar_df1*100), space = 0,col=c("lightyellow","slategray1","lightblue", 
                                                "skyblue", "royalblue3", "darkblue"), 
                names.arg=c("Overall",paste0(levels(MetasubDataPreds[,classTarget_in]),"  (",size,")"), 
                            axes = FALSE) , 
                las =2, cex.names=.6, ylab = "", axisnames = F, axes = F)
  axis(side =2, pos = 0)
  mtext(text = c("Overall",paste0(levels(MetasubDataPreds[,classTarget_in])," (",size,")")), 
        side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
  title(ylab="Proportion of sample predictions %", mgp=c(2,1,0),cex.lab=1)
  legend("topright",inset = c(-0.26,0.4), rev(c(colnames(bar_df1))), 
         fill = rev(c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")),
         bty = 1, cex = 0.8)
  par(mar=c(5, 4, 4, 2) + 0.1)
}


# Part 2 Interface constructed by shiny ------------------------------------------------------------------
  
# ui

library(shiny)

ui <- fluidPage(
  titlePanel('Geographical origin prediction of microbiome'),
  
  tags$head(
    tags$style(HTML("
    .shiny-output-error-validation {
    color: red;
    }
    "))
  ),
  
  # side bar panel part
  sidebarLayout(
    sidebarPanel(
      
      # Data loading tips
      tags$head(tags$style(type="text/css", "
                #loadmessage {
                top: 0px; left: 0px;
                width: 100%; padding: 5px 0px 5px 0px;
                text-align: center; font-weight: bold;
                font-size: 100%; color: #000000;
                background-color: #FFC1C1; z-index: 105;}"),
                ), 
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div("Data loading....",id="loadmessage")),
      
      # Program selection
      selectInput("program","Prediction program",
                  c("Build a new prediction model using mGPS" = "cc",
                    "Build a new prediction model using mGPS and predict new samples" = "bb",
                    "Use an existing model to predict new samples" = "aa"
                    )
                  ),
      
      # Module three : Microbial source prediction based on an exist model
      conditionalPanel(
        condition = "input.program == 'aa'",
        fileInput(inputId = "f_new_test_1",
                  label = "Upload sample(s) abundance file",
                  accept = c("text/csv",
                             ".csv")),
        fileInput(inputId = "model",
                  label = "Upload the prediction model (In .Rda format)",
                  accept = ".Rda"),
        actionButton(inputId = "acbutton_1", 
                     label ="Start")
      ),
      
      # Module two : Build a new prediction model using mGPS and predict new samples
      conditionalPanel(
        condition = "input.program == 'bb'",
        fileInput(inputId = "f_new_test_2",
                  label = "Upload new sample abundance file(s)",
                  accept = c("text/csv",
                             ".csv")),
        radioButtons(inputId = "file_num_2",
                     label = "Upload reference file(s)",
                     c("Merged metadata and abundance file" = "one_file",
                       "Separate metadata and abundance file" = "two_file")
        ),
        conditionalPanel(
          condition = "input.file_num_2 == 'one_file'",
          fileInput(inputId = "f_new_train_2",
                    label = "Upload the reference merged dataset file",
                    accept = c("text/csv",
                               ".csv")),
        ),
        conditionalPanel(
          condition = "input.file_num_2 == 'two_file'",
          fileInput(inputId = "f_metadata_2",
                    label = "Upload the metadata file",
                    accept = c("text/csv",
                               ".csv")),
          fileInput(inputId = "f_abundance_2",
                    label = "Upload the abundance file",
                    accept = c("text/csv",
                               ".csv")),
          textInput(inputId = "by_x_2",
                    label = "Merge column name in metadata file"),
          textInput(inputId = "by_y_2",
                    label = "Merge column name in abundance file")
        ),
        textInput(inputId= "target_2",
                  label ="Enter the main locality level"),
        textAreaInput(inputId = "hierarchy_2",
                      label = "Enter the locality hierarchy",
                      rows = 1),
        textAreaInput(inputId = "abundance_range_2",
                      label = "Column range of abundance data",
                      rows = 1),
        
        checkboxInput('remove_small_2',"Locality sample size cut off (Optional)", value = FALSE),
        conditionalPanel( 
          condition = "input.remove_small_2 == true",
          textAreaInput(inputId = "remove_small_value_2",
                        label = "Cut off of sample number",
                        rows = 1),
        ),
        
        checkboxInput("subsets_2", "Subsets in feature elimination (Optional)", value = FALSE),
        conditionalPanel( 
          condition = "input.subsets_2 == true",
          textAreaInput(inputId = "subsets_value_2",
                        label = "Subsets size",
                        rows = 1),
        ),
        
        actionButton(inputId = "acbutton_2", 
                     label ="Start")
      ),
      
      # Module one : Build a new prediction model using mGPS and test accuracy
      conditionalPanel(
        condition = "input.program == 'cc'",
        radioButtons(inputId = "file_num_3",
                     label = "Input file",
                     c("Merged metadata and abundance file" = "one_file",
                       "Separate metadata and abundance file" = "two_file")
                      ),
        conditionalPanel(
          condition = "input.file_num_3 == 'one_file'",
          fileInput(inputId = "f_new_train_3",
                    label = "Upload the reference merged dataset file",
                    accept = c("text/csv",
                               ".csv")),
        ),
        conditionalPanel(
          condition = "input.file_num_3 == 'two_file'",
          fileInput(inputId = "f_metadata_3",
                    label = "Upload the metadata file",
                    accept = c("text/csv",
                               ".csv")),
          fileInput(inputId = "f_abundance_3",
                    label = "Upload the abundance file",
                    accept = c("text/csv",
                               ".csv")),
          textInput(inputId = "by_x_3",
                    label = "Merge column name in metadata file"),
          textInput(inputId = "by_y_3",
                    label = "Merge column name in abundance file")
        ),
        textInput(inputId= "target_3",
                  label ="Enter the main locality level"),
        textAreaInput(inputId = "hierarchy_3",
                      label = "Enter the locality hierarchy",
                      rows = 1),
        textAreaInput(inputId = "abundance_range_3",
                      label = "Column range of abundance data",
                      rows = 1),
        
        checkboxInput('remove_small_3',"Locality sample size cut off (Optional)", value = FALSE),
        conditionalPanel( 
          condition = "input.remove_small_3 == true",
          textAreaInput(inputId = "remove_small_value_3",
                        label = "Cut off of sample number",
                        rows = 1),
        ),
        
        checkboxInput("subsets_3", "Subsets in feature elimination (Optional)", value = FALSE),
        conditionalPanel( 
          condition = "input.subsets_3 == true",
          textAreaInput(inputId = "subsets_value_3",
                        label = "Subsets size",
                        rows = 1),
        ),
        
        actionButton(inputId = "acbutton_3", 
                     label ="Start")
     )
  ),
    
  # main panel part
    mainPanel(
      
      conditionalPanel(
        condition = "input.program == 'aa'",
        tabsetPanel(
          type = "tabs",
          tabPanel("HELP",
                   h3("Use an existing model to predict new samples"),
                   h4("_______________________________________________"),
                   br(),
                   h4("Function description"),
                   br(),
                   p("In this mode, user can predict new sample origin based on an exsiting prediction model."),
                   h4("_______________________________________________"),
                   br(),
                   h4("Usage"),
                   br(),
                   p("In left side bar:",style = "font-size:16px"),
                   p("1. Select ",strong("Prediction program")," as ",em("Use an existing model to predict new samples", style = "color:purple")), 
                   p("2. ",strong("Upload sample(s) abundance file:")," Upload data file(s) (in .csv format) containing new microbial sample abundance data."),  
                   p("3. ",strong("Upload the prediction model:")," Upload a constructed origin prediction model in .Rda format. Model can be downloaded in ", strong("Output")," tab of function:",em("Build a new prediction model using mGPS", style = "color:purple"), "or ",em("Build a new prediction model using mGPS and predict new samples", style = "color:purple")),
                   br(),
                   p("In ",strong("Result Plot")," tab:",style = "font-size:16px"),
                   p("4. ",strong("Change longitude/latitude range in output map"),em(" (Optional) ")),
                   p("5. ",strong("Whether pull points to land/marine:"),em(" (Optional) "),"If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder."),
                   p(strong("Start the program:")," Click the",strong("Start"),"bar and then click the ",strong("Result Plot")," tab",style = "font-size:16px"),
                   p('Data processing, please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.',style = "font-size:16px"),
                   h4("_______________________________________________"),
                   br(),
                   h4("Result plot and output"),
                   br(),
                   p(strong("Result Plot")," tab:"),
                   p("The exsiting model will be used to predict the origin of new sample."),
                   p("* World map: new samples' prediction origins are plotted on the world map"),
                   br(),
                   p(strong("Output")," tab:"),
                   p("The results of using the exsiting prediction model to predict the source coordinates of new sample(s)."),
                   p("For more introduction of output results, view the tutorial file in ",
                     a("mGPS interface Gitbub",href = "https://github.com/YaliZhang98/mGPS_interface"))
                   ,h4("_______________________________________________"),
                   ),
          tabPanel("Result Plot",
                   fluidRow( column(6, sliderInput("xlim_1", "longitude range:",
                                                   min = -165, max = 168,
                                                   value = c(-165,168))), 
                             column(6, sliderInput("ylim_1", "latitude range:",
                                                   min = -90, max = 90,
                                                   value = c(-90,90)))) ,
                   radioButtons("pull_1", "Whether pull points to land/marine",
                                choices = c("Pull to land" = "land_1",
                                            "Pull to waterbody" = "marine_1",
                                            "Default" = "none"),
                                selected = "none"), 
                   plotOutput(outputId = "predicted_map_1"),
                   downloadButton("downloadMap_1",
                                  "DownloadMap")
                   ),
          
          tabPanel("Output",
                   helpText("Here you can download the predicted original coordinates of samples in your uploaded file. 
                            (Reference microbiome locations are based on MetaSub project microbiome database)"),
                   
                   downloadButton("downloadData_1",
                                  "DownloadData")
                   ),
        )
      ),

      conditionalPanel(
        condition = "input.program == 'bb'",
        tabsetPanel(
          type = "tabs",
          tabPanel("HELP",
                   h3("Build a new prediction model using mGPS and predict new samples"),
                   h4("_______________________________________________"),
                   br(),
                   h4("Function description"),
                   br(),
                   p("In this mode, user can train the microbial origin prediction model based on the reference data set uploaded by the user. The constructed prediction model will be used to predict the new sample to be tested provided by the user and report the prediction result of the sample source. (If user want to visualize the accuracy of the model, please use function:",em("Build a new prediction model using mGPS", style = "color:purple"),")"),
                   h4("_______________________________________________"),
                   br(),
                   h4("Usage"),
                   br(),
                   p("In left side bar:",style = "font-size:16px"),
                   p("1. Select ",strong("Prediction program")," as ",em("Build a new prediction model using mGPS and predict new samples", style = "color:purple")), 
                   p("2. ",strong("Upload new sample(s) abundance file:")," Upload file (in .csv format) containing abundance data of new sample(s)."),  
                   
                   p("3. ",strong("Upload reference file(s):")," Upload data file(s) (in .csv format) containing microbial abundance data and metadata."),  
                   p("   In metadata, at least one locality (eg. continent, city) and coordinates (necessary) data columns should be included. The metadata and abundance data of the sample can be merged into one file (",em("Merged metadata and abundance data", style = "color:purple"),"), or uploaded as two files (",em("Separate metadata and abundance data", style = "color:purple"),")"),
                   p("When ",em("Separate metadata and abundance file", style = "color:purple")," is selected. ",strong("Merge column name in metadata/abundance file: "),"Input the header name of column which is the merged column in two files."),
                   p("4. ",strong("Enter the main locality level:")," Input the main locality target. It should same as that column header. (eg. city)"),
                   p("5. ",strong("Enter the locality hierarchy:")," The locality chain used in mGPS to construct the prediction model (same column headers). It should contain one or two locality information, latitude and longitude. Use ',' as the separator. (eg. continent,city,latitude,longitude)"),
                   p("6. ",strong("Columns range of abundance data:")," Input the columns range number of abundance data in the abundance/merged file. Use ':' as separator (eg 44:1000)"),
                   p("7. ",strong("Locality sample size cut off:"),em(" (Optional) "),"Remove locality whose sample size is less than a certain value. If checked, input the cut off number (eg. 8)"),
                   p("8. ",strong("Subsets in feature elimination:"),em(" (Optional) ")," Limit the number of features to a certain value. If unchecked, mGPS will find the optimal subset size of microbiome features. If checked, there are three types of input format: a.input the subsets size with separator as ',' (eg. 50,100,200,300); b. input the subsets size range with separator as '-' (eg. 50-300); c. input a single value."),
                   br(),
                   p("In ",strong("Result Plot")," tab:",style = "font-size:16px"),
                   p("9. ",strong("Change longitude/latitude range in output map"),em(" (Optional) ")),
                   p("10. ",strong("Whether pull points to land/marine:"),em(" (Optional) "),"If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder."),
                   p(strong("Start the program:")," Click the",strong("Start"),"bar and then click the ",strong("Result Plot")," tab",style = "font-size:16px"),
                   p('Data processing, please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.',style = "font-size:16px"),
                   h4("_______________________________________________"),
                   br(),
                   h4("Result plot and output"),
                   br(),
                   p(strong("Result Plot")," tab:"),
                   p("The reference datasets will be used to construct a origin prediction model by mGPS.
                     Then this model will be used to predict origin of new samples. The prediction origin map and
                     prediction accuracy plot can be downloaded by bottons under the plots."),
                   p("* World map: samples' prediction origins are plotted on the world map"),
                   br(),
                   p(strong("Output")," tab:"),
                   p("The results of using the constructed prediction model to predict the source coordinates of the new samples. In addition, the constructed prediction model can be downloaded (In Rda format)."),
                   p("For more introduction of output results, view the tutorial file in ",
                     a("mGPS interface Gitbub",href = "https://github.com/YaliZhang98/mGPS_interface")),
                   h4("_______________________________________________"),
                   ),

          tabPanel("Result Plot",
                   fluidRow( column(6, sliderInput("xlim_2", "longitude range:",
                                                   min = -165, max = 168,
                                                   value = c(-165,168))), 
                             column(6, sliderInput("ylim_2", "latitude range:",
                                                   min = -90, max = 90,
                                                   value = c(-90,90)))) ,
                   radioButtons("pull_2", "Whether pull points to land/marine",
                                choices = c("Pull to land" = "land_2",
                                            "Pull to waterbody" = "marine_2",
                                            "Default" = "none_2"),
                                selected = "none_2"),   
                   plotOutput(outputId = "predicted_map_2"),
                   downloadButton("downloadMap_2",
                                  "DownloadMap")),
          
          tabPanel("Output",
                   helpText("Predicted original coordinates of samples can be downloaded. Prediction model in RDA format can be downloaded"),
                   downloadButton("downloadData_2",
                                  "Download prediction data"),
                   downloadButton("download_featuresub_2",
                                  "Download feature subsets accuracy in feature elimination"),
                   downloadButton("download_optfeatures_2",
                                  "Download optimal features in prediction model"),
                   downloadButton("downloadModel_2",
                                  "DownloadModel")
                   )
        )
      ),
      
      conditionalPanel(
        condition = "input.program == 'cc'",
        tabsetPanel(
          type = "tabs",
          
          tabPanel("Welcome",
                   h2("Welcome to mGPS application"),
                   br(),
                   p("This is a web program based on the mGPS application created by 
                      Shiny. It can build a microbial origin prediction model and predict 
                      the origin of microbes.", style = "font-size:16px"),
                   p("To learn more about mGPS, please visit:",
                     a(em("mGPS"),
                       href="https://github.com/eelhaik/mGPS"), style = "font-size:16px"),
                   br(),
                   h3("Function"),
                   br(),
                   p("This program contains three functions. These function can be performed by selected ",strong("Prediction program")," in the left side bar. The detail usage 
                      will be introduced in ",strong("HELP")," tab for each function.", style = "font-size:16px"),
                   br(),
                   p("1. ",strong("Build a new prediction model using mGPS")),
                   p("In this mode, user can use the mGPS tool to build a microbial source prediction model based
                            on the microbial abundance data and coordinates data uploaded by the user.","To learn more about mGPS, please visit:",
                     a(em("mGPS"),
                       href="https://github.com/eelhaik/mGPS"),
                     p("2. ",strong("Build a new prediction model using mGPS and predict new samples")),
                     p("In this mode, user can train the microbial origin
                        prediction model based on the reference data set uploaded by the user. The
                        constructed prediction model will be used to predict the new sample to be tested provided
                        by the user and report the prediction result of the sample source. (If user want
                        to visualize the accuracy of the model, please use function:",em("Build a new prediction model using mGPS"),")"),
                     p("3. ",strong("Use an existing model to predict new samples")),
                     p("In this mode, user can predict new sample origin based on an existing prediction model. Model can be downloaded in ", strong("Output")," tab of function:",em("Build a new prediction model using mGPS", style = "color:purple"), "or ",em("Build a new prediction model using mGPS and predict new samples", style = "color:purple")),
                     br(),
                     p("For more detail introduction and examples, visit the ",
                       a("mGPS interface on Gitbub", 
                         href = "https://github.com/YaliZhang98/mGPS_interface"), style = "font-size:16px"),
                   )),
          
          tabPanel("HELP",
                   h3("Build a new prediction model using mGPS"),
                   h4("_______________________________________________"),
                   br(),
                   h4("Function description"),
                   br(),
                   p("This model can use the mGPS tool to build a microbial source prediction model based
                            on the microbial abundance data and coordinates data uploaded by the user."),
                   h4("_______________________________________________"),
                   br(),
                   h4("Usage"),
                   br(),
                   p("In left side bar:",style = "font-size:16px"),
                   p("1. Select ",strong("Prediction program")," as ",em("Build a new prediction model using mGPS", style = "color:purple")), 
                   p("2. ",strong("Input file(s):")," Upload data file(s) (in .csv format) containing microbial abundance data and metadata."),  
                   p("   In metadata, at least one locality (eg. continent, city) and coordinates (necessary) data columns should be included. The metadata and abundance data of the sample can be merged into one file (",em("Merged metadata and abundance data"),"), or uploaded as two files (",em("Separate metadata and abundance data"),")"),
                   p("When ",em("Separate metadata and abundance file")," is selected. ",strong("Merge column name in metadata/abundance file: "),"Input the header name of column which is the merged column in two files."),
                   p("3. ",strong("Enter the main locality level:")," Input the main locality target. It should same as that column header. (eg. city)"),
                   p("4. ",strong("Enter the locality hierarchy:")," The locality chain used in mGPS to construct the prediction model (same column headers). It should contain one or two locality information, latitude and longitude. Use ',' as the separator. (eg. continent,city,latitude,longitude)"),
                   p("5. ",strong("Columns range of abundance data:")," Input the columns range number of abundance data in the abundance/merged file. Use ':' as separator (eg 44:1000)"),
                   p("6. ",strong("Locality sample size cut off:"),em(" (Optional) "),"Remove locality whose sample size is less than a certain value. If checked, input the cut off number (eg. 8)"),
                   p("7. ",strong("Subsets in feature elimination:"),em(" (Optional) ")," Limit the number of features to a certain value. If unchecked, mGPS will find the optimal subset size of microbiome features. If checked, there are three types of input format: a.input the subsets size with separator as ',' (eg. 50,100,200,300); b. input the subsets size range with separator as '-' (eg. 50-300); c. input a single value."),
                   br(),
                   p("In ",strong("Result Plot")," tab:",style = "font-size:16px"),
                   p("8. ",strong("Change longitude/latitude range in output map"),em(" (Optional) ")),
                   p("9. ",strong("Whether pull points to land/marine:"),em(" (Optional) "),"If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder."),
                   p(strong("Start the program:")," Click the",strong("Start"),"bar and then click the ",strong("Result Plot")," tab",style = "font-size:16px"),
                   p('Data processing, please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.',style = "font-size:16px"),
                   h4("_______________________________________________"),
                   br(),
                   h4("Result plot and output"),
                   br(),
                   p(strong("Result Plot")," tab:"),
                   p("Show the accuracy
                      of the prediction model trained by the mGPS tool and based on the reference microbial database
                      uploaded by the user."),
                   p("The original database will be divided into 5 folds, and mGPS will use 4 of
                      these folds to train the model, and the resulting model will be used to predict the microbial source
                      of the remaining fold. Iteratively obtain the prediction result of the original database and compare
                      it with the actual location of the microorganism."),
                  
                   p("* World map: samples' prediction origins are plotted on the world map"),
                   p("* Accuracy bar plot: model built by mGPS accuracy per site for the original reference dataset. mGPS accuracy is shown per-site as the distances between the predicted and true sampling site for the reference samples. The average prediction accuracy across all samples with each population given equal weight is shown on the left."),
                   br(),
                   p(strong("Output")," tab:"),
                   p("The results of using the constructed prediction model to predict the source coordinates of the original dataset samples. In addition, the constructed prediction model can be downloaded (In Rda format)."),
                   p("For more introduction of output results, view the tutorial file in ",
                    a("mGPS interface Gitbub",href = "https://github.com/YaliZhang98/mGPS_interface")),
                   h4("_______________________________________________"),
          ),
          
          tabPanel("Result Plot",
                   
                   fluidRow( column(6, sliderInput("xlim_3", "longitude range:",
                                                   min = -165, max = 168,
                                                   value = c(-165,168))), 
                             column(6, sliderInput("ylim_3", "latitude range:",
                                                   min = -90, max = 120,
                                                   value = c(-90,120)))),
                   
                   radioButtons("pull_3", "Whether pull points to land/marine",
                                choices = c("Pull to land" = "land_3",
                                            "Pull to waterbody" = "marine_3",
                                            "Default" = "none"),
                                selected = "none"),   
                   plotOutput(outputId = "predicted_map_3"),
                   downloadButton("downloadMap_3",
                                  "DownloadMap"),
                   plotOutput(outputId = "predicted_accuracy_3"),
                   downloadButton("downloadAccuracy_3",
                                  "DownloadPlot")
                   ),
          tabPanel("Output",
                   helpText("Predicted original coordinates of samples can be downloaded. Prediction model in RDA format can be downloaded"),
                   downloadButton("downloadData_3",
                                  "Download prediction data"),
                   downloadButton("download_featuresub_3",
                                  "Download feature subsets accuracy in feature elimination"),
                   downloadButton("download_optfeatures_3",
                                  "Download optimal features in prediction model"),
                   downloadButton("downloadModel",
                                  "DownloadModel"))
      )
    )
  )
))


# -------------------------------------------------------------------

# Server

options(shiny.maxRequestSize=200*1024^2)

server <- function(input,output){
    
# Function used in three modules: Pull the predicted points to land or marine
    pull_MetasubDataPreds <- reactive({
      hierarchy_here <- hierarchy_in()
      if (input$pull_1 == "land_1"){
        pull_MetasubDataPreds <- pull_land(MetasubDataPreds(),hierarchy_here)
      }else if (input$pull_1 == "marine_1"){
        pull_MetasubDataPreds <- pull_marine(MetasubDataPreds(),hierarchy_here)  
      }else{
        pull_MetasubDataPreds <- MetasubDataPreds()
      }
      return(pull_MetasubDataPreds)
    })
  
# Module three implementation
    
    # Get the upload model in Module three 
    model_store <- reactive({
      req(input$model)
      path <- input$model$datapath
      validate(
         need('Rda' %in%  strsplit(input$model$datapath,".",fixed = T)[[1]], 
              "The model file you upload should be in .Rda format.")
      )
      load(path)
        return(model_store)
    })
      
    # Get the upload data file
    test_f <- reactive({
      if(input$program == "aa" ){
        req(input$f_new_test_1)
        test_f <- read.csv(input$f_new_test_1$datapath, 
                            header = TRUE)
        print(input$f_new_test_1$datapath)
        validate(
          need('csv' %in%  strsplit(input$f_new_test_1$datapath,".",fixed = T)[[1]], 
                "The sample abundance file you upload should be in.csv format.")
        )
      }
      if(input$program == "bb" ){
        req(input$f_new_test_2)
        test_f <- read.csv(input$f_new_test_2$datapath, 
                           header = TRUE)
        validate(
            need('csv' %in%  strsplit(input$f_new_test_2$datapath,".",fixed = T)[[1]], 
                "The sample abundance file you upload should be in.csv format.")
        )}
        return(test_f)
    })
    
    # Module three: use an exsit model to predict new sample
    MetasubDataPreds <- reactive({
      MetasubDataPreds <- MetaSub_prediction(model_store(),test_f())
      return(MetasubDataPreds)    
    })

# Module one and Module two implementation
    
   # Get header of merge columns in RSA file and metadata file 
    by_x_in <- reactive({
      if(input$program == "bb" ){
        req(input$by_x_2)
        by_x_in <- input$by_x_2
      }
      if(input$program == "cc" ){
        req(input$by_x_3)
        by_x_in <- input$by_x_3
      }
      return(by_x_in)
    })
    
    by_y_in <- reactive({
      if(input$program == "bb" ){
        req(input$by_y_2)
        by_y_in <- input$by_y_2
      }
      if(input$program == "cc" ){
        req(input$by_y_3)
        by_y_in <- input$by_y_3}
      return(by_y_in)
    })
    
    # Get class target 
    classTarget_in <- reactive({
      if(input$program == "bb" ){
        req(input$target_2)
        classTarget_in <- input$target_2}
      if(input$program == "cc" ){
        req(input$target_3)
        classTarget_in <- input$target_3}
      return(classTarget_in)
    })
    
    # Get hierarchy level
    hierarchy_in <- reactive({
      if(input$program == "bb" ){
        req(input$hierarchy_2)
        text <- input$hierarchy_2
        hierarchy_in <- strsplit(text,",")[[1]]}
      if(input$program == "cc" ){
        req(input$hierarchy_3)
        text <- input$hierarchy_3
        hierarchy_in <- strsplit(text,",")[[1]]}
      return(hierarchy_in)
    })
    
    # Get subset size in feature elimination
    subsets_in <- reactive({
      if(input$program == "bb" ){
        if (input$subsets_2 == T){
          req(input$subsets_value_2)
          text <- input$subsets_value_2
          
          if ('-' %in% strsplit(text,"")[[1]]){
            subsets_range <- as.numeric(strsplit(text,"-")[[1]])
            
            subsets_in <- c(subsets_range[1],
                            round((subsets_range[2]-subsets_range[1])/5+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*2+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*3+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*4+subsets_range[1]),
                            subsets_range[2])
          }else{
            subsets_in <- as.numeric(strsplit(text,",")[[1]])
          }
          
         }else{
          subsets_in <- NULL
          }}
      if(input$program == "cc" ){
        if (input$subsets_3 == T){
          req(input$subsets_value_3)
          text <- input$subsets_value_3
          if ('-' %in% strsplit(text,"")[[1]]){
            subsets_range <- as.numeric(strsplit(text,"-")[[1]])
            
            subsets_in <- c(subsets_range[1],
                            round((subsets_range[2]-subsets_range[1])/5+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*2+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*3+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*4+subsets_range[1]),
                            subsets_range[2])
          }else{
            subsets_in <- as.numeric(strsplit(text,",")[[1]])
          }
        }else{
          subsets_in <- NULL
        }}
      return(subsets_in)
    })
    
    # Get sample size cut off
    remove_small <- reactive({
      if(input$remove_small_2 == T ){
        req(input$remove_small_value_2)
        remove_small <- input$remove_small_value_2
        remove_small <- as.numeric(remove_small)
      }else{
        remove_small <- 0
      }
      
      if(input$remove_small_3 == T ){
        req(input$remove_small_value_3)
        remove_small <- input$remove_small_value_3
        remove_small <- as.numeric(remove_small)
      }else{
        remove_small <- 0
      }
      return(remove_small)
    })
    
    # Get column position of RSA data
    abundance_r <- reactive({
      if(input$program == "bb" ){
        req(input$abundance_range_2)
        text <- input$abundance_range_2
        range_1 <- as.numeric(strsplit(text,":")[[1]][1])
        range_2 <- as.numeric(strsplit(text,":")[[1]][2])
        abundance_r <- list(range_1,range_2)}
      if(input$program == "cc" ){
        req(input$abundance_range_3)
        text <- input$abundance_range_3
        range_1 <- as.numeric(strsplit(text,":")[[1]][1])
        range_2 <- as.numeric(strsplit(text,":")[[1]][2])
        abundance_r <- list(range_1,range_2)}
      return(abundance_r)
    })
    
    # Get upload file in Module one and two
    train_f <- reactive({
        
      if(input$program == "bb" ){
          
        if(input$file_num_2 == "one_file" ){
            req(input$f_new_train_2)
            train_f <- read.csv(input$f_new_train_2$datapath, 
                                header = TRUE)
            validate(
              need('csv' %in%  strsplit(input$f_new_train_2$datapath,".",fixed = T)[[1]], 
                   "The training file you upload should be in.csv format.")
          )
            return(train_f)
        }
          
        if (input$file_num_2 == "two_file"){
            req(input$f_metadata_2)
            req(input$f_abundance_2)
            data_metadata <- read.csv(input$f_metadata_2$datapath, 
                                      header = TRUE)
            validate(
              need('csv' %in%  strsplit(input$f_metadata_2$datapath,".",fixed = T)[[1]], 
                   "The metadata file you upload should be in.csv format.")
            )
            data_abundance <- read.csv(input$f_abundance_2$datapath,
                                       header = T)
            validate(
              need('csv' %in%  strsplit(input$f_abundance_2$datapath,".",fixed = T)[[1]], 
                   "The abundance file you upload should be in.csv format.")
            )
            data_abundance <- unique(data_abundance)
            by_x <- by_x_in()
            by_y <- by_y_in()
            
            metasub_data <- merge(data_metadata,data_abundance,by.x=by_x,by.y=by_y) # merge bacterial and meta data
            train_f <- metasub_data
            return(train_f)
        }
      }
        
      if(input$program == "cc" ){
        
        if(input$file_num_3 == "one_file" ){
          
          req(input$f_new_train_3)
          train_f <- read.csv(input$f_new_train_3$datapath, 
                              header = TRUE)
          validate(
            need('csv' %in%  strsplit(input$f_new_train_3$datapath,".",fixed = T)[[1]], "The training file you upload should be in.csv format.")
          )
          return(train_f)
        }
        
        if (input$file_num_3 == "two_file"){
          
          req(input$f_metadata_3)
          req(input$f_abundance_3)
          data_metadata <- read.csv(input$f_metadata_3$datapath, 
                                    header = TRUE)
          validate(
            need('csv' %in%  strsplit(input$f_metadata_3$datapath,".",fixed = T)[[1]], 
                 "The metadata file you upload should be in.csv format.")
          )
          data_abundance <- read.csv(input$f_abundance_3$datapath,
                                     header = T)
          validate(
            need('csv' %in%  strsplit(input$f_abundance_3$datapath,".",fixed = T)[[1]], 
                 "The abundance file you upload should be in.csv format.")
          )
          data_abundance <- unique(data_abundance)
          metasub_data <- merge(data_metadata,data_abundance,by.x=by_x_in(),by.y=by_y_in()) # merge bacterial and meta data
          train_f <- metasub_data
          return(train_f)
        }
      }
    })
    
    # data pre-process
    metasub_data <- reactive({
        metasub_data <- data_preprocess_f(train_f(),classTarget_in(),hierarchy_in(),remove_small())
        return(metasub_data)
    })
    
    # feature elimination  
    featureElim <- reactive({
        range_1 <- abundance_r()[1]
        range_2 <- abundance_r()[2]
        
        if(input$program == "bb" ){
          
          if(input$file_num_2 == "two_file" ){
            req(input$f_metadata_2)
            data_metadata <- read.csv(input$f_metadata_2$datapath, 
                                      header = TRUE)
            validate(
              need('csv' %in%  strsplit(input$f_metadata_2$datapath,".",fixed = T)[[1]], "The metadata file you upload should be in.csv format.")
            )
            range_1 <- dim(data_metadata)[2] + as.numeric(range_1) - 1
            range_2 <- dim(data_metadata)[2] + as.numeric(range_2) - 1
          }}
        if(input$program == "cc" ){
          
          if(input$file_num_3 == "two_file" ){
            req(input$f_metadata_3)
            data_metadata <- read.csv(input$f_metadata_3$datapath, 
                                      header = TRUE)
            validate(
              need('csv' %in%  strsplit(input$f_metadata_3$datapath,".",fixed = T)[[1]], "The metadata file you upload should be in.csv format.")
            )
            range_1 <- dim(data_metadata)[2] + as.numeric(range_1) - 1
            range_2 <- dim(data_metadata)[2] + as.numeric(range_2) - 1
          }}
        
        featureElim <- featureElimination_f(metasub_data(),classTarget_in(),range_1,range_2, subsets_in()) # 44:3712 for MetaSub
        return(featureElim)
    })
    
    # Original prediction result file  
    prediction_output_before <- reactive({
        if(input$program == "bb" ){
          test_file <- test_f()
          prediction_output_before <- New_model_prediction(train_f(),test_f(),classTarget_in(),featureElim()$optVariables,hierarchy_in())
        }
        if(input$program == "cc" ){
          prediction_output_before <- model_accuracy_f(metasub_data(),featureElim()$optVariables,classTarget_in(),hierarchy_in())
        }
        return(prediction_output_before)
    })
      
    # Prediction result file after move points to land or marine
    prediction_output <- reactive({
        prediction_output <- prediction_output_before()
        hierarchy_here <- hierarchy_in()
        if(input$program == "bb" ){
          if (input$pull_2 == "marine_2"){
            df <- prediction_output[[1]]
            df_new <- pull_marine(df,hierarchy_here)
            prediction_output[[1]] <- df_new
          }else if(input$pull_2 == "land_2"){
            df <- prediction_output[[1]]
            df_new <- pull_land(df,hierarchy_here)
            prediction_output[[1]] <- df_new
          }
        }
        if(input$program == "cc" ){
          if (input$pull_3 == "marine_3"){
            df <- prediction_output[[1]]
            df_new <- pull_marine(df,hierarchy_here)
            prediction_output[[1]] <- df_new
          }else if(input$pull_3 == "land_3"){
            df <- prediction_output[[1]]
            df_new <- pull_land(df,hierarchy_here)
            prediction_output[[1]] <- df_new
          }
        }
        return(prediction_output)
    })
      
    # Output the optimal features found in feature elimination step 
    v <- reactive({
        text <- featureElim()
        v <- v_f(text)
        return(v)
    })
      
    # Output the accuracy of each feature subset size  
    git_subset <- reactive({
        text <- featureElim()
        git_subset <- git_subset_f(text)
        return(git_subset)
    })
      
# Visualization of results and provide output files
    
    # Output of Module one
    observeEvent(input$acbutton_3,{
        output$predicted_map_3 <- renderPlot({
          x_ran <- input$xlim_3
          y_ran <- input$ylim_3
          df <- prediction_output()[[1]]
          hierarchy_in <- hierarchy_in()
          classTarget_in <- classTarget_in()
          
          plot_map(df,hierarchy_in,classTarget_in,x_ran,y_ran)
        })
        
        output$predicted_accuracy_3 <- renderPlot({
          classTarget_in <- classTarget_in()
          plot_accuracy(prediction_output()[[1]],classTarget_in)
        })
        
        output$downloadData_3 <- downloadHandler(
          filename = function(){
            return("Prediction_results.csv")
          },
          content = function(file){
            write.csv(prediction_output()[[1]],file)
          })
        
        output$download_optfeatures_3 <- downloadHandler(
          filename = function(){
            return("Optimal_features.csv")
          },
          content = function(file){
            write.csv(v(),file)
          })
        
        output$download_featuresub_3 <- downloadHandler(
          filename = function(){
            return("Features_subsets_accuracy.csv")
          },
          content = function(file){
            write.csv(git_subset(),file)
          })
      
        output$downloadMap_3 <- downloadHandler(
          filename = function(){
            return("Prediction_map.png")
          },
          content = function(file){

            MetasubDataPreds <- prediction_output()[[1]]
            hierarchy_in <- hierarchy_in()
            classTarget_in <- classTarget_in()
            x_ran <- input$xlim_3
            y_ran <- input$ylim_3
            
            png(file, width = 13,height = 8, units = 'in', res = 600)
            par(mai=c(2,1,0.5,0.5), mar=par()$mar+c(3,0,0,0))
            
            map <- rworldmap::getMap(resolution = "coarse")
            palette <-c( "darkorchid4","gold2","dodgerblue3","brown","orangered2","mediumspringgreen","deeppink2")
            
            plot(map, xlim = x_ran,ylim = y_ran,col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
            title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
            
            for (i in 1:(length(hierarchy_in)-2)){
              MetasubDataPreds[,hierarchy_in[i]] <- factor(MetasubDataPreds[,hierarchy_in[i]])
            }
            
            MetasubDataPreds$cityPred <- factor(MetasubDataPreds$cityPred)
            
            # get the coordination of continent
            continent_list <- c("east_asia","europe","middle_east","north_america","oceania","south_america", "sub_saharan_africa",'Australia','North_America','Europe','Africa','South_america','Asia')
            continent_lats <- c(55,69,8,40,-40,-10,-5,-40,40,69,-5,-10,55)
            continent_longs <- c(125,0,60,-130,140,-80,5,140,-130,0,5,-80,125)
            continent_position <- data.frame(cbind(continent_list,continent_lats,continent_longs))
            
            label_continent <- c()
            flag <- 0
            for ( i in 1:length(levels(MetasubDataPreds[,hierarchy_in[1]]))){
              this_continent <- levels(MetasubDataPreds[,hierarchy_in[1]])[i]
              label_continent <- c(label_continent,this_continent)
              find_lats <- MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,][,"latPred"]
              find_longs <- MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,][,"longPred"]
              
              #plot predicted co-ordinates
              points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.2,xlim = c(-165,168))
              
              #plot city prediction accuravy by continent as pies
              correctly_pred <-  mean(as.numeric(as.character(MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,"cityPred"])== 
                                                                as.character(MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,classTarget_in]))) 
              incorrectly_pred <- (1 - correctly_pred)
              
              if (this_continent %in% continent_position$continent_list){
                add.pie(z = c(correctly_pred, incorrectly_pred), x = as.numeric(continent_position[continent_position$continent_list==this_continent,3]), y = as.numeric(continent_position[continent_list==this_continent,2])
                        ,edges=200,
                        radius=10,
                        col=c(palette[i],"black") , labels = ""
                )
              }else{
                add.pie(z = c(correctly_pred, incorrectly_pred), x = as.numeric(-150+flag), y = 105
                        ,edges=200,
                        radius=10,
                        col=c(palette[i],"black") , labels = ""
                )
                flag <- flag + 25
              }
            }
            legend(xpd=T,'bottom',inset=c(0,-0.25), label_continent, pch = "+", col = palette[1:length(label_continent)], cex = 1,n=5)
            box( col = 'black')
            map.axes()
            dev.off()
          })
        
        output$downloadAccuracy_3 <- downloadHandler(
          filename = function(){
            return("Prediction_accuracy.png")
          },
          content = function(file){

            MetasubDataPreds <- prediction_output()[[1]]

            classTarget_in <- classTarget_in()
            
            MetasubDataPreds[,classTarget_in] <- factor(MetasubDataPreds[,classTarget_in])
            
            bar_df1 <- data.frame(row.names = c("Overall",levels(MetasubDataPreds[,classTarget_in])))
            
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 100)
              bar_df1[i+1,"0 - 100km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] < 100)
              bar_df1[ 1,"0 - 100km"] <- overall_prop
            }
            
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 100 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 500)
              bar_df1[i+1,"100 - 500km"] <- prop
              
              overall_prop <-mean(MetasubDataPreds[,"Distance_from_origin"] > 100 & MetasubDataPreds[,"Distance_from_origin"] < 500)
              bar_df1[ 1,"100 - 500km"] <- overall_prop
            }
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 500 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 1000)
              bar_df1[i+1,"500 - 1000km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 500 & MetasubDataPreds[,"Distance_from_origin"] < 1000)
              bar_df1[ 1,"500 - 1000km"] <- overall_prop
            }
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 1000 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 2000)
              bar_df1[i+1,"1000 - 2000km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 1000 & MetasubDataPreds[,"Distance_from_origin"] < 2000)
              bar_df1[ 1,"1000 - 2000km"] <- overall_prop
            }
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 2000 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 3000)
              bar_df1[i+1,"2000 - 3000km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 2000 & MetasubDataPreds[,"Distance_from_origin"] < 3000)
              bar_df1[1,"2000 - 3000km"] <- overall_prop
            }
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 3000 )
              bar_df1[i+1,"> 3000km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 3000)
              bar_df1[ 1,"> 3000km"] <- overall_prop
            }
            
            size <- c()
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              size[i] <- length(which(MetasubDataPreds[,classTarget_in] == this_city))
            }
     
            png(file, width = 13,height = 8, units = 'in', res = 600)
            par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)

            bp <- barplot(t(bar_df1*100), space = 0,col=c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), 
                          names.arg=c("Overall",paste0(levels(MetasubDataPreds[,classTarget_in]),"  (",size,")"), axes = FALSE) , 
                          las =2, cex.names=.6, ylab = "", axisnames = F, axes = F)
            axis(side =2, pos = 0)
            mtext(text = c("Overall",paste0(levels(MetasubDataPreds[,classTarget_in])," (",size,")")), side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
            title(ylab="Proportion of sample predictions %", mgp=c(2,1,0),cex.lab=1)
            legend("topright",inset = c(-0.12,0.4), rev(c(colnames(bar_df1))), 
                   fill = rev(c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) ,
                   bty = 1, cex = 0.8)
           
            par(mar=c(5, 4, 4, 2) + 0.1)
            dev.off()
          })
        
        output$downloadModel <- downloadHandler(
          filename = function(){
            return("Prediction_model.Rda")
          },
          content = function(file){
            output_model <- prediction_output()[[2]]
            save(output_model,file=file)
          })
    })
    
    # Output of Module two
    observeEvent(input$acbutton_2,{
      output$predicted_map_2 <- renderPlot({
        df <- prediction_output()[[1]]
        x_ran <- input$xlim_2
        y_ran <- input$ylim_2
        plot_prediction(df,x_ran,y_ran)
      })
      
      output$downloadMap_2 <- downloadHandler(
        
        filename = function(){
          return("Prediction_map.png")
        },
        content = function(file){
          
          MetasubDataPreds <- prediction_output()[[1]]
          x_ran <- input$xlim_2
          y_ran <- input$ylim_2
          
          png(file, width = 13,height = 8, units = 'in', res = 600)
          par(mai=c(1,1,0.5,0.5))
          
          map <- rworldmap::getMap(resolution = "coarse")
          plot(map, xlim = x_ran,ylim = y_ran, col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
          title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
          find_lats <- MetasubDataPreds$latPred
          find_longs <- MetasubDataPreds$longPred
          points(find_longs,find_lats, type = "p",col = "purple", pch = "+", cex = 1.3)
          map.axes()
          dev.off()
        })
      
      output$downloadData_2 <- downloadHandler(
        
        filename = function(){
          return("Prediction_results.csv")
        },
        content = function(file){
          write.csv(prediction_output()[[1]],file)
        })
      
      output$download_optfeatures_2 <- downloadHandler(
        
        filename = function(){
          return("Optimal_features.csv")
        },
        content = function(file){
          write.csv(v(),file)
        })
      
      output$download_featuresub_2 <- downloadHandler(
        
        filename = function(){
          return("Features_subsets_accuracy.csv")
        },
        content = function(file){
          write.csv(git_subset(),file)
        })
      
      output$downloadModel_2 <- downloadHandler(
        filename = function(){
          return("Prediction_model.Rda")
        },
        content = function(file){
          output_model <- prediction_output()[[2]]
          save(output_model,file=file)
        })
    })
    
    # output of Module three
    observeEvent(input$acbutton_1,{
      output$predicted_map_1 <- renderPlot({
        x_ran <- input$xlim_1
        y_ran <- input$ylim_1
        plot_prediction(pull_MetasubDataPreds(),x_ran,y_ran)
      })
      
      output$downloadMap_1 <- downloadHandler(
        
        filename = function(){
          return("Prediction_map.png")
        },
        content = function(file){
          
          MetasubDataPreds <- pull_MetasubDataPreds()
          
          x_ran <- input$xlim_1
          y_ran <- input$ylim_1
          
          png(file, width = 13,height = 8, units = 'in', res = 600)
          
          par(mai=c(1,1,0.5,0.5))
          
          map <- rworldmap::getMap(resolution = "coarse")
          
          plot(map, xlim = x_ran,ylim = y_ran, col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
          title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
          
          find_lats <- MetasubDataPreds$latPred
          find_longs <- MetasubDataPreds$longPred
          
          points(find_longs,find_lats, type = "p",col = "purple", pch = "+", cex = 1.3)
          
          map.axes()
          dev.off()
        })
      
      output$downloadData_1 <- downloadHandler(
        filename = function(){
          return("Prediction_results.csv")
        },
        content = function(file){
          write.csv(pull_MetasubDataPreds(),file)
        }
      )
    })
}

shinyApp(ui=ui,server=server)
```

## 2 Microorganisms analysis in the urban environments  
### 2.1 Taxonomical profile annotation of taxa in MetaSub
Downloaded the raw csv file of microbial taxonomy profile from Microbe Directory v2.0 (https://github.com/dcdanko/MD2) and GTDB (https://data.gtdb.ecogenomic.org/releases/release202/202.0/), used python to extract the annotation information of taxa to the file "MetaSub_taxonomic_level.txt". Then manually fixed the heterotypic synonym of taxa in the taxonomical profile.  
Resaved the "MetaSub_taxonomic_level.txt" into CSV format "MetaSub_taxonomic_level.csv"

```python

from collections import defaultdict

taxa = open("metasub_taxa_abundance.csv",'r',encoding='UTF-8')
directory = open("microbe-directory.csv",'r',encoding='UTF-8')
fileout = open("MetaSub_taxonomic_level.txt",'w',encoding='UTF-8')

bac120 = open("bac120_taxonomy.tsv",'r',encoding='UTF-8')
bac120_dic =  defaultdict(list)
for line in bac120:
    level = line.strip().split("\t")[1]
    taxa_l = level.split(";")
    species = taxa_l[6].split("__")[1] # species level
    if species in bac120_dic.keys():
        continue
    for i in taxa_l:
        t = i.split("__")[1]
        bac120_dic[species].append(t) # d p c o f g s

bac120_dic2 = {}
for line in bac120:
    level = line.strip().split("\t")[1]
    species = level.split(";")[6] # species level
    species = species.split("__")[1]
    if species in bac120_dic2:
        continue
    bac120_dic2[species] = level

ar122 = open("ar122_taxonomy.tsv",'r',encoding='UTF-8')
ar122_dic = defaultdict(list)
for line in ar122:
    level = line.strip().split("\t")[1]
    taxa_l = level.split(";")
    species = taxa_l[6].split("__")[1] # species level
    if species in ar122_dic:
        continue
    for i in taxa_l:
        t = i.split("__")[1]
        ar122_dic[species].append(t) # d p c o f g s

ssu = open("ssu_silva_taxonomy",'r',encoding='UTF-8')
ssu_dic = defaultdict(list)
for line in ssu:
    taxa_l = line.strip().split(";")
    if taxa_l[-1] in ssu_dic:
        continue
    for i in taxa_l:
        ssu_dic[taxa_l[-1]].append(i)

ncbi = open("ncbi_taxonomy",'r',encoding='UTF-8')
ncbi_dic = defaultdict(list)
for line in ncbi:
    taxa_l = line.strip().split(";")
    s_ncbi = taxa_l[-1]
    if "__" in s_ncbi:
        s_ncbi = s_ncbi.split("__")[1]
        if s_ncbi in ncbi_dic:
            continue
        for i in taxa_l:
            t = i.split("__")[1]
            ncbi_dic[s_ncbi].append(t)

ar_ssu = open("ar122_ssu_silva_taxonomy",'r',encoding='UTF-8')
ar_ssu_dic = defaultdict(list)
for line in ar_ssu:
    taxa_l = line.strip().split(";")
    species = taxa_l[-1]
    if species in ar_ssu_dic:
        continue
    for i in taxa_l:
        ar_ssu_dic[species].append(i)

ar_ncbi = open("ar122_ncbi_taxonomy",'r',encoding='UTF-8')
ar_ncbi_dic = defaultdict(list)
for line in ar_ncbi:
    taxa_l = line.strip().split(";")
    s_ncbi = taxa_l[-1]
    if "__" in s_ncbi:
        s_ncbi = s_ncbi.split("__")[1]
        if s_ncbi in ar_ncbi_dic:
            continue
        for i in taxa_l:
            t = i.split("__")[1]
            ar_ncbi_dic[s_ncbi].append(t)

directory_dic =  defaultdict(list)
for line in directory:
    line = line.strip().split(",")
    if line[7] in directory_dic:
            continue
    for n in range(1,8):
        directory_dic[line[7]].append(line[n]) # kingdom: directory_dic[taxa][0]
    directory_dic[line[7]].append(line[22]) # animal_pathogen
    directory_dic[line[7]].append(line[28]) # plant_pathogen

taxa_list = taxa.readline().strip().split(',')[1:]

anno_file = [directory_dic,bac120_dic,ssu_dic,ncbi_dic,ar122_dic,ar_ssu_dic,ar_ncbi_dic]
flag = {1:'directory_dic',2:'bac120_dic',3:'ssu_dic',4:'ncbi_dic',5:'ar122_dic',6:'ar_ssu_dic',7:'ar_ncbi_dic'}

for t in taxa_list:
    flag_count = 0
    for i in anno_file:
        flag_count += 1
        if t in i:
            anno = ''
            for feature in i[t]:
                if anno:
                    anno = "{}\t{}".format(anno,feature)
                else:
                    anno = feature
            print('{}\t{}\t{}'.format(flag[flag_count],t,anno),file=fileout)
            break

taxa.close()
directory.close()
fileout.close()
bac120.close()
ar122.close()
ssu.close()
ncbi.close()
ar_ssu.close()
ar_ncbi.close()
```

### 2.2 Pathogenicity annotation of taxa
Used R to extract pathogen in MetaSUB and annotated the host of them (Flag in pathogen_type: Animal pathogen-1, plant pathogen-2, dual-kingdom-3 )  
```r
annotated_taxa <- read.csv('MetaSub_taxonomic_level.csv')

for (i in 1:length(annotated_taxa$taxa)){
  animal <- annotated_taxa$Animal_pathoge[i]
  plant <- annotated_taxa$Plant_pathoge[i]
  
  if(is.na(animal)){
    animal <- 999
  }
  if(is.na(plant)){
    plant <- 999
  }
  
  if(animal == 1 && plant == 1){
    annotated_taxa[i,'pathogen_type'] <- 3
  }else if(animal ==1 && plant == 0){
    annotated_taxa[i,'pathogen_type'] <- 1
  }else if(animal ==1 && plant == 999){
    annotated_taxa[i,'pathogen_type'] <- 1
  }else if(animal ==0 && plant == 1){
    annotated_taxa[i,'pathogen_type'] <- 2
  }else if(animal ==999 && plant == 1){
    annotated_taxa[i,'pathogen_type'] <- 2
  }
}

pathogen_taxa <- annotated_taxa[which(annotated_taxa$pathogen_type != 'NA'),]
write.csv(pathogen_taxa,'pathogen_MetaSub_in_Microbe_Directory_new.csv',row.names = F)
```

### 2.3 Microorganisms distribution in MetaSUB
Used R to detect the microorganisms (taxa and pathogen) distribution in MetaSUB.  
```r
# Loading the required packages
library(ggplot2)
library(forcats)
library(ggpmisc)
library(Rcpp)
library(ggrepel)
library(ggforce)
library(scales)
library(vegan)
library(gtable)
library(grid)
library(RColorBrewer)

#### 1 Pre-process on dataset 

# These input files are provided in Github
complete_meta <- read.csv(file = "complete_metadata.csv", header = TRUE)
taxa_abund <-read.csv("metasub_taxa_abundance.csv", header = T)
taxa_abund <- unique(taxa_abund)

# Merge bacterial and meta data
metasub_data <- merge(complete_meta,taxa_abund,by.x="uuid",by.y="uuid")

# Remove control samples
control_samples <- c( which(metasub_data$city %in% c("control", "other_control","neg_control","other","pos_control")), which(metasub_data$control_type %in% c("ctrl cities","negative_control","positive_control"))) 
metasub_data <- droplevels(metasub_data[-c(control_samples), ])

# Re-label london boroughs 
metasub_data$city[metasub_data$city %in% c("kensington","islington")] <- "london" 
metasub_data <- droplevels(metasub_data)

# Remove dubiously labelled samples. 
remove_samples <- which(metasub_data$city %in%  "antarctica")
metasub_data <- droplevels(metasub_data[-c(remove_samples), ])

# Correction of identified misslabelling of data 
metasub_data$latitude[metasub_data$city == "kyiv"] <- metasub_data$city_latitude[metasub_data$city == "kyiv"]
metasub_data$longitude[metasub_data$city == "kyiv"] <- metasub_data$city_longitude[metasub_data$city == "kyiv"]
metasub_data$continent[metasub_data$city == "porto"] <- "europe"

metasub_data[is.na(metasub_data$latitude),]$latitude <- metasub_data[is.na(metasub_data$latitude),]$city_latitude
metasub_data[is.na(metasub_data$longitude),]$longitude <- metasub_data[is.na(metasub_data$longitude),]$city_longitude

# Correction to some incorrect city co-ords for a few london samples
metasub_data[metasub_data$city == "london",]$city_latitude <- 51.50853
metasub_data[metasub_data$city == "london",]$city_longitude <- -0.12574

empty <- which(metasub_data$continent == "")
metasub_data1 <- metasub_data 
metasub_data <- metasub_data[-c(empty),]

metasub_data$city <- factor(metasub_data$city)

levels(metasub_data$city)


#### 2 Calculate the occupancy of each taxa

# Create a data form to store the number of times that taxa appears in each city
bar_df <- data.frame(row.names = c(levels(metasub_data$city)))

taxa_list <- colnames(metasub_data)[43:3711]

for (i in 1:length(levels(metasub_data$city))){
  for (j in 1: length(taxa_list)){
    city <- levels(metasub_data$city)[i]
    taxa <- taxa_list[j]
    count <- sum(metasub_data[metasub_data$city == city,][,taxa] != 0)
    bar_df[city,taxa] <- count
  }
}

# Calculate the global total occupancy times
Overall <- c()
for (i in taxa_list){
  count_sum <- sum(bar_df[,i])
  Overall <- c(Overall,count_sum)
}

taxa_count <- rbind(bar_df,Overall)
row.names(taxa_count) <- c(levels((metasub_data$city)),"Overall")

metasub_occupancy_file <- taxa_count

# Calculate percentage of occupancy for each taxa
bar_df_t <- data.frame(t(taxa_count),stringsAsFactors = F)
bar_df_t <- as.data.frame(lapply(bar_df_t,as.numeric))
bar_df_t <- data.frame(taxa_list,bar_df_t)


for (i in 1:dim(bar_df_t)[1]){
  # the percentage of sample that taxa appeared
  bar_df_t[i,'Percentage_overall'] <- round(bar_df_t[i,'Overall']/4135*100,2)
  # the number of city that taxa appeared
  bar_df_t[i,'Overall_city'] <- sum(bar_df_t[i,][,2:54] != 0)
  # the percentage of city that taxa appeared
  bar_df_t[i,'Percentage_city'] <- round(bar_df_t[i,'Overall_city']/53*100,2)
}

# Remove taxa with 0 global relative abundance
empty <- which(bar_df_t$Overall == 0)
order_taxa <- bar_df_t[-c(empty),] # 3660 taxa

t_metasub_occupancy_file <- order_taxa

#### 3 Calculate and plot the occupancy distribution of taxa

# Occupancy distribution in samples
count_sample2 <- order_taxa[order(order_taxa$Percentage_overall,decreasing=T),]

percentage2 <- matrix(ncol=2,nrow=11)

percentage2[1,1] <- '≥95%'
percentage2[2,1] <- '90-95%'
percentage2[3,1] <- '80-90%'
percentage2[4,1] <- '60-80%'
percentage2[5,1] <- '40-60%'
percentage2[6,1] <- '20-40%'
percentage2[7,1] <- '10-20%'
percentage2[8,1] <- '5-10%'
percentage2[9,1] <- '1-5%'
percentage2[10,1] <- '0.02-1%'
percentage2[11,1] <- '=0.02%\n(1 sample)'

percentage2[1,2] <- sum(count_sample2[,'Percentage_overall'] >= 95)
percentage2[2,2] <- sum(count_sample2[,'Percentage_overall'] >= 90 & count_sample2[,"Percentage_overall"] < 95)
percentage2[3,2] <- sum(count_sample2[,"Percentage_overall"] >= 80 & count_sample2[,"Percentage_overall"] < 90)
percentage2[4,2] <- sum(count_sample2[,"Percentage_overall"] >= 60 & count_sample2[,"Percentage_overall"] < 80)
percentage2[5,2] <- sum(count_sample2[,"Percentage_overall"] >= 40 & count_sample2[,"Percentage_overall"] < 60)
percentage2[6,2] <- sum(count_sample2[,"Percentage_overall"] >= 20 & count_sample2[,"Percentage_overall"] < 40)
percentage2[7,2] <- sum(count_sample2[,"Percentage_overall"] >= 10 & count_sample2[,"Percentage_overall"] < 20)
percentage2[8,2] <- sum(count_sample2[,"Percentage_overall"] >= 5 & count_sample2[,"Percentage_overall"] < 10)
percentage2[9,2] <- sum(count_sample2[,"Percentage_overall"] >= 1 & count_sample2[,"Percentage_overall"] < 5)
percentage2[10,2] <- sum(count_sample2[,"Percentage_overall"] > 0.02 & count_sample2[,"Percentage_overall"] < 1)
percentage2[11,2] <- sum(count_sample2[,"Percentage_overall"] == 0.02 )

percentage2 <- data.frame(percentage2)
percentage2$X2 <- as.numeric(percentage2$X2)
percentage2$X1 <- fct_inorder(percentage2$X1)

ggplot(percentage2,aes(x=X1,y=X2))+
  geom_bar(stat='identity',fill = "lightblue")+
  geom_text(aes(label= paste0(X2,' (',round(X2/3660*100),'%)'),vjust = -0.5),
            position=position_dodge(width=1))+
  ggtitle('Taxa distribution over samples (3660 taxa, 4135 samples)')+
  ylim(0,900)+
  labs(x = 'Percentage of samples (%)', y = 'Number of taxa')+
  theme_bw()+
  theme(axis.title =  element_text(size=12,face = "bold"),
        axis.text =   element_text(size=12))


# Occupancy distribution in city  
percentage2 <- matrix(ncol=2,nrow=9)

percentage2[1,1] <- '100%'
percentage2[2,1] <- '90-100%'
percentage2[3,1] <- '80-90%'
percentage2[4,1] <- '60-80%'
percentage2[5,1] <- '40-60%'
percentage2[6,1] <- '20-40%'
percentage2[7,1] <- '10-20%'
percentage2[8,1] <- '2-10%\n(2-6 cities)'
percentage2[9,1] <- '<2%\n(1 city)'

percentage2[1,2] <- sum(count_sample2[,'Percentage_city'] == 100)
percentage2[2,2] <- sum(count_sample2[,'Percentage_city'] >= 90 & count_sample2[,"Percentage_city"] < 100)
percentage2[3,2] <- sum(count_sample2[,"Percentage_city"] >= 80 & count_sample2[,"Percentage_city"] < 90)
percentage2[4,2] <- sum(count_sample2[,"Percentage_city"] >= 60 & count_sample2[,"Percentage_city"] < 80)
percentage2[5,2] <- sum(count_sample2[,"Percentage_city"] >= 40 & count_sample2[,"Percentage_city"] < 60)
percentage2[6,2] <- sum(count_sample2[,"Percentage_city"] >= 20 & count_sample2[,"Percentage_city"] < 40)
percentage2[7,2] <- sum(count_sample2[,"Percentage_city"] >= 10 & count_sample2[,"Percentage_city"] < 20)
percentage2[8,2] <- sum(count_sample2[,"Percentage_city"] >= 2 & count_sample2[,"Percentage_city"] < 10)
percentage2[9,2] <- sum(count_sample2[,"Percentage_city"] < 2 )

percentage2 <- data.frame(percentage2)
percentage2$X2 <- as.numeric(percentage2$X2)
percentage2$X1 <- fct_inorder(percentage2$X1)

ggplot(percentage2,aes(x=X1,y=X2))+
  geom_bar(stat='identity',fill = "lightblue")+
  geom_text(aes(label=paste0(X2,' (',round(X2/3660*100),'%)'),vjust = -0.5),
            position=position_dodge(width=1))+
  ylim(0,550)+
  ggtitle('Taxa distribution over cities (3660 taxa, 53 cities)')+
  labs(x = 'Percentage of cities (%)', y = 'Number of taxa')+
  theme_bw()+
  theme(axis.title = element_text(size=12,face = "bold"),
        axis.text = element_text(size=12)) 


#### 4 Find the most common taxa and phylum

count_sample <- bar_df_t[order(bar_df_t$Percentage_overall,decreasing=T),]

# Most common top 50 microbiome in MetaSub dataset (sample and city)
count_sample_10 <- count_sample[1:50,]

# Original name of each taxa
taxa_name <- read.csv("metasub_taxa_abundance.csv",check.names = F) 
taxa_name <- as.data.frame(colnames(taxa_name)[2:3670])
colnames(taxa_name) <- 'uuid'
taxa_name$make_id <- make.names(taxa_name$uuid)

for (i in 1:length(count_sample_10$taxa_list)){
  taxa_new <- count_sample_10$taxa_list[i]
  count_sample_10[i,'original_id'] <- taxa_name[taxa_name$make_id == taxa_new,'uuid']
}

count_sample_10$original_id <- factor(count_sample_10$original_id)
count_sample_10$original_id <- fct_inorder(count_sample_10$original_id)
count_sample_10$original_id <- fct_relevel(count_sample_10$original_id,
                                           rev(levels(count_sample_10$original_id)))     
count_sample_10$Percentage_city <- as.numeric(count_sample_10$Percentage_city)

ggplot(count_sample_10,aes(x=original_id,y=Percentage_overall))+
  geom_bar(stat='identity',aes(fill=factor(Percentage_city)),position = "dodge",width=0.5)+
  xlab('Taxa')+
  ylab('Percentage of samples (%)')+
  labs(fill="Percentage of cities (%)")+
  ggtitle('Most common taxa in MetaSub (top 50)')+
  coord_flip()+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_fill_manual(values = c("darkblue","lightblue"))+
  theme_bw()


# Phylum occupancy distribution

annotated_taxa <- read.csv('MetaSub_taxonomic_level.csv')

metasub_data$city <- factor(metasub_data$city)
annotated_taxa$Phylum <- factor(annotated_taxa$Phylum)

annotated_taxa$Order <- factor(annotated_taxa$Order)
annotated_taxa$Family <- factor(annotated_taxa$Family)

levels(annotated_taxa$Order)
levels(annotated_taxa$Family)

phylum_list <- levels(annotated_taxa$Phylum)

# Create a data form to store the number of times that phylum appears in each city
bar_df <- data.frame(row.names = c(levels(metasub_data$city)))

for (j in 1:length(levels(annotated_taxa$Phylum))){
  phylum <- levels(annotated_taxa$Phylum)[j]
  phylum_taxa <- make.names(annotated_taxa[annotated_taxa$Phylum==phylum,'taxa'])
  position <- which(colnames(metasub_data) %in% phylum_taxa)
  
  for (i in 1:dim(bar_df)[1]){
    city <- row.names(bar_df)[i]
    count <- 0
    list0 <- metasub_data[metasub_data$city==city,position]
    
    if (length(position) < 2){
      bar_df[city,phylum] <- sum(list0 > 0)
      
    }else{
      length0 <- dim(list0)[1]
      
      for (k in 1:length0){
        for (m in list0[k,]){
          if (m > 0){
            count <- count + 1
            break
          }
        }
      }
      bar_df[city,phylum] <- count
    }
  }
}

Overall <- c()
for (i in levels(annotated_taxa$Phylum)){
  count_sum <- sum(bar_df[,i])
  Overall <- c(Overall,count_sum)
}

bar_df2 <- rbind(bar_df,Overall)
row.names(bar_df2) <- c(row.names(bar_df),"Overall")

bar_df_t <- data.frame(t(bar_df2))
bar_df_t <- data.frame(taxa_list = row.names(bar_df_t),bar_df_t)

for (i in 1:dim(bar_df_t)[1]){
  bar_df_t[i,'Percentage_overall'] <- round(bar_df_t[i,'Overall']/4135*100,2)
  bar_df_t[i,'Overall_city'] <- sum(bar_df_t[i,][,2:54] != 0)
  bar_df_t[i,'Percentage_city'] <- round(bar_df_t[i,'Overall_city']/53*100,2)
}

count_sample <- bar_df_t[order(bar_df_t$Percentage_overall,bar_df_t$Percentage_city,decreasing=T),]

count_sample_10 <- count_sample
count_sample_10$taxa_list <- factor(count_sample_10$taxa_list)
count_sample_10$taxa_list <- fct_inorder(count_sample_10$taxa_list)
count_sample_10$taxa_list <- fct_relevel(count_sample_10$taxa_list,
                                         rev(levels(count_sample_10$taxa_list)))     
count_sample_10$Percentage_city <- as.numeric(count_sample_10$Percentage_city)

ggplot(count_sample_10,aes(x=taxa_list,y=Percentage_overall))+
  geom_bar(stat='identity',aes(fill=Percentage_city),position = "dodge",width=0.5)+
  xlab('Taxa')+
  ylab('Percentage of samples (%)')+
  labs(fill="Percentage of cities (%)")+
  ggtitle('Phylum occupancy in MetaSub (43 phylum)')+
  coord_flip()+
  theme_bw()


#### 5 Calculate city average relative abundance and global average relative abundance of each taxa

# Calculate city average relative abundance 
bar_df_abundance <- data.frame(row.names = c(levels(metasub_data$city)))

for (i in 1:length(levels(metasub_data$city))){
  for (j in 1: length(taxa_list)){
    city <- levels(metasub_data$city)[i]
    taxa <- taxa_list[j]  
    count <- sum(metasub_data[metasub_data$city == city,][,taxa])/length(metasub_data[metasub_data$city == city,][,taxa])
    bar_df_abundance[city,taxa] <- count
  }
}

bar_df_abundance <- bar_df_abundance[1:53,]
abundance_bar <- bar_df_abundance

# Calculate global average relative abundance
global_proportion <- c()

taxa_list <- colnames(abundance_bar)

for (i in 1:(length(taxa_list))){
  count <- sum(abundance_bar[,i])/length(abundance_bar[,i])
  global_proportion <- c(global_proportion,as.numeric(count))
}

row_name <- rownames(abundance_bar)
abundance_bar <- rbind.data.frame(global_proportion,abundance_bar)
rownames(abundance_bar) <- c('Overall',row_name)

metasub_abundance_file <- abundance_bar

#### 6 Relative abundance of common and rare taxa groups for each city

common_taxa_1 <- order_taxa$taxa_list[which(order_taxa[,'Percentage_overall'] >= 95)] # most common taxa
rare_taxa_1 <- order_taxa$taxa_list[which(order_taxa[,'Percentage_overall'] <= 5)] # rarest taxa

for (i in 1:dim(bar_df_abundance)[1]){
  bar_df_abundance[i,'common_95_p'] <- sum(as.numeric(bar_df_abundance[i,common_taxa_1]))
  bar_df_abundance[i,'rare_5_p'] <- sum(as.numeric(bar_df_abundance[i,rare_taxa_1]))
}

# calculate taxa number in each city
city_taxa_n <- c()
for (i in 1:length(abundance_bar$X)){
  city_taxa_n <- c(city_taxa_n, sum(abundance_bar[i,][,taxa_list] != 0))
}

# Change city name
country_name <- c()
for (country in bar_df_abundance$X){
  if ('_' %in% strsplit(country,'')[[1]]){
    country <- strsplit(country,'_')[[1]]

    for (a in 1:length(country)){
      country_a <- strsplit(country[a],'')[[1]]
      country_a[1] <- toupper(country_a[1])
      country_a <- paste(country_a, collapse='')
      country[a] <- country_a}
    country <- paste(country, collapse=' ')
  }else{
    country <- strsplit(country,'')[[1]]
    country[1] <- toupper(country[1])
    country <- paste(country, collapse='')
  }
  country_name <- c(country_name,country)
}
country_name[1] <- 'Global'

bar_df2 <- bar_df_abundance[,c(3671,3672)]

# png("common_rare_abundance_MetaSub_0308.png", width = 13,height = 8, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(3,0,4,0))
bp <- barplot(t(bar_df2*100), col=c("skyblue","red"),
              names.arg= bar_df_abundance$X,
              args.legend = list(x = "bottom", inset=c(-0.5,0)), las =2,
              cex.names=.7,ylab = "", axisnames = F,axes = F, space =0)
axis(side =2, pos = -1)
mtext(text = paste0(country_name,' (',city_taxa_n,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 0.8)
title(ylab="Taxa relative abundance %", mgp=c(1,0,0),cex.lab=1.2)
legend("top",inset = c(0,-0.2), c('Rare taxa (bottom 45.3%) - found in <=5% samples','Common taxa (top 0.9%) - found in >=95% samples'), fill = rev(c("skyblue","red")) , bty = 1, cex = 1,ncol=2)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()


## correlation: number of taxa and abundance (most common taxa)

n_a_cor <- cbind.data.frame(city_taxa_n,bar_df2$common_95_p)
n_a_cor <- n_a_cor[-1,]

ggplot(data = n_a_cor,aes(x=city_taxa_n,y=(bar_df2$common_95_p[-1])*100))+
  geom_point(size=2)+
  ylab('Relative abundance of most common taxa %')+
  xlab('Number of taxa in city')+
  geom_smooth (method = lm,linetype=1,se=FALSE,span=1)+
  theme_bw()+
  ggpubr::stat_cor(aes(color = Pathogen), label.x = 2050,label.y = 50,color = 'black')


#### 7 Pathogenic microorganisms distribution

# The pathogenicity of taxa are annotated by Microbiome Directory using python script

pathogen <- read.csv('pathogen_MetaSub_in_Microbe_Directory_new.csv',header = TRUE)

# extract pathogens
pathogen$Species <- make.names(pathogen$Species)
pathogen_list <- pathogen$Species

taxa_name_all <- colnames(taxa_count)
pathogen_position <- c(which(taxa_name_all %in% pathogen_list))
taxa_pathogen <- taxa_count[,c(1,pathogen_position)] # datasets just contains pathogens

bar_df_t <- data.frame(t(taxa_pathogen),stringsAsFactors = F)
taxa_list <- row.names(bar_df_t)[2:311] # pathogen name
colnames(bar_df_t) <- bar_df_t[1,]
bar_df_t <- bar_df_t[-1,]

bar_df_t <- as.data.frame(lapply(bar_df_t,as.numeric))
bar_df_t <- data.frame(taxa_list,bar_df_t)

for ( i in 1:length(bar_df_t$taxa_list)){
  taxa <- bar_df_t$taxa_list[i]
  p_type <- pathogen[make.names(pathogen$Species) == taxa, 'pathogen']
  bar_df_t[i,'pathogen_type'] <- p_type
  bar_df_t[i,'Percentage_overall'] <- round(bar_df_t[i,'Overall']/4135*100,2)
  bar_df_t[i,'Overall_city'] <- sum(bar_df_t[i,][,2:54] != 0)
  bar_df_t[i,'Percentage_city'] <- round(bar_df_t[i,'Overall_city']/53*100,2)
}

count_sample2 <- bar_df_t[order(bar_df_t$Percentage_overall,decreasing=T),]

percentage2 <- matrix(ncol=5,nrow=11)

percentage2[1,1] <- '≥95%'
percentage2[2,1] <- '90-95%'
percentage2[3,1] <- '80-90%'
percentage2[4,1] <- '60-80%'
percentage2[5,1] <- '40-60%'
percentage2[6,1] <- '20-40%'
percentage2[7,1] <- '10-20%'
percentage2[8,1] <- '5-10%'
percentage2[9,1] <- '1-5%'
percentage2[10,1] <- '0.02-1%'
percentage2[11,1] <- '=0.02%\n(1 sample)'

# calculate occupancy of pathogens with different hosts
for (i in 1:3){
  percentage2[1,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,'Percentage_overall'] >= 95)
  percentage2[2,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,'Percentage_overall'] >= 90 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 95)
  percentage2[3,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] >= 80 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 90)
  percentage2[4,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] >= 60 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 80)
  percentage2[5,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] >= 40 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 60)
  percentage2[6,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] >= 20 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 40)
  percentage2[7,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] >= 10 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 20)
  percentage2[8,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] >= 5 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 10)
  percentage2[9,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] >= 1 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 5)
  percentage2[10,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] > 0.02 & count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] < 1)
  percentage2[11,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_overall"] == 0.02 )
}

for (i in 1:11){
  percentage2[i,5] <- sum(as.numeric(percentage2[i,2:4]))
}

percentage2 <- data.frame(percentage2)
percentage2$X2 <- as.numeric(percentage2$X2)
percentage2$X3 <- as.numeric(percentage2$X3)
percentage2$X4 <- as.numeric(percentage2$X4)
percentage2$X5 <- as.numeric(percentage2$X5)

bar_df <- rbind.data.frame(cbind(name = percentage2$X1,type = rep('animal',11),c = percentage2$X2,a=percentage2$X5),
                           cbind(name = percentage2$X1,type = rep('both',11),c = percentage2$X4,a=percentage2$X5),
                           cbind(name = percentage2$X1,type = rep('plant',11),c = percentage2$X3,a=percentage2$X5))

bar_df$name <- fct_inorder(bar_df$name)
bar_df$type <- factor(bar_df$type)
bar_df$c <- as.numeric(bar_df$c)
bar_df$a <- as.numeric(bar_df$a)

ggplot(bar_df,aes(x=name,y=c,fill=type))+
  geom_bar(stat = 'identity')+
  ylim(0,60)+
  geom_text(aes(label= c(paste0(percentage2$X5,' (',round(percentage2$X5/310*100),'%)'),rep('',22)),vjust = -0.5,y = c(percentage2$X5,rep(0,22)) + 0.5),
            position=position_dodge(0))+
  ggtitle('Distribution of Metasub pathogen in samples (310 pathogens, 4135 samples)')+
  labs(x = 'Percentage of samples (%)', y = 'Number of pathogen')+
  theme_bw()+
  theme(axis.title =  element_text(size=12,face = "bold"),
        axis.text =   element_text(size=12))+ # adjust size of label
  scale_fill_discrete(name="Host",
                      breaks=c("animal", "both", "plant"),
                      labels=c("Animal", 'Both', "Plant"))
# ggsave('pathogen distribution in MetaSub.png',p)


## taxa distribution in city proportion 

percentage2 <- matrix(ncol=5,nrow=9)

percentage2[1,1] <- '100%'
percentage2[2,1] <- '90-100%'
percentage2[3,1] <- '80-90%'
percentage2[4,1] <- '60-80%'
percentage2[5,1] <- '40-60%'
percentage2[6,1] <- '20-40%'
percentage2[7,1] <- '10-20%'
percentage2[8,1] <- '2-10%\n(2-6 cities)'
percentage2[9,1] <- '<2%\n(1 city)'

for (i in 1:3){
  percentage2[1,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,'Percentage_city'] == 100)
  percentage2[2,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,'Percentage_city'] >= 90 & count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] < 100)
  percentage2[3,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] >= 80 & count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] < 90)
  percentage2[4,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] >= 60 & count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] < 80)
  percentage2[5,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] >= 40 & count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] < 60)
  percentage2[6,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] >= 20 & count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] < 40)
  percentage2[7,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] >= 10 & count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] < 20)
  percentage2[8,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] >= 2 & count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] < 10)
  percentage2[9,i+1] <- sum(count_sample2[count_sample2$pathogen_type == i,"Percentage_city"] < 2 )
}

for (i in 1:9){
  percentage2[i,5] <- sum(as.numeric(percentage2[i,2:4]))
}

percentage2 <- data.frame(percentage2)
percentage2$X2 <- as.numeric(percentage2$X2)
percentage2$X3 <- as.numeric(percentage2$X3)
percentage2$X4 <- as.numeric(percentage2$X4)
percentage2$X5 <- as.numeric(percentage2$X5)

bar_df <- rbind.data.frame(cbind(name = percentage2$X1,type = rep('animal',9),c = percentage2$X2,a=percentage2$X5),
                           cbind(name = percentage2$X1,type = rep('both',9),c = percentage2$X4,a=percentage2$X5),
                           cbind(name = percentage2$X1,type = rep('plant',9),c = percentage2$X3,a=percentage2$X5))

bar_df$name <- fct_inorder(bar_df$name)
bar_df$type <- factor(bar_df$type)
bar_df$c <- as.numeric(bar_df$c)
bar_df$a <- as.numeric(bar_df$a)

ggplot(bar_df,aes(x=name,y=c,fill=type))+
  geom_bar(stat = 'identity')+
  ylim(0,52)+
  geom_text(aes(label= c(paste0(percentage2$X5,' (',round(percentage2$X5/310*100),'%)'),rep('',18)),vjust = -0.4,y = c(percentage2$X5,rep(0,18)) + 0.5),
            position=position_dodge(0))+
  ggtitle('Distribution of Metasub pathogen in cities (310 pathogens, 53 cities)')+
  labs(x = 'Percentage of samples (%)', y = 'Number of pathogen')+
  theme_bw()+
  theme(axis.title =  element_text(size=12,face = "bold"),
        axis.text =   element_text(size=12))+ # adjust size of label
  scale_fill_discrete(name="Host",
                      breaks=c("animal", "both", "plant"),
                      labels=c("Animal", 'Both', "Plant"))
# ggsave('pathogen distribution in MetaSub city.png',p)


#### 8 Most common and rarest Pathogens

## most common pathogen

count_sample <- bar_df_t[order(bar_df_t$Percentage_overall,decreasing=T),]
count_sample_10 <- count_sample[1:20,] # top 20

count_sample_10$taxa_list <- factor(count_sample_10$taxa_list)
count_sample_10$taxa_list <- fct_inorder(count_sample_10$taxa_list)
count_sample_10$taxa_list <- fct_relevel(count_sample_10$taxa_list,
                                           rev(levels(count_sample_10$taxa_list)))     
count_sample_10$Percentage_city <- as.numeric(count_sample_10$Percentage_city)

ggplot(count_sample_10,aes(x=taxa_list,y=Percentage_overall))+
  geom_bar(stat='identity',aes(fill=factor(Percentage_city)),position = "dodge",width=0.5)+
  geom_text(stat='identity', mapping = aes(y = Percentage_overall+15, label = Percentage_overall),size = 3.5) +
  expand_limits(y=c(0,120))+
  xlab('Taxa')+
  ylab('Percentage of samples (%)')+
  labs(fill="Percentage of cities (%)")+
  ggtitle('Most common MetaSub pathogens (top 20)')+
  coord_flip()+
  theme_bw()+
  scale_fill_brewer(palette = "Pastel1")
# ggsave('metasub most common pathogen 20.png',p)

## rarest pathogen

count_sample <- bar_df_t[order(bar_df_t$Percentage_overall,bar_df_t$Percentage_city,decreasing=F),]
count_sample_10 <- count_sample[1:20,]

count_sample_10$taxa_list <- factor(count_sample_10$taxa_list)
count_sample_10$taxa_list <- fct_inorder(count_sample_10$taxa_list)
count_sample_10$taxa_list <- fct_relevel(count_sample_10$taxa_list,
                                           rev(levels(count_sample_10$taxa_list)))     
count_sample_10$Percentage_city <- as.numeric(count_sample_10$Percentage_city)

# find city with unique-taxon
for (i in 1:length(count_sample_10$taxa_list)){
  exsit_city <- c()
  for (j in city_name){
    if (count_sample_10[i,j] != 0){
      j2 <- country_name[which(city_name == j)]
      exsit_city <- c(exsit_city,j2)
    }
  }
  if (length(exsit_city) != 0){
    fill <- exsit_city[1]
    if (length(exsit_city) != 1){
      for (k in 2:length(exsit_city)){
        fill <- paste(fill,exsit_city[k],sep=',')
      }
    }
    count_sample_10[i,'exsit_city'] <- fill
  }else{
    count_sample_10[i,'exsit_city'] <- NA
  }
}

ggplot(count_sample_10,aes(x=taxa_list,y=Percentage_overall))+
  geom_bar(stat='identity',aes(fill=factor(Percentage_city)),position = "dodge",width=0.5)+
  geom_text(stat='identity', mapping = aes(y = Percentage_overall+0.01, label = exsit_city),hjust = 0,size=3.5) +
  expand_limits(y=c(0,0.08))+
  xlab('Taxa')+
  ylab('Percentage of samples (%)')+
  labs(fill="Percentage of cities (%)")+
  ggtitle('Rarest MetaSub pathogens (top 20)')+
  coord_flip()+
  theme_bw()+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_fill_manual(values = c("darkblue"))
# ggsave('metasub rarest pathogen 20.png',p)


#### 9 city pathogen abundance

bar_df_abundance <- abundance_bar

for (i in 1:dim(bar_df_abundance)[1]){
  bar_df_abundance[i,'pathogen_p'] <- sum(as.numeric(bar_df_abundance[i,pathogen_list]))
  bar_df_abundance[i,'non_pathogen_p'] <- 1 - bar_df_abundance[i,'pathogen_p']
  bar_df_abundance[i,'pathogen_n'] <- sum(bar_df_abundance[i,pathogen_list] > 0)
  bar_df_abundance[i,'non_pathogen_n'] <- 3660 - bar_df_abundance[i,'pathogen_n']
}

animal_pathogen <- pathogen[pathogen$pathogen == 1,'Species']
plant_pathogen <- pathogen[pathogen$pathogen == 2,'Species']
both_pathogen <- pathogen[pathogen$pathogen == 3,'Species']

for (i in 1:dim(bar_df_abundance)[1]){
  bar_df_abundance[i,'animal_pathogen_n'] <- sum(bar_df_abundance[i,animal_pathogen] > 0)
  bar_df_abundance[i,'animal_pathogen_p'] <- sum(as.numeric(bar_df_abundance[i,animal_pathogen]))
  bar_df_abundance[i,'plant_pathogen_n'] <- sum(bar_df_abundance[i,plant_pathogen] > 0)
  bar_df_abundance[i,'plant_pathogen_p'] <- sum(as.numeric(bar_df_abundance[i,plant_pathogen]))
  bar_df_abundance[i,'both_pathogen_n'] <- sum(bar_df_abundance[i,both_pathogen] > 0)
  bar_df_abundance[i,'both_pathogen_p'] <- sum(as.numeric(bar_df_abundance[i,both_pathogen]))
}

bar_df2 <- bar_df_abundance[,c('plant_pathogen_p','both_pathogen_p','animal_pathogen_p')]

# png("pathogen_abundance_separate_0307.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(3,0,0,0))

bp <- barplot(t(bar_df2*100), col=c("#619CFF","#00BA38","#F8766D"), 
              ylim = c(0,50),
              names.arg= paste0(country_name,' (',bar_df_abundance$pathogen_n,')') ,
              args.legend = list(x = "top", inset=c(-0.5,0)), las =2, cex.names=.75,
              ylab = "Pathogen relative abundance %",
              border = NA,)
# legend("top",inset = c(-0.32,-0), 
#        c('Animal Pathogen','Both Pathogen','Plant Pathogen'), 
#        fill = rev(c("#619CFF","#00BA38","#F8766D")) , bty = 1, cex = 1,
#        border = NA,n=3)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()


## pathogen abundance of each city on map

metadata <- metasub_data

for (i in 2: length(bar_df_abundance$X)){
  city <- bar_df_abundance$X[i]
  lat <- sum(as.numeric(metadata[metadata$city == city,'latitude']))/length(metadata[metadata$city == city,'latitude'])
  long <- sum(as.numeric(metadata[metadata$city == city,'longitude']))/length(metadata[metadata$city == city,'longitude'])
  bar_df_abundance[i,'lat'] <- lat
  bar_df_abundance[i,'long'] <- long
  bar_df_abundance[i,'city_origin_name'] <- country_name[which(city_name == city)]
}

ggplot()+
  borders("world",colour = "darkgrey",fill="grey")+
  ylim(-55,80)+
  xlim(-180,180)+
  ylab('Latitude')+
  xlab("Longitude")+
  geom_point(bar_df_abundance,mapping = aes(x=long,y=lat,col=pathogen_p),size=1.5) +
  geom_text_repel(bar_df_abundance,mapping = aes(x=long,y=lat,label=city_origin_name,col=pathogen_p),size=3.5,max.overlaps = 20)+
  scale_colour_gradient(low = 'blue', high = "red")+
  theme_bw() + 
  labs(color = "Pathogen\nabundance") +
  theme(panel.grid=element_blank(), panel.background = element_rect(fill="lightskyblue1"))

# England insets
ggplot()+
  borders("world",colour = "darkgrey",fill="grey")+
  ylim(40,60)+
  xlim(-10,5)+
  geom_point(bar_df_abundance,mapping = aes(x=long,y=lat,col=pathogen_p*100),size=2) +
  geom_text_repel(bar_df_abundance,mapping = aes(x=long,y=lat,label=city_origin_name,col=pathogen_p*100),size=3.5)+
  scale_colour_gradient(low = 'blue', high = "red")+
  theme_bw() + 
  labs(color = "Pathogen relative abundance (%)") +
  theme(panel.grid=element_blank(), panel.background = element_rect(fill="lightskyblue1"))

## correlation of pathogen size and abundance
ggplot(data = bar_df_abundance[-1,],aes(x=pathogen_n,y=(pathogen_p)*100))+
  geom_point(size=2)+
  ylab('Relative abundance of pathogen %')+
  xlab('Number of pathogen in city')+
  stat_smooth(method=lm,formula=y~x,se=FALSE)+
  geom_smooth (method = lm,linetype=1,se=FALSE,span=1)+
  theme_bw()+
  ggpubr::stat_cor(label.x = 170,label.y = 55,color = 'black')
```

### 2.4 Prediction result of MetaSUB using mGPS
The prediction result files (MetaSUB_prediction_result.csv) of MetaSUB used in following were produced by mGPS and have been provided in the Github. These files can be produced directly on mGPS interface.  

#### 2.4.1 Pre-process of MetaSUB input file
```r
complete_meta <- read.csv(file = "complete_metadata.csv", header = T)
taxa_abund <-read.csv("metasub_taxa_abundance.csv", header = T)
taxa_abund <- unique(taxa_abund)

# Merge RSA and metadata
metasub_data <- merge(complete_meta,taxa_abund,by.x="uuid",by.y="uuid")

# Remove control samples
control_samples <- c( which(metasub_data$city %in% c("control", "other_control","neg_control","other","pos_control")), which(metasub_data$control_type %in% c("ctrl cities","negative_control","positive_control"))) 
metasub_data <- droplevels(metasub_data[-c(control_samples), ])

# Re-label london boroughs 
metasub_data$city[metasub_data$city %in% c("kensington","islington")] <- "london" 
metasub_data <- droplevels(metasub_data)

# Remove sparse samples locations and dubiously labelled samples. 
small_cities <-  names(which(summary(metasub_data$city) < 8))
remove_samples <- which(metasub_data$city %in%  c("antarctica", small_cities))
metasub_data <- droplevels(metasub_data[-c(remove_samples), ])

# Correction of identified misslabelling of data 
metasub_data$latitude[metasub_data$city == "kyiv"] <- metasub_data$city_latitude[metasub_data$city == "kyiv"]
metasub_data$longitude[metasub_data$city == "kyiv"] <- metasub_data$city_longitude[metasub_data$city == "kyiv"]
metasub_data$continent[metasub_data$city == "porto"] <- "europe"

metasub_data[is.na(metasub_data$latitude),]$latitude <- metasub_data[is.na(metasub_data$latitude),]$city_latitude
metasub_data[is.na(metasub_data$longitude),]$longitude <- metasub_data[is.na(metasub_data$longitude),]$city_longitude

# correction to some incorrect city co-ords
metasub_data[metasub_data$city == "london",]$city_latitude <- 51.50853
metasub_data[metasub_data$city == "london",]$city_longitude <- -0.12574
write.csv(metasub_data,'metasub_data.csv',row.names=F)
```
#### 2.4.2 mGPS interface implementation by MetaSUB data
OPen mGPS_interface.r, click *RunApp*.  
Select `Prediction program` as `Build a new prediction model using mGPS`  
Select `Merged metadata and abundance data`  
`Input file`: metasub_data.csv  
`Enter the main locality level`: city  
`Enter the locality hierarchy`: continent,city,latitude,longitude  
`Column range of abundance data`: 43:3711  
`Subsets in feature elimination`: 50,100,200,300,500,1500  
Click `Start` and `Result Plot` tab  
When application done, click Output tab to download the result files: `Download optimal features in prediction model`: global_GITS.csv; `Download prediction data`: MetaSUB_prediction_result.csv; `DownloadModel`: MetaSub_model.Rda  

## 3 Characterizing the GITs used for MetaSUB biolocalization

### 3.1 The hallmarks of GITs
Used R to detect the hallmarks of GITs and the abundance of pathogenic GITs in MetaSUB.  
The variable in this part was connected with the last part.

```r
### 1. Complete GITs information (occupancy, relative abundance, taxonomic level)

metasub_git_file <- read.csv('global_GITS.csv',check.names = F)
git <- metasub_git_file

taxonomic_level1 <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
taxonomic_level2 <- c('kingdom','phylum','class','order','family','genus','species')

for (i in 1:dim(git)[1]){
  taxa <- git$taxa[i]
  git[i,'global_percentage'] <- metasub_abundance_file[1,taxa] # global average relative abundance of gits
  git[i,'Percentage_overall'] <- t_metasub_occupancy_file[taxa,'Percentage_overall'] # sample occupancy of gits
  git[i,'Percentage_city'] <- t_metasub_occupancy_file[taxa,'Percentage_city'] # city occupancy of gits
  
  for (j in 1:length(taxonomic_level1)){
    if (length(annotated_taxa[annotated_taxa$taxa==taxa,taxonomic_level1[j]]) != 0){
      git[i,taxonomic_level2[j]] <- annotated_taxa[annotated_taxa$taxa==taxa,taxonomic_level1[j]]
    }else{
      git[i,taxonomic_level2[j]] <- 'Non-annotated GITs'
    }
  }
}

sum(git$kingdom == 'Non-annotated GITs') # 15 GITs has not been annotated

git$phylum <- factor(git$phylum)
git$class <- factor(git$class)


### 2. occupancy distribution of GITs and non-GITs

# occupancy of GITs
git_sample <- git$Percentage_overall
git_city <- git$Percentage_city

# occupancy of non-GITs
git_list <- git$taxa

taxa_name_all <- rownames(t_metasub_occupancy_file)
git_position <- c(which(taxa_name_all %in% git_list))

non_git_sample <- t_metasub_occupancy_file[-c(git_position),'Percentage_overall'] # sample occupancy of non-GITs
non_git_city <- t_metasub_occupancy_file[-c(git_position),'Percentage_city'] # city occupancy of non-GITs

# create plot data form
plot_df <- rbind.data.frame(cbind(type = rep('git_sample',length(git_sample)),ubiquity = git_sample),
                            cbind(type = rep('git_city',length(git_city)),ubiquity = git_city),
                            cbind(type = rep('non-git_sample',length(non_git_sample)),ubiquity = non_git_sample),
                            cbind(type = rep('non-git_city',length(non_git_city)),ubiquity = non_git_city))

plot_df$ubiquity <- as.numeric(plot_df$ubiquity) 

ggplot(plot_df, aes(x = ubiquity, fill=type)) +
  geom_density(alpha = 0.3)+
  xlab("GITs ubiquity (%)")+
  theme_bw()+
  scale_fill_discrete(name="Sample type",
                      breaks=c("git_city","git_sample","non-git_city","non-git_sample"),
                      labels=c("GITs in sites","GITs in cities","non-GITs in sites","non-GITs in samples"))


### 3. occupancy versus importance value (informative)
ggplot(git)+
  geom_point(aes(x=Percentage_overall,y=Imp,color='red'))+
  geom_point(aes(x=Percentage_city,y=Imp,color='blue'))+
  xlab("GITs ubiquity (%)")+
  ylab('Importance')+
  labs(color = "Sample type") +
  theme_bw()+
  scale_color_manual(labels = c("Percentage of cities","Percentage of samples"), values = c("red", "blue")) 


### 4. occupancy versus log(abundance)

# The fitting formula is calculated by Excel: log10(relative abundance) on the x-axis and occupancy on the y-axis
x <- log10(git$global_percentage)
y <- (-4.8604)*(log10(git$global_percentage)^2)-10.555*log10(git$global_percentage)+92.935 

ggplot(git,aes(x=log10(global_percentage),y=Percentage_overall))+
  geom_point()+
  ylab("Sample occupancy")+
  xlab("log (mean global relative abundance of GITs)")+
  geom_line(mapping=aes(x=x,y=y),color = 'blue')+
  theme_bw()+
  # ggpubr::stat_cor(label.x = -15,label.y = 100,color = 'black')+
  theme(axis.text.x = element_blank())


### 5. relative abundance versus importance

git_abun <- git$global_percentage

git_position2 <- c(which(colnames(all_taxa) %in% git_list))
non_git_abun <- metasub_abundance_file[1,-c(git_position2)]

plot_df2 <- rbind.data.frame(cbind(type = rep('git',length(git_abun)),abundance = git_abun),
                             cbind(type = rep('non-git',length(non_git_abun)),abundance = non_git_abun))

plot_df2$abundance <- as.numeric(plot_df2$abundance)

p1 = ggplot(plot_df2, aes(x = log10(abundance*100), fill=type)) +
  geom_density(alpha = 0.3)+
  xlab("log (mean global relative abundance %)")+
  xlim(-9,2)+
  scale_fill_discrete(name="Sample type",
                      breaks=c("git","non-git"),
                      labels=c("GITs","non-GITs"))

p2 = ggplot(git,aes(x=log10(global_percentage*100),y=Imp))+
  geom_point()+
  ylab('Importance')+
  xlab("log (mean global relative abundance %)")+
  xlim(-9,2)

# combine two figures
p1 = p1 +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank(),
        legend.position="top")

p2 = p2+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x =element_blank(),
        axis.title.y =element_blank())%+replace% 
  theme(panel.background = element_rect(fill = NA))

g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# overlap to get two layer
pp <- c(subset(g1$layout, name == "panel", se = t:r))
g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], 
                     pp$t,pp$l, pp$b, pp$l)
grid.newpage()

# axis adjustment
ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths))
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
grid.newpage() 
grid.draw(g)


### 6. Relative abundance of GIT and non-GIT per city.

city_mGPS <- c("auckland","baltimore","barcelona","berlin","bogota","brisbane","denver","doha",
               "fairbanks","hamilton","hanoi","hong_kong","ilorin","kuala_lumpur","kyiv",
               "lisbon","london","marseille","minneapolis","naples","new_york_city","offa","oslo",
               "paris","porto","rio_de_janeiro","sacramento","san_francisco","santiago","sao_paulo",
               "sendai","seoul","singapore","sofia","stockholm","taipei","tokyo","vienna","yamaguchi","zurich")

bar_df_abundance <- metasub_abundance_file[c(city_mGPS),]
metasub_abundance_mGPS_file <- bar_df_abundance

global_proportion <- c('')
for (i in 1:dim(bar_df_abundance)[2]){
  count <- sum(bar_df_abundance[,i])/length(bar_df_abundance[,i])
  global_proportion <- c(global_proportion,as.numeric(count))
}

bar_df_abundance <- rbind.data.frame(global_proportion,bar_df_abundance)
rownames(bar_df_abundance) <- c('Overall',city_mGPS)

country_name <- c('Global')
for (country in city_mGPS){
  if ('_' %in% strsplit(country,'')[[1]]){
    country <- strsplit(country,'_')[[1]]
    for (a in 1:length(country)){
      country_a <- strsplit(country[a],'')[[1]]
      country_a[1] <- toupper(country_a[1])
      country_a <- paste(country_a, collapse='')
      country[a] <- country_a}
    country <- paste(country, collapse=' ')
  }else{
    country <- strsplit(country,'')[[1]]
    country[1] <- toupper(country[1])
    country <- paste(country, collapse='')
  }
  country_name <- c(country_name,country)
}

# calculate total relative abundance of GIT
for (i in 1:dim(bar_df_abundance)[1]){
  uniq_git <- c()
  uniq_n <- 0
  bar_df_abundance[i,'git_p'] <- sum(as.numeric(bar_df_abundance[i,git$taxa]))
  bar_df_abundance[i,'non_git_p'] <- 1 - bar_df_abundance[i,'git_p']
  bar_df_abundance[i,'git_n'] <- sum(bar_df_abundance[i,git$taxa] > 0)

  if (i == 1){
    next
  }
  # find city-unique taxa
  for (j in 1:dim(bar_df_abundance)[2]){
    taxa_git <- colnames(bar_df_abundance)[j]
    if (sum(as.numeric(bar_df_abundance[-i,][,j])) == 0){
      uniq_git <- c(uniq_git,taxa_git)
      uniq_n <- uniq_n + 1
    }
  }
  bar_df_abundance[i,'uniq_git_n'] <- uniq_n
  bar_df_abundance[i,'uniq_git'] <- paste(uniq_git,sep = ',')
}

bar_df_abundance$X <- rownames(bar_df_abundance)

## plot relative abundance of GIT and non-GIT per city
bar_df0 <- bar_df_abundance[,c('git_p','non_git_p')]

par(xpd = T, mar = par()$mar + c(1,0,0,4))
bp <- barplot(t(bar_df0*100), col=c("royalblue3","slategray1"), names.arg= paste0(country_name,' (',bar_df_abundance$git_n,')' ),
              args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.6,ylab = "taxa abundance %")
legend("topright",inset = c(-0.1,0.4), c('Non-Gits','Gits'), fill = rev(c("royalblue3","slategray1")) , bty = 1, cex = 0.8)
par(mar=c(5, 4, 4, 2) + 0.1)

# png("MetaSub_git_abundance.png", width = 13,height = 8, units = 'in', res = 600)
par(xpd = T, mar = par()$mar + c(6,0,0,4), mgp = c(0,0.7,0), las=2)

bp <- barplot(t(bar_df0*100), col=c("royalblue3","slategray1"),
              names.arg=paste0(country_name,' (',bar_df_abundance$git_n,')' ),
              args.legend = list(x = "topright", inset=c(-0.5,0)), las =2,
              cex.names=.6,ylab = "", axisnames = F,axes = F, space =0)

axis(side =2, pos = 0)
mtext(text = paste0(country_name,' (',bar_df_abundance$git_n,')' ), side = 1, at = bp, line = 0, padj = 1, cex = 0.8)
title(ylab="Taxa relative abundance %", mgp=c(0,0,0),cex.lab=1.2)
legend("topright",inset = c(-0.08,0.4), c('Non-GITs','GITs'), fill = rev(c("royalblue3","slategray1")) , bty = 1, cex = 1)

par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()

# calculate city average relative abundance of each taxa
for (i in 1:dim(git)[1]){
  taxa <- git$taxa[i]
  git[i,'global_percentage'] <- bar_df_abundance[1,taxa] # global average relative abundance of gits
}

git$global_percentage <- as.numeric(git$global_percentage)


## GITs on phylum level

# calculate phylum relative abundance
phylum_df <- data.frame(row.names = c(levels(factor(bar_df_abundance$X))))

for (i in 1:length(levels(factor(bar_df_abundance$X)))){
  
  city <- levels(factor(bar_df_abundance$X))[i]
  
  for (j in 1:length(levels(git$phylum))){
    phylum <- levels(git$phylum)[j]
    taxa_list <- git[git$phylum == phylum,][,'taxa']
    # Calculate the sum of relative abundance values of taxa in the same phylum
    percentage <- sum(as.numeric(bar_df_abundance[bar_df_abundance$X==city,taxa_list]))
    phylum_df[city,phylum] <- percentage
  }
  phylum_df[city,'others'] <- 1 - sum(phylum_df[city,levels(git$phylum)])
}

po <- which(row.names(phylum_df) == 'Overall')
phylum_df_all <- phylum_df[po,]
phylum_df <- phylum_df[-po,]
phylum_df <- rbind.data.frame(phylum_df_all,phylum_df)

# calculate phylum amount per city
phylum_n <- c()
for (i in 1:length(row.names(phylum_df))){
  phylum_n <- c(phylum_n,sum(phylum_df[i,-c(9,11)] > 0))
}

phylum_df2 <- cbind.data.frame(phylum_df,'phylum_n'= phylum_n)

palette <- colorRampPalette(brewer.pal(12, "Paired"))(11)[1:10]

bar_df2 <- phylum_df[,1:10]

# plot figure
par(xpd = T, mai=c(1.9,1,0.5,0.5))
bp <- barplot(t(bar_df2*100), col=palette, names.arg= paste0(row.names(bar_df2),' (',phylum_n,')' ),args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.6,ylab = "taxa abundance %")
legend("bottom",inset = c(0,-0.95), rev(colnames(bar_df2)), fill = rev(palette) , bty = 1, cex = 0.8,ncol=4)
par(mar=c(5, 4, 4, 2) + 0.1)

# png("git_phylum_abundance.png", width = 13,height = 8, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(6,0,0,0))
bp <- barplot(t(bar_df2*100), col=palette,
              names.arg= paste0(country_name,' (',phylum_n,')' ),
              args.legend = list(x = "bottom", inset=c(-0.5,0)), las =2,
              cex.names=.7,ylab = "", axisnames = F, axes = T, space =0,
              ylim = c(0,100))
mtext(text = paste0(country_name,' (',phylum_n,')' ), side = 1, at = bp, line = 0.5, padj = 1, cex = 0.8)
title(ylab="GITs relative abundance % (Phylum level)",mgp=c(2,0,0),cex.lab=1.2)
legend("bottom",inset = c(0,-0.38), rev(colnames(bar_df2)),  fill = rev(palette) , bty = 1, cex = 1,ncol=5)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()


## GITs on class level

git$class <- as.character(git$class)
git[137,'class'] <- 'Kinetoplastea'
git$class <- factor(git$class)

level_o <- levels(factor(bar_df_abundance$X))
level_n <- c('Overall',level_o[1:23],level_o[25:41])

# calculate class relative abundance
class_df <- data.frame(row.names = level_n)
for (i in 1:length(level_n)){
  
  city <- level_n[i]
  
  for (j in 1:length(levels(git$class))){
    class <- levels(git$class)[j]
    taxa_list <- git[git$class == class,][,'taxa']
    percentage <- sum(as.numeric(bar_df_abundance[bar_df_abundance$X==city,taxa_list]))
    class_df[city,class] <- percentage
  }
  class_df[city,'others'] <- 1 - sum(class_df[city,levels(git$class)])
}

# number of class in each city
class_n <- c()
for (i in 1:length(row.names(class_df))){
  class_n <- c(class_n,sum(class_df[i,-c(12,17)] > 0))
}

class_df2 <- cbind.data.frame(class_df,'class_n'= class_n)

for (i in 1:dim(class_df)[1]){
  po <- which(class_df[i,]==0)
  blank <- colnames(class_df)[po]
  blank <- paste(blank,collapse = ',')
  class_df2[i,'lack'] <- blank
}

# plot

palette <- colorRampPalette(brewer.pal(12, "Paired"))(17)[1:16]
palette[11] <- "#FE982C"
# palette[17] <- "#F0EB99"

# png("git_class_abundance.png", width = 13,height = 8, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(6,0,0,0))
bp <- barplot(t(bar_df2*100), col=palette,
              names.arg= paste0(country_name,' (',class_n,')' ),
              args.legend = list(x = "bottom", inset=c(-0.5,0)), las =2,
              cex.names=.7,ylab = "", axisnames = F,axes = T, space =0)
# axis(side =2, pos = -1)
mtext(text = paste0(country_name,' (',class_n,')' ), side = 1, at = bp, line = 0.5, padj = 1, cex = 0.8)
title(ylab="GITs relative abundance % (Class level)", mgp=c(2,0,0),cex.lab=1.2)
legend("bottom",inset = c(0,-0.43), rev(colnames(bar_df2)),  fill = rev(palette) , bty = 1, cex = 0.8,ncol=6)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()
```
### 3.2 The abundance of pathogenic GITs
Connected with last part of R script

```r
### 7. Relative abundance of pathogenic GITs

# pathogen git analysis

pathogen_position <- c(which(colnames(bar_df_abundance) %in% pathogen$Species))
git_position <- c(which(colnames(bar_df_abundance) %in% git$taxa))

git_pathogen_position <- c()
for (i in git_position){
  if (i %in% pathogen_position){
    git_pathogen_position <- c(i,git_pathogen_position)
  }
}

for (i in 1:dim(bar_df_abundance)[1]){
  bar_df_abundance[i,'git_pathogen'] <- sum(as.numeric(bar_df_abundance[i,git_pathogen_position]))
  bar_df_abundance[i,'git_non-pathogen'] <- sum(as.numeric(bar_df_abundance[i,git$taxa])) - bar_df_abundance[i,'git_pathogen']
  bar_df_abundance[i,'git_n'] <- sum(bar_df_abundance[i,git$taxa] > 0)
  bar_df_abundance[i,'git_pathogen_n'] <- sum(bar_df_abundance[i,git_pathogen_position] > 0)
}

animal_pathogen_position <- c(which(colnames(bar_df_abundance) %in% pathogen[pathogen$pathogen == 1,'Species']))
plant_pathogen_position <- c(which(colnames(bar_df_abundance) %in% pathogen[pathogen$pathogen == 2,'Species']))
both_pathogen_position <- c(which(colnames(bar_df_abundance) %in% pathogen[pathogen$pathogen == 3,'Species']))

for (i in 1:dim(bar_df_abundance)[1]){
  
  git_pathogen_a <- c()
  git_pathogen_p <- c()
  git_pathogen_b <- c()
  for (j in git_position){
    if (j %in% animal_pathogen_position){
      git_pathogen_a <- c(j,git_pathogen_a)
    }else if (j %in% plant_pathogen_position) {
      git_pathogen_p <- c(j,git_pathogen_p)
    }else if (j %in% both_pathogen_position) {
      git_pathogen_b <- c(j,git_pathogen_b)
    }
  }
  
  bar_df_abundance[i,'git_animal_pathogen_n'] <- sum(bar_df_abundance[i,git_pathogen_a] > 0)
  bar_df_abundance[i,'git_animal_pathogen_p'] <- sum(as.numeric(bar_df_abundance[i,git_pathogen_a]))
  bar_df_abundance[i,'git_plant_pathogen_n'] <- sum(bar_df_abundance[i,git_pathogen_p] > 0)
  bar_df_abundance[i,'git_plant_pathogen_p'] <- sum(as.numeric(bar_df_abundance[i,git_pathogen_p]))
  bar_df_abundance[i,'git_both_pathogen_n'] <- sum(bar_df_abundance[i,git_pathogen_b] > 0)
  bar_df_abundance[i,'git_both_pathogen_p'] <- sum(as.numeric(bar_df_abundance[i,git_pathogen_b]))
}

# Plot the RSA of git in each city

bar_df_abundance2 <- bar_df_abundance[order(bar_df_abundance$continent,decreasing=T),][-41,] # sample_size
bar_df_abundance2 <- rbind.data.frame(bar_df_abundance[order(bar_df_abundance$sample_size,decreasing=T),][41,],bar_df_abundance2)
bar_df_abundance3 <- bar_df_abundance2[,c('X','continent')]
bar_df2 <- bar_df_abundance2[,c('git_plant_pathogen_p','git_both_pathogen_p','git_animal_pathogen_p','git_non-pathogen')]

country_name2 <- c()

for (country in bar_df_abundance2$X){
  if ('_' %in% strsplit(country,'')[[1]]){
    country <- strsplit(country,'_')[[1]]
    for (a in 1:length(country)){
      country_a <- strsplit(country[a],'')[[1]]
      country_a[1] <- toupper(country_a[1])
      country_a <- paste(country_a, collapse='')
      country[a] <- country_a}
    country <- paste(country, collapse=' ')
  }else{
    country <- strsplit(country,'')[[1]]
    country[1] <- toupper(country[1])
    country <- paste(country, collapse='')
  }
  country_name2 <- c(country_name2,country)
}

# png("city_pathogen_git.png", width = 13,height = 8, units = 'in', res = 600)
par(xpd = T, mar = par()$mar + c(2.5,0,0,1))
bp <- barplot(t(bar_df2*100), col=c("#619CFF","#00BA38","#F8766D","slategray2"), 
              ylim = c(0,100), xlim = c(0,60),
              names.arg= paste0(country_name2,' (',bar_df_abundance2$git_pathogen_n,'/',bar_df_abundance2$git_n,')' ) ,
              args.legend = list(x = "top", inset=c(-0.5,0)), las =2, cex.names=.75,
              ylab = "Taxa relative abundance %",
              border = NA,)
legend("top",inset = c(-0.2,-0), 
       c('Non Pathogen GITs','Animal Pathogen GITs','Both Pathogen GITs','Plant Pathogen GITs'), 
       fill = rev(c("#619CFF","#00BA38","#F8766D","slategray2")) , bty = 'n', cex = 1,
       border = NA,n=3)

par(xpd = T, new=T)

plot(bp,y=(bar_df_abundance2$pathogen_p)*100,type = "b", col="black",axes = F,
     ylim = c(0,100),xlim = c(0,60),xlab = "",ylab = "",pch=19,bg="black",cex=0.8)
legend(35,104,inset = c(-0.35,0),legend = c("Pathogen abundance in all taxa"),pch = 19,col = "black",bty = "n",cex = 1,
       pt.cex = 0.8,lty = 1)

par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()

# Correlation between city population size and pathogen's total relative abundance

# population in each city
for (i in 2:(dim(bar_df_abundance)[1])){
  city <- bar_df_abundance$X[i]
  sample_size <- sum(metasub_data[,'city'] == city)
  bar_df_abundance[i,'sample_size'] <- sample_size
  for (j in 1:(dim(metasub_data)[1])){
    city_meta <- metasub_data[j,'city']
    if (city_meta == city && !is.na(metasub_data[j, "city_total_population"])){
      bar_df_abundance[i,'population'] <- metasub_data[j,'city_total_population']
      break
    }
  }
  for (j in 1:(dim(metasub_data)[1])){
    city_meta <- metasub_data[j,'city']
    if (city_meta == city && !is.na(metasub_data[j, "continent"])){
      bar_df_abundance[i,'continent'] <- metasub_data[j,'continent']
      break
    }
  }
}

bar_df_sub1 <- bar_df_abundance[,c('population','git_pathogen')]
colnames(bar_df_sub1) <- c('population','abundance')
bar_df_sub1[,'Pathogen'] <- 'Gits(200)'
bar_df_sub2 <- bar_df_abundance[,c('population','pathogen_p')]
colnames(bar_df_sub2) <- c('population','abundance')
bar_df_sub2[,'Pathogen'] <- 'All_taxa'
bar_df3 <- rbind.data.frame(bar_df_sub1,bar_df_sub2)

library(ggpmisc)

bar_df3$Pathogen<-factor(bar_df3$Pathogen,
                         levels = c('All_taxa','Gits(200)'),
                         labels = c("Pathogen in all taxa","Pathogen in GITs"))

options(scipen = 3)

ggplot(data = bar_df3,aes(x=population,y=abundance*100,fill=Pathogen))+
  geom_point(size=2)+
  ylab('Pathogen relative abundance %')+
  xlab('City population')+
  stat_smooth(method=lm,formula=y~x,se=FALSE)+
  geom_smooth (method = lm,linetype=1,se=FALSE,span=1)+
  ggpubr::stat_cor(aes(color = Pathogen), label.x = 7000000,label.y = 45,color = 'black')+
  facet_wrap(~Pathogen) +
  theme_bw()+
  theme(legend.position="none")
```

## 4 Distinguishing local and non-local taxa using mGPS

### 4.1 Mobility among mGPS predicted local and non-local taxa
In "MetaSUB_prediction_result.csv", according to the predicted source coordinates of each sample, those whose predicted source sites were in the same continent as the sampling sites were labeled as local samples. For the samples whose predicted source sites were not in the same continent as the sampling sites, they were labeled as "sampling site-predicted site".   
In the following, demonstrated that mGPS localizes the non-local samples to their local regions and that local and non-local taxa exhibit different abundance curves.

```R
# differentiate local and non-local 
library(ggplot2)
library(ggsignif)
library(forcats)

prediction_result <- read.csv('MetaSUB_prediction_result.csv')
git <- read.csv('global_GITS.csv')

# Find taxa with RSA as 0 according to the global RSA of each taxon in previous analysis.
zero_taxa <- c('Erwinia.phage.vB_EamM_Phobos',
               'Listeria.phage.A500',
               'Pseudomonas.phage.F10',
               'Staphylococcus.phage.B236',
               'Staphylococcus.phage.DW2',
               'Staphylococcus.phage.ROSA',
               'Staphylococcus.virus.29',
               'Staphylococcus.virus.X2',
               'Thermoanaerobacter.kivui')

taxa_name <- read.csv("metasub_taxa_abundance.csv",check.names = F) 
taxa_name <- as.data.frame(colnames(taxa_name)[2:3670])
colnames(taxa_name) <- 'uuid'
taxa_name$make_id <- make.names(taxa_name$uuid)

non_gits <- colnames(prediction_result)[43:3711] 
non_gits <- non_gits[-which(non_gits %in% zero_taxa)]
non_gits <- non_gits[-which(non_gits %in% git$taxa)]

# Function to do the test
migration_test <- function(sampling_c, predicted_c, migration){
  
  continent_1 <- prediction_result[prediction_result$continent == sampling_c,]
  continent_1_native <- continent_1[continent_1$migration=='local',]
  continent_2 <- prediction_result[prediction_result$continent == predicted_c,]
  continent_2_native <- continent_2[continent_2$migration=='local',]
  continent_1_to_2 <- continent_1[continent_1$migration==migration,]
  
  continent_1_native$migration <- 'continent_1'
  continent_2_native$migration <- 'continent_2'
  continent_1_to_2$migration <- '1-2'
  total_continent <- rbind.data.frame(continent_1_native,continent_2_native,continent_1_to_2)
  
  # Test RSA of top 15 informative GITs
  continent_git <- total_continent[,git$taxa[1:15]]
  continent_git$migration <- total_continent$migration
  reshape2::melt(continent_git,id.vars="migration") -> dfa1_1
  
  for (i in 1:length(dfa1_1$variable)){
    taxa_new <- dfa1_1$variable[i]
    dfa1_1[i,'variable2'] <- taxa_name[taxa_name$make_id == taxa_new,'uuid']
  }
  
  dfa1_1$variable2 <- factor(dfa1_1$variable2)
  dfa1_1$variable2 <- fct_inorder(dfa1_1$variable2)
  
  p = ggplot(dfa1_1, aes(x=migration,y=value,fill=migration)) + 
    geom_boxplot() +
    ylim(0,2)+
    ylab('Relative abundance')+
    ggtitle(paste0('Sampling site: ',sampling_c,'\nPredicted source: ',predicted_c))+
    guides(fill=FALSE) +
    facet_wrap(~variable2,ncol = 5) +
    geom_signif(comparisons = list(c("1-2", "continent_1"),
                                   c("1-2","continent_2"),
                                   c("continent_1","continent_2")),
                step_increase = 0.3,
                map_signif_level = T,
                test = wilcox.test)+
    theme_bw()+
    theme(axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 90))+
    scale_x_discrete("Sample source", labels = c("1-2" = "Non-local","continent_1" = "\nSampling site", "continent_2" = "Predicted site"))
  ggsave(paste0(sampling_c,' to ',predicted_c,'.png'),p)
  
  # Control one: Mix non-local and local samples in sampling site
  continent_random <- rbind.data.frame(continent_1_native,continent_1_to_2)
  continent_random2 <- continent_random[sample(nrow(continent_random), ), ]
  continent_random3 <- continent_random2[sample(nrow(continent_random2), 20), ]
  
  continent_random3$migration <- 'Random sample'
  for (i in 1:length(continent_random3$uuid)){
    uuid_new <- paste0(continent_random3$uuid[i],'-2')
    continent_random3[i,'uuid'] <- uuid_new
  }
  
  total_continent <- rbind.data.frame(continent_1_native,continent_2_native,continent_random3)
  continent_git <- total_continent[,git$taxa[1:15]]
  continent_git$migration <- total_continent$migration
  reshape2::melt(continent_git,id.vars="migration") -> dfa1_2
  
  dfa1_2$migration <- factor(dfa1_2$migration)
  dfa1_2$migration<- fct_relevel(dfa1_2$migration,c("Random sample","continent_1","continent_2"))  
  
  for (i in 1:length(dfa1_2$variable)){
    taxa_new <- dfa1_2$variable[i]
    dfa1_2[i,'variable2'] <- taxa_name[taxa_name$make_id == taxa_new,'uuid']
  }
  
  dfa1_2$variable2 <- factor(dfa1_2$variable2)
  dfa1_2$variable2 <- fct_inorder(dfa1_2$variable2)
  
  p = ggplot(dfa1_2, aes(x=migration,y=value,fill=migration)) + 
    geom_boxplot() +
    ylim(0,2)+
    ylab('Relative abundance')+
    ggtitle(paste0('Control: Random sample\nSampling site: ',sampling_c,'\nPredicted source: ',predicted_c))+
    guides(fill=FALSE) +
    facet_wrap(~variable2,ncol = 5) +
    geom_signif(comparisons = list(c("Random sample", "continent_1"),
                                   c("Random sample","continent_2"),
                                   c("continent_1","continent_2")),
                step_increase = 0.3,
                map_signif_level = T,
                test = wilcox.test)+
    theme_bw()+
    theme(axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 90))+
    scale_x_discrete("Sample source", labels = c("1-2" = "Non-native","continent_1" = "\nSampling site", "continent_2" = "Predicted site"))
  ggsave(paste0(sampling_c,' to ',predicted_c,'_control1.png'),p)
  
  # Control two: non-GITs detection (repeat 1000 times)
  repeat_random_15non_git <- c()
  
  for (re in 1:1000){
    continent_git <- total_continent[,sample(non_gits, 15)]
    continent_git$migration <- total_continent$migration
    
    flag <- 0
    for (t in 1:15){
      df_taxa <- continent_git[,c(t,16)]
      con_1 <- df_taxa[df_taxa$migration == 'continent_1',1]
      con_2 <- df_taxa[df_taxa$migration == 'continent_2',1]
      con_1_2 <- df_taxa[df_taxa$migration == '1-2',1]
      compare_1_2 <- wilcox.test(con_1,con_2)
      compare_n_1 <- wilcox.test(con_1,con_1_2)
      compare_n_2 <- wilcox.test(con_2,con_1_2)
      p_1_2 <- compare_1_2$p.value
      p_n_1 <- compare_n_1$p.value
      p_n_2 <- compare_n_2$p.value
      
      if (!is.na(p_1_2) & !is.na(p_n_1) & !is.na(p_n_2)){
        if (p_1_2 < 0.05){
          if (p_n_1 < 0.05){
            if (p_n_1 < p_n_2){
              flag <- flag + 1
              print(t)
            }
          }
        }
      }
    }
    repeat_random_15non_git <- c(repeat_random_15non_git,flag) 
  }
  
  mean_control <- mean(repeat_random_15non_git)
  median_control <- median(repeat_random_15non_git)
  
  return(list(mean_control,median_control))
}

# Sampling site: Europe Predicted original site: sub_saharan_africa
E_SSA <- migration_test('europe','sub_saharan_africa','E-SSA')

# Sampling site: Europe Predicted original site: North America
E_NA <- migration_test('europe','north_america','E-NA')

# Sampling site: Oceania Predicted original site: East Asia
O_EA <- migration_test('oceania','east_asia','O-EA')
```

### 4.2 Rank abundance curve of local and non-local samples differentiated by mGPS
```r
library(BiodiversityR)
library(ggplot2)
library(forcats)
library(ggrepel)
library(scales)

metasub_data <- read.csv('MetaSUB_prediction_result.csv')
metasub_data$city <- factor(metasub_data$city)

taxa_list <- colnames(metasub_data)[43:3711]

bar_df_abundance <- data.frame(row.names = c(levels(metasub_data$city),paste0(levels(metasub_data$city),'_migrate')))

for (i in 1:length(levels(metasub_data$city))){
  for (j in 1: length(taxa_list)){
    city <- levels(metasub_data$city)[i]
    taxa <- taxa_list[j]  
    city_subset <- subset(metasub_data[metasub_data$city == city,])
    origin_same <- subset(city_subset[city_subset$city == city_subset$cityPred,])
    origin_diff <- subset(city_subset[city_subset$city != city_subset$cityPred,])
    
    abun_same <- sum(origin_same[,taxa])/length(origin_same[,taxa])
    abun_diff <- sum(origin_diff[,taxa])/length(origin_diff[,taxa])
    
    bar_df_abundance[city,taxa] <- abun_same
    
    migrate <- paste0(city,'_migrate')
    bar_df_abundance[migrate,taxa] <- abun_diff
  }
}

## Calculate RAC for local and non-local samples on global level

europe_city <- levels(factor(as.character(metasub_data[,'city'])))
europe_m <- paste0(europe_city,'_migrate')
europe_subset <- bar_df_abundance[c(europe_city,europe_m),]

# yamaguchi has no local samples
if ('yamaguchi' %in% europe_city){
  europe_city <- europe_city[europe_city != 'yamaguchi']
}

europe_avg <- colMeans(bar_df_abundance[europe_city,])
europe_m_avg <- colMeans(bar_df_abundance[europe_m,])
  
europe_subset2 <- rbind.data.frame(europe_avg,europe_m_avg)
rownames(europe_subset2) <- c('native','non-native')

rank_metasub <- data.frame()
for (i in rownames(europe_subset2)){
  rank_metasub_i <- data.frame(rankabundance(subset(europe_subset2,rownames(europe_subset2) == i),digits = 6))[1:2] # digits = 6/10
  rank_metasub_i$sample <- paste0(i,' (',sum(rank_metasub_i$abundance != 0),')')
  rank_metasub_i$number <- sum(rank_metasub_i$abundance != 0)
  rank_metasub <- rbind(rank_metasub,rank_metasub_i)
}

# Calculate the proportion of local and non-local samples
correctly_pred <-  round( mean(metasub_data[,"cityPred"]==metasub_data[,"city"]) * 100 ,1)
incorrectly_pred <- round( (100 - correctly_pred) ,1)
label_o <- paste0('Local: ',correctly_pred,'%\n','Non-local: ',incorrectly_pred,'%')
rank_metasub3[1,'text2'] <- label_o
rank_metasub3[,'Migration'] <- rank_metasub3[,'sample']

# k-s test of local and non-local samples
native_s <- rank_metasub3[rank_metasub3$sample == levels(factor(rank_metasub3$sample))[1],][,'abundance']
non_native_s <- rank_metasub3[rank_metasub3$sample == levels(factor(rank_metasub3$sample))[2],][,'abundance']

# Record the p value of ks test into file 
ks_result <- ks.test(native_s,non_native_s)
ks_p <- ks_result$p.value

# Label the p-value of ks.test
rank_metasub3[1,'text3'] <- paste0("italic('p-value')~","`: " ,ks_p,"`")

ggplot(rank_metasub3,aes(rank,log(abundance,10)))+ 
    geom_line(aes(linetype=Migration))+
    scale_linetype_manual(values=c( "solid","longdash"))+
    geom_text(rank_metasub3,mapping=aes(x=2500,y=-1.5,label=text2),nudge_x=0.1,nudge_y=0.1,size =3.2, hjust = 1,colour='black')+
    geom_text_repel(rank_metasub3,mapping=aes(x=2400,y=-2.15,label=text3),size =3.2,colour='black',parse = TRUE,force=0,hjust  =  0)+
    scale_size_manual(values=c(2, 2))+
    ggtitle('RAC for Global')+
    labs(x = 'Microbiome abundance rank', y = 'Relative abundance (log)', color = NULL)+
    theme(legend.position="none",panel.grid = element_blank(),panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill='transparent')) 

## Calculate RAC for local and non-local samples on each continent

# Divide local and non-local in each continent
europe_subset2 <- data.frame()
for (c in levels(factor(metasub_data$continent))){
  europe_city <- levels(factor(as.character(metasub_data[metasub_data$continent == c,][,'city'])))
  europe_m <- paste0(europe_city,'_migrate')
  europe_subset <- bar_df_abundance[c(europe_city,europe_m),]
  
  if ('yamaguchi' %in% europe_city){
    europe_city <- europe_city[europe_city != 'yamaguchi']
  }
  europe_avg <- colMeans(bar_df_abundance[europe_city,])
  europe_m_avg <- colMeans(bar_df_abundance[europe_m,])
  
  if (dim(europe_subset2)[1] == 0 ){
    europe_subset2 <- rbind.data.frame(europe_avg,europe_m_avg)
    colnames(europe_subset2) <- colnames(bar_df_abundance)
    rownames(europe_subset2) <- c(c,paste0(c,':migrate'))
    rowname_store <- rownames(europe_subset2)
  }else{
    europe_subset2 <- rbind.data.frame(europe_subset2,europe_avg,europe_m_avg)
    rownames(europe_subset2) <- c(rowname_store,c,paste0(c,':migrate'))
    rowname_store <- rownames(europe_subset2)
  }
}

# Calculation of RSA for each continent
rank_metasub <- data.frame()
for (i in rownames(europe_subset2)){
  rank_metasub_i <- data.frame(rankabundance(subset(europe_subset2,rownames(europe_subset2) == i),digits = 6))[1:2] # digits = 6/10
  rank_metasub_i$sample <- paste0(i,' (',sum(rank_metasub_i$abundance != 0),')')
  rank_metasub_i$number <- sum(rank_metasub_i$abundance != 0)
  rank_metasub <- rbind(rank_metasub,rank_metasub_i)
  sum(rank_metasub_i$abundance != 0)
}

rank_metasub3 <- subset(rank_metasub,abundance != 0)

# Label the continent and migration situation of each sample
for (i in 1:length(rank_metasub3$sample)){
  sam <- rank_metasub3$sample[i]
  country <- strsplit(sam,':')[[1]][1]
  mig <- strsplit(strsplit(sam,':')[[1]][2],' ',fixed = T)[[1]][1]
  if (!is.na(mig)){
    rank_metasub3[i,'Migration'] <- 'non-native'
    
    country <- strsplit(sam,':migrat')[[1]][1]
    if ('_' %in% strsplit(country,'')[[1]]){
      country <- strsplit(country,'_')[[1]]
      
      for (a in 1:length(country)){
        country_a <- strsplit(country[a],'')[[1]]
        country_a[1] <- toupper(country_a[1])
        country_a <- paste(country_a, collapse='')
        country[a] <- country_a}
      country <- paste(country, collapse=' ')
    }else{
      country <- strsplit(country,'')[[1]]
      country[1] <- toupper(country[1])
      country <- paste(country, collapse='')
    }
    rank_metasub3[i,'country'] <- strsplit(sam,':migrat')[[1]][1]
    rank_metasub3[i,'country2'] <- country
  }else{
    rank_metasub3[i,'Migration'] <- 'native'
    country <- strsplit(sam,' ')[[1]][1]
    if ('_' %in% strsplit(country,'')[[1]]){
      country <- strsplit(country,'_')[[1]]
      for (a in 1:length(country)){
        country_a <- strsplit(country[a],'')[[1]]
        country_a[1] <- toupper(country_a[1])
        country_a <- paste(country_a, collapse='')
        country[a] <- country_a}
      country <- paste(country, collapse=' ')
    }else{
      country <- strsplit(country,'')[[1]]
      country[1] <- toupper(country[1])
      country <- paste(country, collapse='')
    }
    rank_metasub3[i,'country'] <- strsplit(sam,' ')[[1]][1]
    rank_metasub3[i,'country2'] <- country
  }
}

# Calculate the proportion of local and non-local samples
city_mig_p  <- data.frame(row.names = levels(fct_inorder(rank_metasub3$country)))
for (city in levels(fct_inorder(rank_metasub3$country))){
  correctly_pred <-  round( mean(metasub_data[metasub_data$continent == city,"cityPred"]== 
                                   metasub_data[metasub_data$continent == city,"city"]) * 100 ,1)
  incorrectly_pred <- round( (100 - correctly_pred) ,1)
  city_mig_p[city,'native'] <- correctly_pred
  city_mig_p[city,'non-native'] <- incorrectly_pred
}

# K-S test
psign <- data.frame(row.names = levels(factor(rank_metasub3$country2)))
for (city_ks in levels(factor(rank_metasub3$country2))){
  city_sample <- rank_metasub3[rank_metasub3$country2 == city_ks,]
  native_s <- city_sample[city_sample$sample == levels(factor(city_sample$sample))[1],][,'abundance']
  non_native_s <- city_sample[city_sample$sample == levels(factor(city_sample$sample))[2],][,'abundance']
  
  ks_result <- ks.test(native_s,non_native_s)
  ks_p <- ks_result$p.value
  psign[city_ks,'p.value'] <- ks_p
}

for (o in 1:dim(city_mig_p)[1]){
  label_o <- paste0('Local: ',city_mig_p[o,'native'],'%\n','Non-local: ',city_mig_p[o,'non-native'],'%')
  
  for (m in 1:length(rank_metasub3$country)){
    if (rank_metasub3$country[m] == rownames(city_mig_p)[o]){
      p <- psign[rank_metasub3$country2,'p.value']
      rank_metasub3[m,'text2'] <- label_o
      rank_metasub3[m,'text3'] <- paste0("italic('p-value')~","`: " ,p,"`")
      break
    }
  }
}

color <- hue_pal()(7)
color2 <- c(color[5],color[1],color[6],color[2],color[4],color[3],color[7])

ggplot(rank_metasub3,aes(rank,log(abundance,10),color = country2))+ 
  geom_line(aes(linetype=Migration,color = country2))+
  scale_linetype_manual(values=c( "solid","longdash"))+
  scale_colour_manual(values = color2)+
  facet_wrap(~country2) +
  xlim(0,2500)+
  ylim(-6,-0.3)+
  geom_text(rank_metasub3,mapping=aes(x=rep(2490,length(country2)),y=rep(-1.5,length(country2)),label=text2),nudge_x=0.1,nudge_y=0.1,size =3.2, hjust = 1,colour='black')+
  geom_text_repel(rank_metasub3,mapping=aes(x=rep(2500,length(country2)),y=rep(-2.7,length(country2)),label=text3),size =3.2,colour='black',parse = TRUE,force=0,hjust  =  0)+
  scale_size_manual(values=c(2, 2))+
  labs(x = 'Abundance rank', y = 'Relative abundance (log)', color = NULL)+
  theme(legend.position="none",panel.grid = element_blank(),panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill='transparent')) 


## Calculate RAC for local and non-local samples in each city

color <- hue_pal()(7)
color <- c(color[5],color[1],color[6],color[2],color[4],color[3],color[7])

# Using loop to analyze cities in each continent
for (num in 1:length(levels(factor(metasub_data$continent)))){
  c <- levels(factor(metasub_data$continent))[num]
  europe_city <- levels(factor(as.character(metasub_data[metasub_data$continent == c,][,'city'])))
  europe_m <- paste0(europe_city,'_migrate')
  europe_subset <- bar_df_abundance[c(europe_city,europe_m),]

  # Calculate RSA for each city
  rank_metasub <- data.frame()
  for (i in rownames(europe_subset)){
    rank_metasub_i <- data.frame(rankabundance(subset(europe_subset,rownames(europe_subset) == i),digits = 6))[1:2] # digits = 6/10
    rank_metasub_i$sample <- paste0(i,' (',sum(rank_metasub_i$abundance != 0),')')
    rank_metasub_i$number <- sum(rank_metasub_i$abundance != 0)
    rank_metasub <- rbind(rank_metasub,rank_metasub_i)
  }
  
  rank_metasub <- subset(rank_metasub,abundance != 0)
  rank_metasub2 <- rank_metasub[order(rank_metasub$number),]
  
  # Label city name and migration situation of each sample
  for (i in 1:length(rank_metasub2$sample)){
    sam <- rank_metasub2$sample[i]
    country <- strsplit(sam,'_migrat')[[1]][1]
    mig <- strsplit(strsplit(sam,'_migrat')[[1]][2],' ',fixed = T)[[1]][1]
    if (!is.na(mig)){
      rank_metasub2[i,'Migration'] <- 'non-native'
      country <- strsplit(sam,'_migrat')[[1]][1]
      if ('_' %in% strsplit(country,'')[[1]]){
        country <- strsplit(country,'_')[[1]]
        for (a in 1:length(country)){
          country_a <- strsplit(country[a],'')[[1]]
          country_a[1] <- toupper(country_a[1])
          country_a <- paste(country_a, collapse='')
          country[a] <- country_a}
        country <- paste(country, collapse=' ')
      }else{
        country <- strsplit(country,'')[[1]]
        country[1] <- toupper(country[1])
        country <- paste(country, collapse='')
      }
      rank_metasub2[i,'country'] <- strsplit(sam,'_migrat')[[1]][1]
      rank_metasub2[i,'country2'] <- country
    }else{
      rank_metasub2[i,'Migration'] <- 'native'
      country <- strsplit(sam,' ')[[1]][1]
      if ('_' %in% strsplit(country,'')[[1]]){
        country <- strsplit(country,'_')[[1]]
        for (a in 1:length(country)){
          country_a <- strsplit(country[a],'')[[1]]
          country_a[1] <- toupper(country_a[1])
          country_a <- paste(country_a, collapse='')
          country[a] <- country_a}
        country <- paste(country, collapse=' ')
      }else{
        country <- strsplit(country,'')[[1]]
        country[1] <- toupper(country[1])
        country <- paste(country, collapse='')
      }
      rank_metasub2[i,'country'] <- strsplit(sam,' ')[[1]][1]
      rank_metasub2[i,'country2'] <- country
    }
  }

  city_mig_p  <- data.frame(row.names = levels(fct_inorder(rank_metasub2$country)))
  for (city in levels(fct_inorder(rank_metasub2$country))){
    correctly_pred <-  round( mean(metasub_data[metasub_data$city == city,"cityPred"]== 
                                     metasub_data[metasub_data$city == city,"city"]) * 100 ,1)
    incorrectly_pred <- round( (100 - correctly_pred) ,1)
    city_mig_p[city,'native'] <- correctly_pred
    city_mig_p[city,'non-native'] <- incorrectly_pred
  }
  
  for (city_ks in levels(factor(rank_metasub2$country2))){
    city_sample <- rank_metasub2[rank_metasub2$country2 == city_ks,]
    native_s <- city_sample[city_sample$sample == levels(factor(city_sample$sample))[1],][,'abundance']
    non_native_s <- city_sample[city_sample$sample == levels(factor(city_sample$sample))[2],][,'abundance']
    
    ks_result <- ks.test(native_s,non_native_s)
    ks_p <- ks_result$p.value
    psign[city_ks,'p.value'] <- ks_p
  }
  
  for (o in 1:dim(city_mig_p)[1]){
    label_o <- paste0('Local: ',city_mig_p[o,'native'],'%\n','Non-local: ',city_mig_p[o,'non-native'],'%')
    for (m in 1:length(rank_metasub2$country)){
      if (rank_metasub2$country[m] == rownames(city_mig_p)[o]){
        p <- psign[psign$city == rank_metasub2$country2[m],'p.value'] # p.sign
        rank_metasub2[m,'text2'] <- label_o
        rank_metasub2[m,'text3'] <- paste0("italic('p-value')~","`: " ,p,"`")
        break
      }
    }
  }
  
  name1 <- paste0("RAC for ",c," cities")
  
  p = ggplot(rank_metasub2,aes(rank,log(abundance,10)))+
    geom_line(aes(linetype=Migration),color = color[num])+ 
    scale_linetype_manual(values=c( "solid","longdash"))+
    facet_wrap(~country2) +
    xlim(0,2500)+
    ylim(-6,-0.3)+
    geom_text(rank_metasub2,mapping=aes(x=rep(2490,length(country2)),y=rep(-1.5,length(country2)),label=text2),nudge_x=0.1,nudge_y=0.1,size =3.2, hjust = 1,colour='black')+
    geom_text_repel(rank_metasub2,mapping=aes(x=rep(2500,length(country2)),y=rep(-2.8,length(country2)),label=text3),size =3.2,colour='black',parse = TRUE,force=0,hjust  =  0)+
    scale_size_manual(values=c(2, 2))+
    ggtitle(name1)+
    labs(x = '', y = '', color = NULL)+
    theme(panel.grid = element_blank(),panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill='transparent'),legend.position="none") 
  name2 <- paste0('rank abundance curve e6 - ',c,' city RAc.png')
  ggsave(name2,p)
}
```

## 5 Antimicrobial resistance genes of local and non-local samples differentiated by mGPS  
The input files "MetaSUB City Metadata - Sheet1.csv" was provided in Github. The file "After_process_megares_amr_class_rpkmg.csv" was provided by Danko, et al., (See in Material and Method part in thesis).  
```r
library(ggplot2)
library(ggsignif)
library(forcats)
library(ggpmisc)

metasub_data <- read.csv('MetaSUB_prediction_result.csv')
AMR_data <- read.csv('After_process_megares_amr_class_rpkmg.csv', check.names = F)

metasub_ARM<- merge(metasub_data,AMR_data,by.x="uuid",by.y="uuid")

correspond <- c()
missing <- c()
for (i in metasub_data$uuid){
  if (i %in% AMR_data$X){
    correspond <- c(correspond,i)
  }else{
    missing <- c(missing,i)
  }
}

AMR_list <- colnames(AMR_data)[2:22]
metasub_ARM$city <- factor(metasub_ARM$city)

# Calculate average AMR relative abundance of local and non-local samples in each city
bar_df_abundance <- data.frame(row.names = c(levels(metasub_ARM$city),paste0(levels(metasub_ARM$city),':migrate')))

for (i in 1:length(levels(metasub_ARM$city))){
  
  # Calculate average relative abundance of each AMR
  for (j in 1: length(AMR_list)){
    
    city <- levels(metasub_ARM$city)[i]
    taxa <- AMR_list[j]  
    city_subset <- subset(metasub_ARM[metasub_ARM$city == city,])
    # Separate local and non-local samples
    origin_same <- subset(city_subset[city_subset$city == city_subset$cityPred,])
    origin_diff <- subset(city_subset[city_subset$city != city_subset$cityPred,])
    
    abun_same <- sum(origin_same[,taxa])/length(origin_same[,taxa])
    abun_diff <- sum(origin_diff[,taxa])/length(origin_diff[,taxa])
    
    bar_df_abundance[city,taxa] <- abun_same
    
    migrate <- paste0(city,':migrate')
    bar_df_abundance[migrate,taxa] <- abun_diff
  }
}

# Calculate average total AMR relative abundance
# Average value of all samples
global_native <- c()
global_non <- c()
# Average value of all cities' value, which is the average value of all samples in that city
global_city_native <- c()
global_city_non <- c()
for (j in 1:length(AMR_list)){
  taxa <- AMR_list[j]  
  # Divide local and non-local samples
  origin_same <- subset(metasub_ARM[metasub_ARM$city == metasub_ARM$cityPred,])
  origin_diff <- subset(metasub_ARM[metasub_ARM$city != metasub_ARM$cityPred,])
  
  abun_same <- mean(origin_same[,taxa])
  abun_diff <- mean(origin_diff[,taxa])
  
  global_native <- c(global_native,abun_same)
  global_non <- c(global_non,abun_diff)
  
  global_city_native <- c(global_city_native, mean(bar_df_abundance[c(1:38,40),taxa]))
  global_city_non <- c(global_city_non,mean(bar_df_abundance[41:80,taxa]))
}

## global analysis

global_avg_df_1 <- cbind.data.frame('AMR'= AMR_list,'Source'=rep('Local',21),'value'=global_native)
global_avg_df_2 <- cbind.data.frame('AMR'= AMR_list,'Source'=rep('Non-local',21),'value'=global_non)
global_avg_df <- rbind.data.frame(global_avg_df_1,global_avg_df_2)

global_avg_df$AMR <- factor(global_avg_df$AMR)
global_avg_df$AMR <- fct_inorder(global_avg_df$AMR)
global_avg_df$AMR <- fct_relevel(global_avg_df$AMR,
                                 rev(levels(global_avg_df$AMR))) 

global_avg_df$Source <- factor(global_avg_df$Source)
global_avg_df$Source <- fct_relevel(global_avg_df$Source,
                                    rev(levels(global_avg_df$Source))) 

ggplot(global_avg_df,aes(x=AMR,y=value,fill=Source))+
  geom_bar(stat="identity",position=position_dodge(0.75))+
  coord_flip()+
  ylab('AMR relative abundance (RPKM)')+
  ggtitle('AMR relative abundance on global level')+
  theme_bw()

## Total AMR distribution of local and non-local on global level

for (i in 1:length(metasub_ARM$uuid)){
  if (metasub_ARM[i,'city'] == metasub_ARM[i,'cityPred']){
    metasub_ARM[i,'Migration'] = 'Local'
  }else{
    metasub_ARM[i,'Migration'] = 'Non-local'
  }
}

ggplot(metasub_ARM, aes(x=Migration,y=`AMR (total)`,fill=Migration)) + 
  geom_boxplot() +
  ylim(0,100)+
  ylab('AMR relative abundance (RPKM)')+
  xlab('Sample origin')+
  guides(fill=FALSE) +
  ggtitle('AMR comparison for native and non-native samples (wilcox.test)')+
  geom_signif(comparisons = list(c("Native", "Non-native")),
              step_increase = 0.3,
              map_signif_level = T,
              test = wilcox.test)+
  theme_bw()

## Relative abundance distribution of each AMR in local and non-local group
# boxplot for each AMR

AMR_df <- subset(metasub_ARM[,c(AMR_list[1:20],'Migration')])
reshape2::melt(AMR_df,id.vars="Migration") -> dfa0

dfa0$Migration <- factor(dfa0$Migration)
levels(dfa0$Migration) <- c('Local','Non-local')

ggplot(dfa0, aes(x=Migration,y=value,fill=Migration)) + 
  geom_boxplot() +
  ylim(0,100)+
  ylab('AMR relative abundance (RPKM)')+
  xlab('Sample source')+
  guides(fill=FALSE) +
  ggtitle('AMR comparison for local and non-local samples (wilcox.test)')+
  facet_wrap(~variable) +
  geom_signif(comparisons = list(c("Local", "Non-local")),
              step_increase = 0.3,
              map_signif_level = T,
              test = wilcox.test)+
  theme_bw()

## The total AMR relative abundance for local and non-local in each city

bar_df_abundance2 <- bar_df_abundance
bar_df_abundance2$city_1 <- row.names(bar_df_abundance)
reshape2::melt(bar_df_abundance2,id.vars="city_1") -> dfa1

# Label city name and migration situation (local and non-local)
for (i in 1:length(dfa1$city_1)){
  city <- dfa1$city_1[i]
  if (':' %in% strsplit(city,'')[[1]]){
    country <- strsplit(city,':')[[1]][1]
    if ('_' %in% strsplit(country,'')[[1]]){
      country <- strsplit(country,'_')[[1]]
      for (a in 1:length(country)){
        country_a <- strsplit(country[a],'')[[1]]
        country_a[1] <- toupper(country_a[1])
        country_a <- paste(country_a, collapse='')
        country[a] <- country_a}
      country <- paste(country, collapse=' ')
    }else{
      country <- strsplit(country,'')[[1]]
      country[1] <- toupper(country[1])
      country <- paste(country, collapse='')
    }
    dfa1[i,'city'] <- country
    dfa1[i,'Sample origin'] <- 'Non-native'
  }else{
    country <- city
    if ('_' %in% strsplit(country,'')[[1]]){
      country <- strsplit(country,'_')[[1]]
      for (a in 1:length(country)){
        country_a <- strsplit(country[a],'')[[1]]
        country_a[1] <- toupper(country_a[1])
        country_a <- paste(country_a, collapse='')
        country[a] <- country_a}
      country <- paste(country, collapse=' ')
    }else{
      country <- strsplit(country,'')[[1]]
      country[1] <- toupper(country[1])
      country <- paste(country, collapse='')
    }
    dfa1[i,'city'] <- country
    dfa1[i,'Sample origin'] <- 'Native'
  }
}

# AMR_Sum analysis for each city
dfa1$variable <- as.character(dfa1$variable)
ARM_sum_df <- subset(dfa1[dfa1$variable == 'AMR (total)',])

city_anno <- read.csv('MetaSUB City Metadata - Sheet1.csv')

# label country name and continent name
for (i in 1:length(ARM_sum_df$city)){
  city <- ARM_sum_df$city_1[i]
  if (':' %in% strsplit(city,'')[[1]]){
    city <- strsplit(city,':')[[1]][1]
  }
  ARM_sum_df[i,'country'] <- city_anno[city_anno$city == city,'country']
  ARM_sum_df[i,'continent'] <- city_anno[city_anno$city == city,'continent']
}

ARM_sum_df <- ARM_sum_df[order(ARM_sum_df$country),]
ARM_sum_df <- ARM_sum_df[order(ARM_sum_df$continent),]

# order the table according to the AMR value
for (i in 1:length(ARM_sum_df$city_1)){
  if (':' %in% strsplit(ARM_sum_df$city_1[i],'')[[1]]){
    ARM_sum_df[i,'order'] <- ARM_sum_df[i,'value']
    city_m <- strsplit(ARM_sum_df$city_1[i],':')[[1]][1]
    ARM_sum_df[ARM_sum_df$city_1 == city_m,'order'] <- ARM_sum_df[i,'value']
  }
}

ARM_sum_df2 <- ARM_sum_df[order(ARM_sum_df$order),]
ARM_sum_df2$city <- factor(ARM_sum_df2$city)
ARM_sum_df2$city <- fct_inorder(ARM_sum_df2$city)
ARM_sum_df2$city <- fct_relevel(ARM_sum_df2$city,rev(levels(ARM_sum_df2$city))) 

ARM_sum_df3 <- ARM_sum_df2
ARM_sum_df3[length(ARM_sum_df3$city_1),'value'] <- 10 # This value is too large to plot

ARM_sum_df2$`Sample origin` <- factor(ARM_sum_df2$`Sample origin`)
ARM_sum_df2$`Sample origin` <- fct_inorder(ARM_sum_df2$`Sample origin`)
ARM_sum_df2$`Sample origin` <- fct_relevel(ARM_sum_df2$`Sample origin`,rev(levels(ARM_sum_df2$`Sample origin`))) 

ggplot(ARM_sum_df3,aes(x=city,y=value,fill=`Sample origin`))+
  geom_bar(stat="identity",position=position_dodge(1))+
  ylab('AMR relative abundance (RPKM)')+
  xlab('City')+
  ggtitle('Sum AMR relative abundance comparison for each city')+
  theme_bw()+
  theme(axis.text.x = element_text( vjust = 0.5, hjust = 1, angle = 90))+
  annotate("text", x=4.7 , y=10,label = "* True value = 13.65", )+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"),
                    breaks = c("Native", "Non-native"))

## Find correlation of tourists and AMR relative abundance of non-local

tour_df <- ARM_sum_df[ARM_sum_df$`Sample origin`=='Non-native',]

for (i in 1:length(tour_df$city_1)){
  city <- strsplit(tour_df$city_1[i],':')[[1]][1]
  tour_df[i,'tourists'] <- city_anno[city_anno$city == city,'number_of_tourists_annually']
}

tour_df$tourists <- as.integer(tour_df$tourists)
tour_df$tourists <- as.numeric(tour_df$tourists)

tour_df <- tour_df[order(tour_df$tourists),]
tour_df$tourists <- factor(tour_df$tourists)
tour_df$tourists <- fct_inorder(tour_df$tourists)
tour_df$tourists <- fct_relevel(tour_df$tourists,
                                rev(levels(tour_df$tourists)))                            

ggplot(data = tour_df,aes(x=tourists,y=value))+
  geom_point(size=2)+
  ylab('AMR relative abundance (RPKM)')+
  xlab('Number of tourists annually')+
  stat_smooth(method=lm,formula=y~x,se=FALSE)+
  geom_smooth (method = lm,linetype=1,se=FALSE,span=1)+
  theme_bw()+
  theme(axis.text.x = element_text( vjust = 0.5, hjust = 1, angle = 90))+
  ggpubr::stat_cor( label.x = 80000000,label.y = 11,color = 'black')+
  stat_poly_eq(
    aes(label = ..eq.label..),
    formula = y~x, parse = TRUE, geom = "text",label.x =80000000,label.y = 12, hjust = 0,color = 'black')
```

## 6 mGPS implementation on monument dataset  

### 6.1 Pre-process of the data and the mGPS input file
The data files (monumentome_read_abundance.csv and Metadata_Monumentome.xlsx) were provided by Monumentome Project (unpublished).  
In Metadata_Monumentome.xlsx, manually extracted following columns: sample_id (the first one), pangea_barcode, project, sample_id2 (the seconde sample_id column and rename it), sample_type, city, country, monument, latitude, longitude, continent. Then changed all the comma to underscore . Resaved the file into CSV format as 'metadata_monumentome.csv'

```r
library(ggplot2)
library(forcats)

monument_read_abun <- read.csv('monumentome_read_abundance.csv',check.names = F)
monument_metadata <- read.csv('metadata_monumentome.csv',check.names = F)

# read filter: Remove read data less than 100
read_posi <- c()
for (i in 1:length(monument_read_abun$name)){
  sample <- monument_read_abun$name[i]
  sample_id <- strsplit(sample,'.',fixed = T)[[1]][1]
  bracken <- strsplit(sample,'.',fixed = T)[[1]][2]
  if (bracken == 'bracken_num'){
    read_posi <- c(read_posi,i)
    rownames(monument_read_abun)[i] <- sample_id
  }
}

monument_read <- monument_read_abun[read_posi,]
monument_read <- monument_read[,-1]
monument_read[monument_read < 100] <- 0

data_normalise <- function(df) {
  return(df/rowSums(df))
}

monument_abun <- monument_read/rowSums(monument_read) # relative abundance data

monument_abun[is.na(monument_abun)] <- 0
monument_abun$sample_id <- rownames(monument_abun)

# metadata process: remove control samples
control_sample_p <- which(monument_metadata$sample_type %in% c("control","negative_control", "positive_control"))
monument_metadata_non_control <- monument_metadata[-c(control_sample_p),]

## create monument data used in mGPS
monument_data <- merge(monument_metadata_non_control,monument_abun,by.x="sample_id",by.y="sample_id")

for (i in 1:length(monument_data$city)){
  if (is.na(monument_data$latitude[i])){
    city <- monument_data$city[i]
    if (city == 'Tel Megiddo'){
      monument_data[i,'latitude'] <- 32.5806
      monument_data[i,'longitude'] <- 35.1795
    }else if (city == 'Ulaanbaatar'){
      monument_data[i,'latitude'] <- 47.8923
      monument_data[i,'longitude'] <- 106.9090
    }else if (city == 'Seoul'){
      monument_data[i,'latitude'] <- 37.5480
      monument_data[i,'longitude'] <- 126.9887
    }
  }
}
write.csv(monument_data,'monument_data_filter.csv',row.names = F)

## relative abundance data table with control samples
monument_data <- merge(monument_metadata,monument_abun,by.x="sample_id",by.y="sample_id")
write.csv(monument_data,'monument_data_filter_with_control.csv',row.names = F)

## read data table with control samples
monument_read$sample_id <- rownames(monument_read)  
monument_read_data <- merge(monument_metadata,monument_read,by.x="sample_id",by.y="sample_id")
write.csv(monument_read_data,'monument_read_data_with_control.csv',row.names = F)
```

### 6.2 Monument sample prediction by mGPS interface
Used mGPS interface to construct microbiome source prediction model for monument dataset and predict the source of samples in original monument dataset.  
OPen mGPS_interface.r, click *RunApp*.  
Select `Prediction program` as `Build a new prediction model using mGPS`  
Select `Merged metadata and abundance data`  
`Input file`: monument_data_filter.csv  
`Enter the main locality level`: city  
`Enter the locality hierarchy`: country,city,latitude,longitude  
`Column range of abundance data`: 12:9173  
`Subsets in feature elimination`: 10,20,30,50,80,100,120,150,180,200,250,300,350,400,450,500,800,1000  
Click `Start` and `Result Plot` tab  
When application done, click `Pull to land`  
Click Output tab to download the result files: `Download prediction data`: Prediction_results.csv; `Download optimal features in prediction model`: Optimal_features.csv; `Download feature subsets accuracy in feature elimination`: Features_subsets_accuracy.csv   
```r
# Re-plot the predicted sites on world map and added the sampling site
library(rworldmap)
library(mapplots)
library(maps)
library(ggrepel)

monument_data <- read.csv('Prediction_results.csv')

monument_data$continent <- factor(monument_data$continent)
monument_data$country <- factor(monument_data$country)
monument_data$city <- factor(monument_data$city)

monument_data$country <- factor(monument_data$country)
png('monument.png', width = 13,height = 8, units = 'in', res = 600)

par(mai=c(2,1,0.5,0.5), mar=par()$mar+c(3,0,0,0))
map <- rworldmap::getMap(resolution = "coarse")
palette <-c("brown","gold2","orangered2","mediumspringgreen","darkorchid4","deeppink2","dodgerblue3")

plot(map, xlim = c(-168,168),ylim = c(-84,78), col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.5)
#find coord preds by region
for ( i in 1:length(levels(monument_data$country))){
  this_continent <- levels(monument_data$country)[i]
  find_lats <- monument_data[monument_data[,"country"] == this_continent,][,"latPred"]
  find_longs <- monument_data[monument_data[,"country"] == this_continent,][,"longPred"]
  #plot predicted co-ordinates
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.2)
}

legend(-186,-26, c("Brazil","Greece","Israel","Japan","Mongolia","South Korea","United States"), pch = 17, col = palette, cex = 1.3, bg ="lightskyblue1",ncol = 1)

#Plot city sampling locations
map.axes(cex.axis = 1.3)
par(fig = c(0.358,0.817,0,0.56), new = T) #par(fig = c(0.45,0.9,0.45,0.5), new = T) # 
plot(map,xlim = c(-160,168), ylim = c(-50,76), col = "grey", border = "darkgrey", bg ="lightskyblue1")
for ( i in 1:length(levels(monument_data$country))){
  this_continent <- levels(monument_data$country)[i]
  find_lats <- monument_data[monument_data$country == this_continent,]$latitude
  find_longs <- monument_data[monument_data$country == this_continent,]$longitude
  
  points(find_longs, find_lats, col = palette[i], pch = 17, cex = 1)
}
box( col = 'black')
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

## Test the correlation between sample size and model accuracy

city_sample_sizes <- data.frame()
for (i in 1:length(levels(monument_data$city))){ 
  this_city <- levels(monument_data$city)[i]
  city_sample_sizes[i,"city_size"] <-  summary(monument_data$city)[this_city]
  city_sample_sizes[i,"within_500"] <-  mean(monument_data[monument_data$city == this_city,]$Distance_from_origin < 500)
}

#plot
lowess_line <- lowess(city_sample_sizes$city_size,100*city_sample_sizes[,"within_500"])

ggplot(data = city_sample_sizes,aes(x=city_size,y=100*within_500))+
  geom_point(size=2)+
  ylim(-10,70)+
  xlim(-5,55)+
  ylab('Predictions within 500km (%)')+
  xlab('Number of samples from city')+
  ggtitle('Hierarchical level: country-city-latitude-longitude')+
  geom_line(aes(x=lowess_line$x,y=lowess_line$y),colour = "red")+
  geom_text_repel(city_sample_sizes,mapping = aes(x=city_size,y=100*within_500,label=levels(MetasubDataPreds$city)),size=3.5)+
  theme_bw()+
  ggpubr::stat_cor(label.x = -5,label.y = 70,color = 'black')

## effect of removing redundant taxa
git_subset <- read.csv("Features_subsets_accuracy.csv")
png(file="subsets accuracy.png",width=600, height=350)
plot(git_subset[,"n_vars"], git_subset[,"accuracy"], type = "l", xlab = "", ylab = "", cex = 1.6, col = "dodgerblue")
title(ylab="Accuracy of random forest classifier",xlab = "Number of taxa", mgp=c(2.5,1.5,1),cex.lab=1.2)
points(git_subset[4,"n_vars"],git_subset[4,"accuracy"], pch = 20, cex = 1.5,col = "red")
dev.off()

## Monumentome GITs (top 25)
v <- read.csv("Optimal_features.csv")
top_species <- v[1:25,"taxa"]
#plot
png('optimal_features.png', width = 12,height = 7, units = 'in', res = 600)
par(font = 3)
dotchart(rev(v[1:25,"Overall"])*100,labels= rev(top_species),cex=1.2,pt.cex = 1.3,
         xlab="Mean decrease in accuracy", mgp = c(2.2,0,0))
dev.off()
```

### 6.3 mGPS model for monument dataset based on MetaSUB GITs
Used monument dataset to train the mGPS model, but the feature used in training was the taxa in monument dataset that also can be found in MetaSUB GITs.
```r
# Load the following functions in mGPS interface script:
# "mGPS", "model_accuracy_f", and "pull_land"

monument_data <- read.csv('monument_data_filter.csv')
metasub_git <- read.csv('global_GITS.csv')

optVars <- c() 
for (i in metasub_git$taxa){
  if (i %in% colnames(monument_data)){
    optVars <- c(optVars,i) 
  }
}

prediction_output_before <- model_accuracy_f(monument_data, optVars,'city',c('country','city','latitude','longitude'))

df <- prediction_output_before[[1]]
prediction_result <- pull_land(df,c('country','city','latitude','longitude'))

## plot the result figures
# predicted site on world map:

library(rworldmap)
library(mapplots)
library(maps)

monument_data <- prediction_result

monument_data$continent <- factor(monument_data$continent)
monument_data$country <- factor(monument_data$country)
monument_data$city <- factor(monument_data$city)

monument_data$country <- factor(monument_data$country)
png('monument_metasub_git.png', width = 13,height = 8, units = 'in', res = 600)

par(mai=c(2,1,0.5,0.5), mar=par()$mar+c(3,0,0,0))

map <- rworldmap::getMap(resolution = "coarse")

palette <-c("brown","gold2","orangered2","mediumspringgreen","darkorchid4","deeppink2","dodgerblue3")

plot(map, xlim = c(-168,168),ylim = c(-84,78), col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.5)
#find coord preds by region
for ( i in 1:length(levels(monument_data$country))){
  this_continent <- levels(monument_data$country)[i]
  find_lats <- monument_data[monument_data[,"country"] == this_continent,][,"latPred"]
  find_longs <- monument_data[monument_data[,"country"] == this_continent,][,"longPred"]
  #plot predicted co-ordinates
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.2)
}

legend(-186,-26, c("Brazil","Greece","Israel","Japan","Mongolia","South Korea","United States"), pch = 17, col = palette, cex = 1.3, bg ="lightskyblue1",ncol = 1)

#Plot city sampling locations
map.axes(cex.axis = 1.3)
par(fig = c(0.358,0.817,0,0.56), new = T) #par(fig = c(0.45,0.9,0.45,0.5), new = T) # 
plot(map,xlim = c(-160,168), ylim = c(-50,76), col = "grey", border = "darkgrey", bg ="lightskyblue1")
for ( i in 1:length(levels(monument_data$country))){
  this_continent <- levels(monument_data$country)[i]
  find_lats <- monument_data[monument_data$country == this_continent,]$latitude
  find_longs <- monument_data[monument_data$country == this_continent,]$longitude
  
  points(find_longs, find_lats, col = palette[i], pch = 17, cex = 1)
}
box( col = 'black')
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

# accuracy bar plot:

bar_df1 <- data.frame(row.names = c("Overall",levels(monument_data$city)))
for (i in 1: length(levels(monument_data$city))){
  this_city <- levels(monument_data$city)[i]
  prop <- mean(monument_data[monument_data$city == this_city,][,"Distance_from_origin"] < 100)
  bar_df1[i+1,"0 - 100km"] <- prop
  
  overall_prop <- mean(monument_data[,"Distance_from_origin"] < 100)
  bar_df1[ 1,"0 - 100km"] <- overall_prop
}
for (i in 1: length(levels(monument_data$city))){
  this_city <- levels(monument_data$city)[i]
  prop <- mean(monument_data[monument_data$city == this_city,][,"Distance_from_origin"] > 100 & monument_data[monument_data$city == this_city,][,"Distance_from_origin"] < 500)
  bar_df1[i+1,"100 - 500km"] <- prop
  
  overall_prop <-mean(monument_data[,"Distance_from_origin"] > 100 & monument_data[,"Distance_from_origin"] < 500)
  bar_df1[ 1,"100 - 500km"] <- overall_prop
}
for (i in 1: length(levels(monument_data$city))){
  this_city <- levels(monument_data$city)[i]
  prop <- mean(monument_data[monument_data$city == this_city,][,"Distance_from_origin"] > 500 & monument_data[monument_data$city == this_city,][,"Distance_from_origin"] < 1000)
  bar_df1[i+1,"500 - 1000km"] <- prop
  
  overall_prop <- mean(monument_data[,"Distance_from_origin"] > 500 & monument_data[,"Distance_from_origin"] < 1000)
  bar_df1[ 1,"500 - 1000km"] <- overall_prop
}
for (i in 1: length(levels(monument_data$city))){
  this_city <- levels(monument_data$city)[i]
  prop <- mean(monument_data[monument_data$city == this_city,][,"Distance_from_origin"] > 1000 & monument_data[monument_data$city == this_city,][,"Distance_from_origin"] < 2000)
  bar_df1[i+1,"1000 - 2000km"] <- prop
  
  overall_prop <- mean(monument_data[,"Distance_from_origin"] > 1000 & monument_data[,"Distance_from_origin"] < 2000)
  bar_df1[ 1,"1000 - 2000km"] <- overall_prop
}
for (i in 1: length(levels(monument_data$city))){
  this_city <- levels(monument_data$city)[i]
  prop <- mean(monument_data[monument_data$city == this_city,][,"Distance_from_origin"] > 2000 & monument_data[monument_data$city == this_city,][,"Distance_from_origin"] < 3000)
  bar_df1[i+1,"2000 - 3000km"] <- prop
  
  overall_prop <- mean(monument_data[,"Distance_from_origin"] > 2000 & monument_data[,"Distance_from_origin"] < 3000)
  bar_df1[1,"2000 - 3000km"] <- overall_prop
}
for (i in 1: length(levels(monument_data$city))){
  this_city <- levels(monument_data$city)[i]
  prop <- mean(monument_data[monument_data$city == this_city,][,"Distance_from_origin"] > 3000 )
  bar_df1[i+1,"> 3000km"] <- prop
  
  overall_prop <- mean(monument_data[,"Distance_from_origin"] > 3000)
  bar_df1[ 1,"> 3000km"] <- overall_prop
}
size <- c()
for (i in 1: length(levels(monument_data$city))){
  
  this_city <- levels(monument_data$city)[i]
  size[i] <- length(which(monument_data$city == this_city))
}

png('accuracy_metasub_git_monument.png', width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(4,1,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df1*100), space = 0,col=c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), 
              names.arg=c("Overall",paste0(city_names,"  (",size,")"), axes = FALSE) , 
              las =2, cex.names=1.2, ylab = "", axisnames = F, axes = F)
axis(side =2, pos = 0)
mtext(text = c("Overall",paste0(levels(MetasubDataPreds$city)," (",size,")")), side = 1, at = bp, line = 0, padj = 1, cex = 1.2)
title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1.2)
legend("topright",inset = c(-0.16,0.4), rev(c(colnames(bar_df1))), 
       fill = rev(c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) ,
       bty = 1, cex = 1.2)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()
```
### 6.4 Source prediction of monument sample based on MetaSUB mGPS model
Used Module three in mGPS interface. The MetaSUB mGPS model had been downloaded when applied the MetaSUB dataset in mGPS.  
OPen mGPS_interface.r, click *RunApp*.  
Select `Prediction program` as `Use an existing model to predict new samples`  
`Upload new sample(s) abundance file`: monument_data_filter.csv  
`Upload the prediction model`: MetaSub_model.Rda  
Click `Start` and `Result Plot` tab  

### 6.5 Exploration of the monument GITs
#### 6.5.1 Taxonomical profile annotation of monument taxa
Used python to annotate the taxonomical profile and pathogenicity of taxa in monument dataset using Microbiome Dictionary v2 and GTDB (See in Material and Method in thesis). Outputted the file "monument_taxonomic_level.txt"  
```python
# Same as the annotation of MetaSUB taxa
from collections import defaultdict

taxa = open("monument_data_filter.csv",'r',encoding='UTF-8')
directory = open("microbe-directory.csv",'r',encoding='UTF-8')
fileout = open("monument_taxonomic_level.txt",'w',encoding='UTF-8')

bac120 = open("bac120_taxonomy.tsv",'r',encoding='UTF-8')
bac120_dic =  defaultdict(list)
for line in bac120:
    level = line.strip().split("\t")[1]
    taxa_l = level.split(";")
    species = taxa_l[6].split("__")[1] # species level
    if species in bac120_dic.keys():
        continue
    for i in taxa_l:
        t = i.split("__")[1]
        bac120_dic[species].append(t) # d p c o f g s
            
bac120_dic2 = {}
for line in bac120:
    level = line.strip().split("\t")[1]
    species = level.split(";")[6] # species level
    species = species.split("__")[1]
    if species in bac120_dic2:
        continue
    bac120_dic2[species] = level
   
ar122 = open("ar122_taxonomy.tsv",'r',encoding='UTF-8')
ar122_dic = defaultdict(list)
for line in ar122:
    level = line.strip().split("\t")[1]
    taxa_l = level.split(";")
    species = taxa_l[6].split("__")[1] # species level
    if species in ar122_dic:
        continue
    for i in taxa_l:
        t = i.split("__")[1]
        ar122_dic[species].append(t) # d p c o f g s

ssu = open("ssu_silva_taxonomy",'r',encoding='UTF-8')
ssu_dic = defaultdict(list)
for line in ssu:
    taxa_l = line.strip().split(";")
    if taxa_l[-1] in ssu_dic:
        continue
    for i in taxa_l:
        ssu_dic[taxa_l[-1]].append(i)
  
ncbi = open("ncbi_taxonomy",'r',encoding='UTF-8')
ncbi_dic = defaultdict(list)
for line in ncbi:
    taxa_l = line.strip().split(";")
    s_ncbi = taxa_l[-1]
    if "__" in s_ncbi:
        s_ncbi = s_ncbi.split("__")[1]
        if s_ncbi in ncbi_dic:
            continue
        for i in taxa_l:
            t = i.split("__")[1]
            ncbi_dic[s_ncbi].append(t)
  
ar_ssu = open("ar122_ssu_silva_taxonomy",'r',encoding='UTF-8')
ar_ssu_dic = defaultdict(list)
for line in ar_ssu:
    taxa_l = line.strip().split(";")
    species = taxa_l[-1]
    if species in ar_ssu_dic:
        continue
    for i in taxa_l:
        ar_ssu_dic[species].append(i)
    
ar_ncbi = open("ar122_ncbi_taxonomy",'r',encoding='UTF-8')
ar_ncbi_dic = defaultdict(list)
for line in ar_ncbi:
    taxa_l = line.strip().split(";")
    s_ncbi = taxa_l[-1]
    if "__" in s_ncbi:
        s_ncbi = s_ncbi.split("__")[1]   
        if s_ncbi in ar_ncbi_dic:
            continue
        for i in taxa_l:
            t = i.split("__")[1]
            ar_ncbi_dic[s_ncbi].append(t)

directory_dic =  defaultdict(list) 
for line in directory:
    line = line.strip().split(",")
    if line[7] in directory_dic:
            continue
    for n in range(1,8):
        directory_dic[line[7]].append(line[n]) # kingdom: directory_dic[taxa][0]
    directory_dic[line[7]].append(line[22]) # animal_pathogen
    directory_dic[line[7]].append(line[28]) # plant_pathogen
    
taxa_list = taxa.readline().strip().split(',')[1:]
anno_file = [directory_dic,bac120_dic,ssu_dic,ncbi_dic,ar122_dic,ar_ssu_dic,ar_ncbi_dic]
flag = {1:'directory_dic',2:'bac120_dic',3:'ssu_dic',4:'ncbi_dic',5:'ar122_dic',6:'ar_ssu_dic',7:'ar_ncbi_dic'}
for t in taxa_list:
    t = t.strip('"')
    flag_count = 0
    for i in anno_file:
        flag_count += 1
        if t in i:
            anno = ''
            for feature in i[t]:
                if anno:
                    anno = "{}\t{}".format(anno,feature)
                else:
                    anno = feature
print('{}\t{}\t{}'.format(flag[flag_count],t,anno),file=fileout)
            break
taxa.close()
directory.close()
fileout.close()
bac120.close()
ar122.close()
ssu.close()
ncbi.close()
ar_ssu.close()
ar_ncbi.close()
```
Opened file monument_taxonomic_level.txt in excel to manually fix the heterotypic synonym of taxa in the taxonomical profile.  
Added a new column named "pathogen_type". According to the pathogen annotation in file, labeled animal pathogen as "1", plant pathogen as "2", dual-kongdom pathogen as "3".
Resaved the revised file into CSV format "monument_taxonomic_level.csv"

#### 6.5.2 Monument GITs exploration
```r
# GITs in control samples 
monument_metadata <- read.csv('metadata_monumentome.csv',check.names = F)
git_file <- read.csv('Optimal_features.csv',header = T)

control_sample_p <- which(monument_metadata$sample_type %in% c("control","negative_control", "positive_control"))
monument_metadata_control <- monument_metadata[c(control_sample_p),]

monument_control_data <- merge(monument_metadata_control,monument_abun,by.x="sample_id",by.y="sample_id")

colnames(monument_control_data) <- make.names(colnames(monument_control_data))
rownames(monument_control_data) <- monument_control_data$sample_id

monument_control_gits_table <- monument_control_data[,git_file$taxa]
monument_control_gits_table[monument_control_gits_table > 0] <- 1 
rowSums(monument_control_gits_table)

control_have_git_number <- data.frame(rowSums(monument_control_gits_table))
# write.csv(control_have_git_number,'control_have_git_number.csv')

git_n_in_18_control <- colSums(monument_control_gits_table)
git_n_in_18_control <- data.frame(git_n_in_18_control)
# write.csv(git_n_in_18_control2,'git_n_in_18_control.csv')
# Use excel to plot the figure: The number of control samples that GITs appeared

## GITs relative abundance per city
monument_data <- read.csv('monument_data_filter.csv',check.names = F,header = T)

rownames(monument_data) <- monument_data$sample_id
po <- which(make.names(colnames(monument_data)) %in% git_file$taxa)
monument_gits <- monument_data[,po]

monument_gits$gits_RA <- rowSums(monument_gits)
monument_gits$non_gits_RA <- 1- monument_gits$gits_RA

for (i in 1:dim(monument_gits)[1]){
  sample <- rownames(monument_gits)[i]
  monument_gits[i,'city'] <- monument_data[monument_data$sample_id==sample,'city']
  monument_gits[i,'continent'] <- monument_data[monument_data$sample_id==sample,'continent']
}

monument_gits <-monument_gits[order(monument_gits$continent),]

monument_gits$city <- factor(monument_gits$city)
monument_gits$city <- fct_inorder(monument_gits$city)  

city_git <- data.frame(row.names = levels(monument_gits$city))
for (i in 1:length(levels(monument_gits$city))){
  city <- levels(monument_gits$city)[i]
  city_git[i,'gits_RA'] <- mean(monument_gits[monument_gits$city==city,'gits_RA'])
  city_git[i,'non_gits_RA'] <- mean(monument_gits[monument_gits$city==city,'non_gits_RA'])
  city_git[i,'git_n'] <- sum(colSums(monument_gits[monument_gits$city==city,1:50]) > 0)
}

global_gits_RA <- mean(city_git[,'gits_RA'])
global_non_gits_RA <- mean(city_git[,'non_gits_RA'])

bar_df <- cbind.data.frame(c(global_gits_RA,city_git$gits_RA),c(global_non_gits_RA,city_git$non_gits_RA))
colnames(bar_df) <- c('GITs','Non GITs')
rownames(bar_df) <- c('Global',rownames(city_git))

# png("monument_git_abundance.png", width = 13,height = 8, units = 'in', res = 600)
par(xpd = T, mar = par()$mar + c(4,0,0,5), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df*100), col=c("royalblue3","slategray1"),
              names.arg=paste0(rownames(bar_df),' (',c(50,city_git$git_n),')' ),
              args.legend = list(x = "topright", inset=c(-0.5,0)), las =2,
              cex.names=1.2,ylab = "", axisnames = F,axes = F, space =0)

axis(side =2, pos = 0)
mtext(text = paste0(rownames(bar_df),' (',c(50,city_git$git_n),')' ), side = 1, at = bp, line = 0, padj = 1, cex = 1.2)
title(ylab="Taxa abundance %", mgp=c(0,0,0),cex.lab=1.3)
legend("topright",inset = c(-0.1,0.4), c('Non-Gits','Gits'), fill = rev(c("royalblue3","slategray1")) , bty = 1, cex = 1.2)

par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()


## city unique GITs ----------

city_git_avg <- data.frame(row.names = colnames(monument_gits)[1:50])
for (i in 1:length(levels(monument_gits$city))){
  city <- levels(monument_gits$city)[i] 
  city_1 <- data.frame(colSums(monument_gits[monument_gits$city==city,1:50]))
  city_git_avg <- cbind.data.frame(city_git_avg,city_1)
  colnames(city_git_avg)[i] <- city
}

city_git_avg[city_git_avg > 0] <- 1

git_share_in_city_n <- data.frame(rowSums(city_git_avg))
# write.csv(git_share_in_city_n,'git_share_in_city_n.csv')
# Use excel to plot the figure: The number of cities that GITs appeared in

## pathogen annotation for GITs 
git_anno_file <- read.csv('monument_taxonomic_level.csv',header = T)

git_name_o <- colnames(monument_gits)[1:50]
git_name_m <- make.names(git_name_o)

for (i in 1:length(git_file$taxa)){
  git_in_file <- git_file$taxa[i]
  git_o_file <- git_name_o[which(git_name_m %in% git_in_file)]
  git_file[i,'taxa_name'] <- git_o_file
  
  if (git_o_file %in% git_anno_file$taxa){
    git_file[i,'Kingdom'] <- git_anno_file[git_anno_file$taxa==git_o_file,'Kingdom']
    git_file[i,'Phylum'] <- git_anno_file[git_anno_file$taxa==git_o_file,'Phylum']
    git_file[i,'Pathogen_type'] <- git_anno_file[git_anno_file$taxa==git_o_file,'pathogen_type']
  }
}

write.csv(git_file,'git_annotation.csv',row.names = F)
```
## 7 Microbiome diversity in Tel Megiddo

### 7.1 Krakenuniq data of Tel Megiddo samples
The Krakenuniq output data of Tel Megiddo samples were provided.  
Filtered the species in each sample with 100 *taxReads* and 1000 *kmers*.  

```bash
ls *.output* | while read line; do awk '{if($3>100 && $4>1000) print $0}' ${line} > filter/${line}.filter; done

# merge files
cd filter
ls | while read line; do cat ${line} | grep -v '#' | grep -v '%' >> ../merge_krakenuniq_filter; done
cd ..
cat merge_krakenuniq_filter | cut -f 9 | sed 's/^ *//' > krakeuniq_species

cat krakeuniq_species | sort | uniq > krakeuniq_species_uniq 
cat krakeuniq_species_uniq | wc -l # 4533 species for all samples; 3917 species for samples in monument dataset
```
### 7.2 Tel Megiddo taxa overview
The file 'TelMegiddo - sampling_final.xlsx' which records the meta data of sample was provided by Eran Elhaik.  
``` r
library(RColorBrewer)
library(forcats)
library(ggplot2)
library(ggpubr)
library(Hmisc) 

monument_file <- read.csv('monument_data_filter_with_control.csv',header = T,check.names = F)

monument_file_filter <- monument_file[-c(which(rowSums(monument_file[,12:7193]) == 0)),-c(which(colSums(monument_file[,12:7193]) == 0) + 11)]

git_file <- read.csv('git_annotation.csv',header = T, check.names = F)

telm_file <- monument_file[monument_file$city == "Tel Megiddo",]

# remove taxa with zero abundance
taxa_all_list <- colnames(telm_file)[12:length(colnames(telm_file))]
taxa_abun_sum <- colSums(telm_file[12:length(colnames(telm_file))])
posi <- which(taxa_abun_sum!=0)
telm_file <- telm_file[,c(colnames(telm_file)[1:11],taxa_all_list[posi])]

telm_git_list <- c()
for (i in git_file$taxa_name){
  if (i %in% taxa_all_list[posi]){
    telm_git_list <- c(telm_git_list,i)
  }else{
    telm_miss_git <- i
  }
}

## taxa count per city

# Create a data form to store the number of times each taxon appears in each city
monument_file$city <- factor(monument_file$city)

taxa_count_per_city <- data.frame(row.names = c(levels(monument_file$city)))

for (i in 1:length(levels(monument_file$city))){
  for (j in 1: length(taxa_all_list)){
    city <- levels(monument_file$city)[i]
    taxa <- taxa_all_list[j]
    count <- sum(monument_file[monument_file$city == city,][,taxa] != 0)
    taxa_count_per_city[city,taxa] <- count
  }
}

telm_uniq_taxa <- c()
for (i in colnames(taxa_count_per_city)){
  if (sum(taxa_count_per_city[,i])>0){
    if (sum(taxa_count_per_city[,i])==taxa_count_per_city["Tel Megiddo",i]){
      telm_uniq_taxa <- c(telm_uniq_taxa,i)
    }
  }
}

# Tel Megiddo uniq taxa
telm_uniq_taxa_file <- data.frame(telm_uniq_taxa)

## Calculate average relative abundance in city and global 
monument_abun <- data.frame(row.names = c(levels(monument_file$city)))

for (i in 1:length(levels(monument_file$city))){
  for (j in 1: length(taxa_all_list)){
    city <- levels(monument_file$city)[i]
    taxa <- taxa_all_list[j]  
    abun <- sum(monument_file[monument_file$city == city,][,taxa])/length(monument_file[monument_file$city == city,][,taxa])
    monument_abun[city,taxa] <- abun
  }
}

# Calculate global average relative abundance
global_proportion <- c()

for (i in 1:(length(taxa_all_list))){
  abun <- sum(monument_abun[,i])/length(monument_abun[,i])
  global_proportion <- c(global_proportion,abun)
}

monument_abun <- rbind.data.frame(global_proportion,monument_abun)
rownames(monument_abun)[1] <- 'Global'

# Global average relative abundance of annotated taxa
taxa_anno_file <- read.csv('monument_taxonomic_level.csv',check.names = F,header = T)

anno_taxa_abun <- sum(monument_abun['Global',taxa_anno_file$taxa])
sum(monument_abun['Global',]) # 1

## Annotate the site, surface material and access information to each sample
library(readxl)
telm_anno_file_2 <- read_excel('TelMegiddo - sampling_final.xlsx',sheet = 1, na = 'NA')

for (i in 1:length(telm_file$pangea_barcode)){
  s_id <- telm_file$pangea_barcode[i]
  s_id <- strsplit(s_id,'-')[[1]][2]
  s_id <- substr(s_id,2,length(strsplit(s_id,'')[[1]]))
  telm_file[i,'Site'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Site_merged']
  telm_file[i,'Material'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Material']
  telm_file[i,'Access'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Access']
}

telm_file2 <- telm_file[,c(1:11,4271,4272,4273,12:4270)]
rownames(telm_file2) <- telm_file2$sample_id

## Filter samples with zero taxa RSA and filter taxa with zero RSA

# Inpur read data of monument samples
monument_read_file <- read.csv('monument_read_data_with_control.csv',check.names = F,header = T)

monument_read_file_filter <- monument_read_file[-c(which(rowSums(monument_read_file[,12:7193]) == 0)),-c(which(colSums(monument_read_file[,12:7193]) == 0)+11)]

# Extract data of samples from Tel Megiddo 
telm_read_file <- monument_read_file_filter[monument_read_file_filter$city == "Tel Megiddo",]

# remove taxa with zero abundance
taxa_all_list2 <- colnames(telm_read_file)[12:length(colnames(telm_read_file))]
taxa_read_sum <- colSums(telm_read_file[12:length(colnames(telm_read_file))])

posi2 <- which(taxa_read_sum!=0)

telm_read_file <- telm_read_file[,c(colnames(telm_read_file)[1:11],taxa_all_list2[posi2])]

for (i in 1:length(telm_read_file$pangea_barcode)){
  s_id <- telm_read_file$pangea_barcode[i]
  s_id <- strsplit(s_id,'-')[[1]][2]
  s_id <- substr(s_id,2,length(strsplit(s_id,'')[[1]]))
  telm_read_file[i,'Site'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Site_merged']
  telm_read_file[i,'Material'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Material']
  telm_read_file[i,'Access'] <- telm_anno_file_2[telm_anno_file_2$number == s_id,'Access']
}

telm_read_file2 <- telm_read_file[,c(1:11,4271,4272,4273,12:4270)] # after read filter
rownames(telm_read_file2) <- telm_read_file2$sample_id

telm_read_data <- telm_read_file2[,15:dim(telm_read_file2)[2]] 
rownames(telm_read_data) <- telm_read_file2$sample_id

telm_taxa_list <- colnames(telm_file2)[15:4273]

## Microbiome diversity analysis

### alpha diversity
library(vegan)

alpha_diversity <- function(otu){
  Read_number <- rowSums(otu)
  Observed_species <- estimateR(otu)[1, ] # Richness index
  Chao1  <- estimateR(otu)[2, ] # Chao 1 index
  ACE  <- estimateR(otu)[4, ] # ACE index
  Shannon <- diversity(otu,'shannon') # shannon index
  Simpson <- diversity(otu,'simpson') # Gini-Simpson index
  Goods_coverage <- 1 - rowSums(otu == 1) / rowSums(otu)
  output <- data.frame(Read_number, Observed_species, Chao1, ACE, Shannon, Simpson, Goods_coverage)
  return(output)
}

telm_alpha_diversity <- alpha_diversity(telm_read_data)

# Remove samples with less than 1000 read number and less than 10 observed species
telm_alpha_diversity <- telm_alpha_diversity[-c(which(telm_alpha_diversity$Read_number<1000)),]
telm_alpha_diversity <- telm_alpha_diversity[-c(which(telm_alpha_diversity$Observed_species<10)),]

# Annotate environmental variables of each sample
for (i in 1:length(rownames(telm_alpha_diversity))){
  id <- rownames(telm_alpha_diversity)[i]
  telm_alpha_diversity[i,'Site'] <- telm_read_file2[telm_read_file2$sample_id == id, 'Site']
  telm_alpha_diversity[i,'Material'] <- telm_read_file2[telm_read_file2$sample_id == id, 'Material']
  telm_alpha_diversity[i,'Access'] <- telm_read_file2[telm_read_file2$sample_id == id, 'Access']
}

sample_id_list <- rownames(telm_alpha_diversity)

# Filter RSA data and read data
telm_file2 <- telm_file2[sample_id_list,] 
telm_read_file2 <- telm_read_file2[sample_id_list,]

telm_read_data <- telm_read_file2[,15:dim(telm_read_data)[2]] # just contain read data
telm_abun_data <- telm_file2[,15:dim(telm_file2)[2]] # just contain abundance data

## taxa richness distribtion in samples
taxa_Richness <- as.data.frame(specnumber(telm_read_data,MARGIN = 2))

p =ggplot(taxa_Richness, aes(x = `specnumber(telm_read_data, MARGIN = 2)`)) +
  geom_density(alpha = 0.3)+
  xlab("Number of samples")+
  theme_bw()
# ggsave('taxa richness in samples.png',p)  

## Kingdom relative abundance in total Tel Megiddo dataset

anno_num <- which(colnames(telm_file2) %in% taxa_anno_file$taxa)
sum(colSums(telm_file2[,which(colnames(telm_file2) %in% taxa_anno_file$taxa)])/dim(telm_file2)[1])

# Taxa average RSA in samples
taxa_avg_abun <- colSums(telm_file2[,15:dim(telm_file2)[2]])/dim(telm_file2)[1]
sum(taxa_avg_abun) 
taxa_avg_abun <- as.data.frame(taxa_avg_abun)

# Extract annotation table for Tel Megiddo samples
taxa_anno_file_telm <- taxa_anno_file[which(taxa_anno_file$taxa %in% telm_taxa_list),]

taxa_anno_file_telm <- droplevels(taxa_anno_file_telm)
taxa_anno_file_telm$Phylum <- factor(taxa_anno_file_telm$Phylum)
taxa_anno_file_telm$Kingdom <- factor(taxa_anno_file_telm$Kingdom)

# Create a table to preserve the average RSA of each kingdom
kingdom_df <- data.frame(row.names = c(levels(taxa_anno_file_telm$Kingdom),'No-annotation'))

annotation_p <- 0
for (i in 1:length(levels(taxa_anno_file_telm$Kingdom))){
  kingdom <- levels(taxa_anno_file_telm$Kingdom)[i]
  taxa_list <- taxa_anno_file_telm[taxa_anno_file_telm$Kingdom == kingdom,][,'taxa']
  percentage <- sum(as.numeric(taxa_avg_abun[taxa_list,'taxa_avg_abun']))
  kingdom_df[i,'abundance'] <- percentage
  annotation_p <- annotation_p + percentage
}
kingdom_df['No-annotation','abundance'] <- 1 - annotation_p

kingdom_df$Kingdom <- rownames(kingdom_df)
kingdom_df$Kingdom <- factor(kingdom_df$Kingdom)
kingdom_df$Kingdom <- fct_relevel(kingdom_df$Kingdom,rev(c('Viruses','Archaea','No-annotation','Bacteria'))) 

# Laber the RSA of each kingdom
label_value <- paste(round(kingdom_df$abundance * 100, 1), '%', sep = '')
label_value <- sort(label_value,decreasing = T)

ggplot(data = kingdom_df, aes(x ='Kingdom', y = abundance*100, fill = Kingdom)) + 
  geom_bar(stat = 'identity', position = 'stack', width = 1)+ 
  coord_polar(theta = 'y') + 
  labs(x = '', y = '', title = '') + 
  labs(fill="Kingdom (RSA%)")+
  ggtitle('RSA of annotated taxa (3461) in Tel Meggido (4259) \n - 81% annotated taxa count for 77% global RSA')+
  theme_bw()+
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank()) +
  scale_fill_discrete(labels = paste0(rev(c('Viruses','Archaea','No-annotation','Bacteria')),' (', c(label_value[1:3],'0.0006%'),")"))
  
## total pathogen RSA in Tel Megiddo samples
taxa_avg_abun <- data.frame(colSums(telm_abun_data)/dim(telm_abun_data)[1])

pathogen_abun <- data.frame(row.names = 'total')

pathogen_abun[1,'animal_pathogen_n'] <- sum(taxa_avg_abun[animal_pathogen,1] > 0)
pathogen_abun[1,'animal_pathogen_p'] <- sum(as.numeric(taxa_avg_abun[animal_pathogen,1]))
pathogen_abun[1,'plant_pathogen_n'] <- sum(taxa_avg_abun[plant_pathogen,1] > 0)
pathogen_abun[1,'plant_pathogen_p'] <- sum(as.numeric(taxa_avg_abun[plant_pathogen,1]))
pathogen_abun[1,'both_pathogen_n'] <- sum(taxa_avg_abun[both_pathogen,1] > 0)
pathogen_abun[1,'both_pathogen_p'] <- sum(as.numeric(taxa_avg_abun[both_pathogen,1]))

for (i in 1:dim(taxa_avg_abun)[1]){
  taxa <- rownames(taxa_avg_abun)[i]
  if (taxa %in% taxa_anno_file_telm2$taxa){
    taxa_avg_abun[i,'pathogen_type'] <- taxa_anno_file_telm2[taxa_anno_file_telm2$taxa==taxa,'pathogen_type']
  }else{
    taxa_avg_abun[i,'pathogen_type'] <- 'No_annotation'
  }
}
```
### 7.3 Microbiome diversity and composition detection in different variable groups  
Used the same code to analyze different groups of variables (site, surface material and access), changing only the input of the variable *variable*. The script was connected with the last part
```r
variable <- 'Site'
# variable <- 'Material'
# variable <- 'Access'

telm_file2[,variable] <- factor(telm_file2[,variable])

# Annotate the group of each sample
for (i in 1:length(rownames(telm_alpha_diversity))){
  id <- rownames(telm_alpha_diversity)[i]
  telm_alpha_diversity[i,'group'] <- telm_read_file2[telm_read_file2$sample_id == id, variable]
}

# Calculate the amonumt of smaples in each group
sample_n <- table(telm_alpha_diversity$group)

# Produce a comparision list in Wilcoxon test
mycompare <- list()
for (i in 1:(length(levels(factor(telm_alpha_diversity$group)))-1)){
  k1 <- levels(factor(telm_alpha_diversity$group))[i]
  for (j in (i+1):length(levels(factor(telm_alpha_diversity$group)))){
    k2 <- levels(factor(telm_alpha_diversity$group))[j]
    k <- c(k1,k2)
    mycompare <- c(mycompare,list(k))
  }
}

# Use kruskal test to check if the observed species is significantly different between different group
kruskal.test(Observed_species~group,data=telm_alpha_diversity)

# Plot observed species in each group
ggplot(telm_alpha_diversity,aes(x=group,y=Observed_species,color=group))+
  geom_boxplot(alpha=1,outlier.size = 1,size = 0.9,width=0.7,fill='transparent')+
  geom_jitter(position = position_jitter(0.1),size=3,alpha=0.6)+
  theme_classic()+
  labs(x=paste0('Tel Meggido ',variable),y='Observed species')+
  geom_signif(comparisons = mycompare,
              step_increase = 0.12,
              map_signif_level = T,
              test = wilcox.test)+
  theme(title = element_text(size=14))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text = element_text(colour = 'black',size = 9))+
  theme(legend.position="none")+
  scale_x_discrete(breaks = c(levels(factor(telm_alpha_diversity$group))),
                   labels = c(paste0(capitalize(levels(factor(telm_alpha_diversity$group))),' (',sample_n,')')))

# Use kruskal test to check if the Shannon index is significantly different between different group
kruskal.test(Shannon~group,data=telm_alpha_diversity)

# Plot Shannon index in each group
ggplot(telm_alpha_diversity,aes(x=group,y=Shannon,color=group))+
  geom_boxplot(alpha=1,outlier.size = 1,size = 0.9,width=0.7,fill='transparent')+
  geom_jitter(position = position_jitter(0.1),size=3,alpha=0.6)+
  theme_classic()+
  labs(x=paste0('Tel Meggido ',variable),y='Shannon ')+
  geom_signif(comparisons = mycompare,
              step_increase = 0.13,
              map_signif_level = T,
              test = wilcox.test)+
  theme(title = element_text(size=14))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text = element_text(colour = 'black',size = 9))+
  theme(legend.position="none")+
  scale_x_discrete(breaks = c(levels(factor(telm_alpha_diversity$group))),
                   labels = c(paste0(capitalize(levels(factor(telm_alpha_diversity$group))),' (',sample_n,')')))

# Use kruskal test to check if the Gini-Simpson is significantly different between different group
kruskal.test(Simpson~group,data=telm_alpha_diversity)

# Plot Gini-Simpson index in each group
ggplot(telm_alpha_diversity,aes(x=group,y=Simpson,color=group))+
  geom_boxplot(alpha=1,outlier.size = 1,size = 0.9,width=0.7,fill='transparent')+
  geom_jitter(position = position_jitter(0.1),size=3,alpha=0.6)+
  theme_classic()+
  labs(x=paste0('Tel Meggido ',variable),y='Simpson')+
  geom_signif(comparisons = mycompare,
              step_increase = 0.12,
              map_signif_level = T,
              # color = 'black',
              test = wilcox.test)+
  theme(title = element_text(size=14))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,angle = 45,hjust = 1))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text = element_text(colour = 'black',size = 9))+
  theme(legend.position="none")+
  scale_x_discrete(breaks = c(levels(factor(telm_alpha_diversity$group))),
                   labels = c(paste0(capitalize(levels(factor(telm_alpha_diversity$group))),' (',sample_n,')')))

## Beta diversity analysis

# Load packages used in beta diversity analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(amplicon)

# PCoA analysis
sample_group <- as.data.frame(telm_alpha_diversity$group)
colnames(sample_group) <- 'Group'
rownames(sample_group) <- rownames(telm_alpha_diversity)

pcoa_result <- BetaDiv(otu=t(telm_read_data),map=sample_group,group='Group',
                       dist='bray',method='PCoA',Micromet='adonis')

pcoa_points <- pcoa_result[[2]]
pcoa_method = 'PCoA'
pcoa_adonis <- pcoa_result[[5]]

# Re-plot the PCoA figure
ggplot(pcoa_points, aes(x=x, y=y, fill=Group)) +
  geom_point(alpha=.7, size=2, pch=21) +
  labs(x=paste0(pcoa_method," 1 (33.65%)"), # according to pcoa_result[[1]]
       y=paste0(pcoa_method," 2 (17.34%)"),
       title='adonis: R=0.031 p=0.245') + # according to pcoa_adonis
  stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))+
  guides(color=F)+
  theme_bw()+
  geom_hline(aes(yintercept=0), colour="black", linetype=2) +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_text(vjust=-8.5,hjust=0.97,size = 10),
        text=element_text(size=10),
        axis.text =element_text(colour="black",size=9),
        legend.text=element_text(size=9))+
  scale_fill_discrete(name = variable, 
                      labels =  c(paste0(capitalize(levels(factor(pcoa_points$Group))),' (',sample_n,')')))
# ggsave(paste0('telm_',variable,'_PCoA_Bray_Crutis.png'),p)

# NMDS analysis

nmds_result <- BetaDiv(otu=t(telm_read_data),map=sample_group,group='Group',
                       dist='bray',method='NMDS',Micromet='adonis')

nmds_points <- nmds_result[[2]]
nmds_method = 'NMDS'
nmds_adonis <- nmds_result[[5]]

# Re-plot the NMDS figure
ggplot(nmds_points, aes(x=x, y=y, fill=Group)) +
  geom_point(alpha=.7, size=2, pch=21) +
  labs(x=paste(nmds_method,"1", sep=""),
       y=paste(nmds_method,"2",sep=""),
       title='adonis: R=0.031 p=0.251\nStress=0.13')+ # according to pcoa_adonis and nmds_result[[1]]
  stat_ellipse( linetype=2,level=0.65,aes(group  =Group, colour= Group))+
  guides(color=F)+
  theme_bw()+
  geom_hline(aes(yintercept=0), colour="black", linetype=2) +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_text(vjust=-13.5,hjust=0.97,size = 10),
        text=element_text(size=10),
        axis.text =element_text(colour="black",size=9),
        legend.text=element_text(size=9))+
  scale_fill_discrete(name = variable, 
                      labels =  c(paste0(capitalize(levels(factor(nmds_points$Group))),' (',sample_n,')')))

## analysis microbiome composition  

telm_file2[,variable] <- factor(telm_file2[,variable])

# Calculate average RSA of each taxa in each group
abun_df <- data.frame(row.names = c(levels(telm_file2[,variable])))
for (i in 1:length(levels(telm_file2[,variable]))){
  for (j in 1: length(telm_taxa_list)){
    v <- levels(telm_file2[,variable])[i]
    taxa <- telm_taxa_list[j]  
    abun <- sum(telm_file2[telm_file2[,variable] == v,][,taxa])/length(telm_file2[telm_file2[,variable] == v,][,taxa])
    abun_df[v,taxa] <- abun
  }
}

sum(abun_df[1,])

# Phylum relative abundance
phylum_df <- data.frame(row.names = c(levels(telm_file2[,variable])))

for (i in 1:length(levels(telm_file2[,variable]))){
  city <- levels(telm_file2[,variable])[i]
  for (j in 1:length(levels(taxa_anno_file_telm$Phylum))){
    phylum <- levels(taxa_anno_file_telm$Phylum)[j]
    taxa_list <- taxa_anno_file_telm[taxa_anno_file_telm$Phylum == phylum,][,'taxa']
    # Calculate the sum of relative abundance values of taxa in the same phylum
    percentage <- sum(as.numeric(abun_df[city,taxa_list]))
    phylum_df[city,phylum] <- percentage
  }
  phylum_df[city,'No-annotation'] <- 1 - sum(phylum_df[city,levels(taxa_anno_file_telm$Phylum)])
}

# unique_phylum in each group
phylum_df2 <- phylum_df
phylum_df2[phylum_df2 > 0] <- 1
phylum_occupancy <- colSums(phylum_df2)
phylum_occupancy <- phylum_occupancy[order(phylum_occupancy)]

uniq_p_city <- phylum_df2[,c(names(phylum_occupancy[1:8]))]

phylum_average_abun <- colSums(phylum_df)/dim(phylum_df)[1]
phylum_average_abun <- phylum_average_abun[order(phylum_average_abun,decreasing = T)]

# Plot the most abundant phylum
select_phylum <- c('Actinobacteria',       
                   'Proteobacteria',
                   'Firmicutes',
                   'Bacteroidetes',
                   'Gemmatimonadetes',
                   'Myxococcota',
                   'Planctomycetes',
                   'No-annotation'
)

select_phylum_df <- phylum_df[,select_phylum]

for (i in 1:dim(select_phylum_df)[1]){
  select_phylum_df[i,'other'] <- sum(phylum_df[i,-which(colnames(phylum_df) %in% select_phylum)])
}

select_phylum_df <- select_phylum_df %>% dplyr::select(select_phylum[1:7], 'other','No-annotation')

select_phylum_df <- select_phylum_df[order(select_phylum_df$Actinobacteria),]

sample_n <- table(telm_file2[,variable])
sample_n2 <- sample_n[rownames(select_phylum_df)]

palette <- colorRampPalette(brewer.pal(9, "Set3"))(9) 

# Figure for site variable
# png("TelM_site_phylum_abundance.png", width = 6,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(5,0,0,8.2))
bp <- barplot(t(select_phylum_df*100), col=palette,
              # names.arg= ,
              args.legend = list(x = "bottom", inset=c(-0.5,0)), las =2,
              cex.names=.7,ylab = "", axisnames = F, axes = T, space =0,
              ylim = c(0,100))
mtext(text = paste0(capitalize(rownames(select_phylum_df)),' (',sample_n2,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 1)
title(ylab="Relative abundance % (Phylum level)",mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-0.65,0), colnames(select_phylum_df),  fill = palette , bty = 1, cex = 1)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()

# Figure for surface material variable
# png("TelM_material_phylum_abundance.png", width = 6,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(5,0,0,8.2))
bp <- barplot(t(select_phylum_df*100), col=palette,
              # names.arg= ,
              args.legend = list(x = "bottom", inset=c(-0.5,0)), las =2,
              cex.names=.7,ylab = "", axisnames = F, axes = T, space =0,
              ylim = c(0,100))
mtext(text = paste0(capitalize(rownames(select_phylum_df)),' (',sample_n2,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 1)
title(ylab="Relative abundance % (Phylum level)",mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-0.65,0), colnames(select_phylum_df),  fill = palette , bty = 1, cex = 1)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()


## taxa number in each phylum appeared in each group
phylum_taxa <- data.frame(row.names = c(levels(telm_file2[,variable])))

for (i in 1:length(levels(telm_file2[,variable]))){
  city <- levels(telm_file2[,variable])[i]
  for (j in 1:length(levels(taxa_anno_file_telm$Phylum))){
    phylum <- levels(taxa_anno_file_telm$Phylum)[j]
    taxa_list <- taxa_anno_file_telm[taxa_anno_file_telm$Phylum == phylum,][,'taxa']
    # Calculate the sum of relative abundance values of taxa in the same phylum
    num <- sum(abun_df[city,taxa_list]>0)
    phylum_taxa[city,phylum] <- num
  }
  
  phylum_taxa[city,'No-annotation'] <- sum(abun_df[city,]>0) - sum(phylum_taxa[city,levels(taxa_anno_file_telm$Phylum)])
  phylum_taxa[city,'Total_taxa_number'] <- sum(abun_df[city,]>0)
}


## Pathogen analysis 

taxa_anno_file_telm2 <- taxa_anno_file_telm
taxa_anno_file_telm2[is.na(taxa_anno_file_telm2)] <- 0

animal_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 1,'taxa']
plant_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 2,'taxa']
both_pathogen <- taxa_anno_file_telm2[taxa_anno_file_telm2$pathogen_type == 3,'taxa']

pathogen_abun <- data.frame(row.names = rownames(abun_df))

for (i in 1:dim(pathogen_abun)[1]){
  pathogen_abun[i,'animal_pathogen_n'] <- sum(abun_df[i,animal_pathogen] > 0)
  pathogen_abun[i,'animal_pathogen_p'] <- sum(as.numeric(abun_df[i,animal_pathogen]))
  pathogen_abun[i,'plant_pathogen_n'] <- sum(abun_df[i,plant_pathogen] > 0)
  pathogen_abun[i,'plant_pathogen_p'] <- sum(as.numeric(abun_df[i,plant_pathogen]))
  pathogen_abun[i,'both_pathogen_n'] <- sum(abun_df[i,both_pathogen] > 0)
  pathogen_abun[i,'both_pathogen_p'] <- sum(as.numeric(abun_df[i,both_pathogen]))
}

bar_df2 <- pathogen_abun[,c('plant_pathogen_p','both_pathogen_p','animal_pathogen_p')]

# Figure of pathogen abundance for site and material varibale
# png("pathogen_abundance_site.png", width = 6,height = 6, units = 'in', res = 600)
# png("pathogen_abundance_material.png", width = 6,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(5,0,0,8.2))
bp <- barplot(t(bar_df2*100), col=c("#619CFF","#00BA38","#F8766D"), 
              ylim = c(0,20),
              names.arg= paste0(capitalize((rownames(pathogen_abun))),' (',sample_n,')') ,
              args.legend = list(x = "right", inset=c(-0.5,0)), las =2, cex.names=.7,
              ylab = "",axisnames = F, axes = T)
mtext(text = paste0(capitalize((rownames(pathogen_abun))),' (',sample_n,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 1)
title(ylab="Pathogen relative abundance %",mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-0.65,0),
       c('Animal Pathogen','Dual-kingdom\nPathogen','Plant Pathogen'),
       fill = rev(c("#619CFF","#00BA38","#F8766D")) , bty = 1, cex = 1)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()

# Figure of pathogen abundance for access varibale
# png("pathogen_abundance_access.png", width = 4,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(3.5,0,0,8.5))
bp <- barplot(t(bar_df2*100), col=c("#619CFF","#00BA38","#F8766D"), 
              ylim = c(0,10),
              names.arg= paste0(capitalize((rownames(pathogen_abun))),' (',sample_n,')') ,
              args.legend = list(x = "right", inset=c(-0.5,0)), las =2, cex.names=.7,
              ylab = "",axisnames = F, axes = T)
mtext(text = paste0(capitalize((rownames(pathogen_abun))),' (',sample_n,')'), side = 1, at = bp, line = 0.5, padj = 1, cex = 1)
title(ylab="Pathogen relative abundance %",mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-1.95,0),
       c('Animal Pathogen','Dual-kingdom\nPathogen','Plant Pathogen'),
       fill = rev(c("#619CFF","#00BA38","#F8766D")) , bty = 1, cex = 1)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()
```

### 7.4 The mGPS prediction for public and non-public samples in Tel Megiddo

```r
library(readxl)
library(rworldmap)
library(maps)

pre_file <- read.csv('Prediction_results.csv')
anno_file <- read_excel('TelMegiddo - sampling_final.xlsx',sheet = 1, na = 'NA')

merge_file <- pre_file[,c('pangea_barcode','longitude','latitude',"cityPred","latPred","longPred",'Distance_from_origin')]

for (i in 1:length(merge_file$pangea_barcode)){
  s_id <- merge_file$pangea_barcode[i]
  s_id <- strsplit(s_id,'-')[[1]][2]
  s_id <- substr(s_id,2,length(strsplit(s_id,'')[[1]]))
  merge_file[i,'Access'] <- anno_file[anno_file$number == s_id,'Access']
}
merge_file$Access <- factor(merge_file$Access)

# png('TelMegiddo_access_prediction.png', width = 13,height = 8, units = 'in', res = 600)
par(mai=c(2,1,0.5,0.5), mar=par()$mar+c(3,0,0,0))
map <- rworldmap::getMap(resolution = "coarse")
palette <-c( "darkorchid4","gold2","dodgerblue3","brown","orangered2")
plot(map, xlim = c(-75,80),ylim = c(-15,50), col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.6)
#find coord preds by region
for ( i in 1:length(levels(merge_file$Access))){
  this_continent <- levels(merge_file$Access)[i]
  find_lats <- merge_file[merge_file[,"Access"] == this_continent,][,"latPred"]
  find_longs <- merge_file[merge_file[,"Access"] == this_continent,][,"longPred"]
  #plot predicted co-ordinates
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 2)
}
legend(55,0, c('non-public','public'), pch = "+", col = palette, cex = 1.5, bg ="lightskyblue1",ncol = 1)
map.axes(cex.axis = 1.4)
par(mar=c(5, 4, 4, 2) + 0.1)
# dev.off()

# Accuracy bar plot
bar_df1 <- data.frame(row.names = c("Overall",levels(merge_file$Access)))
for (i in 1: length(levels(merge_file$Access))){
  this_city <- levels(merge_file$Access)[i]
  prop <- mean(merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] < 100)
  bar_df1[i+1,"0 - 100km"] <- prop
  
  overall_prop <- mean(merge_file[,"Distance_from_origin"] < 100)
  bar_df1[ 1,"0 - 100km"] <- overall_prop
}
for (i in 1: length(levels(merge_file$Access))){
  this_city <- levels(merge_file$Access)[i]
  prop <- mean(merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] > 100 & merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] < 500)
  bar_df1[i+1,"100 - 500km"] <- prop
  
  overall_prop <-mean(merge_file[,"Distance_from_origin"] > 100 & merge_file[,"Distance_from_origin"] < 500)
  bar_df1[ 1,"100 - 500km"] <- overall_prop
}
for (i in 1: length(levels(merge_file$Access))){
  this_city <- levels(merge_file$Access)[i]
  prop <- mean(merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] > 500 & merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] < 1000)
  bar_df1[i+1,"500 - 1000km"] <- prop
  
  overall_prop <- mean(merge_file[,"Distance_from_origin"] > 500 & merge_file[,"Distance_from_origin"] < 1000)
  bar_df1[ 1,"500 - 1000km"] <- overall_prop
}
for (i in 1: length(levels(merge_file$Access))){
  this_city <- levels(merge_file$Access)[i]
  prop <- mean(merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] > 1000 & merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] < 2000)
  bar_df1[i+1,"1000 - 2000km"] <- prop
  
  overall_prop <- mean(merge_file[,"Distance_from_origin"] > 1000 & merge_file[,"Distance_from_origin"] < 2000)
  bar_df1[ 1,"1000 - 2000km"] <- overall_prop
}
for (i in 1: length(levels(merge_file$Access))){
  this_city <- levels(merge_file$Access)[i]
  prop <- mean(merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] > 2000 & merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] < 3000)
  bar_df1[i+1,"2000 - 3000km"] <- prop
  
  overall_prop <- mean(merge_file[,"Distance_from_origin"] > 2000 & merge_file[,"Distance_from_origin"] < 3000)
  bar_df1[1,"2000 - 3000km"] <- overall_prop
}
for (i in 1: length(levels(merge_file$Access))){
  this_city <- levels(merge_file$Access)[i]
  prop <- mean(merge_file[merge_file$Access == this_city,][,"Distance_from_origin"] > 3000 )
  bar_df1[i+1,"> 3000km"] <- prop
  
  overall_prop <- mean(merge_file[,"Distance_from_origin"] > 3000)
  bar_df1[ 1,"> 3000km"] <- overall_prop
}
size <- c()
for (i in 1: length(levels(merge_file$Access))){
  
  this_city <- levels(merge_file$Access)[i]
  size[i] <- length(which(merge_file$Access == this_city))
}

# png('accuracy_tel.png', width = 4,height = 6, units = 'in', res = 600)
par(xpd = T, mgp = c(0,0.7,0), las=2,mar = par()$mar + c(3.5,0,0,8.5))
bp <- barplot(t(bar_df1*100), space = 0,col=c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), 
names.arg=c("Overall",paste0(levels(merge_file$Access),"  (",size,")"), axes = F) , 
              las =2, cex.names=1.2, ylab = "", axisnames = F, axes = F)
axis(side =2, pos = 0)
mtext(text = c("Overall",paste0(levels(merge_file$Access)," (",size,")")), side = 1, at = bp, line = 0.5, padj = 1, cex = 1.2)
title(ylab="Proportion of sample predictions %", mgp=c(2,0,0),cex.lab=1.2)
legend("right",inset = c(-1.95,0), rev(c(colnames(bar_df1))), 
       fill = rev(c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) ,
       bty = 1, cex = 1.2)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()
```
