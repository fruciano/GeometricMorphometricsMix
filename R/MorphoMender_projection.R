MorphoMender_projection=function(data, models, method=c("Fruciano2020", "Valentin2008")){
  
}

model_input_validation=function(models){
  # Models should be provided as a list
  if(!is.list(models)){
    stop("models must be a list")
  }
  # If the list has a single element, return a warning that usually this is not recommended
  if(length(models)==1){
    warning("A single model has been provided, this is usually not recommended")
  }
  # Check the object type of each element in the list to see whether they are
  # three dimensional arrays, data frames or matrices. If they are none of these
  # return an error
  for(i in 1:length(models)){
    if(!is.array(models[[i]]) & !is.data.frame(models[[i]]) & !is.matrix(models[[i]])){
      stop("Each element in the list must be a three dimensional array, a data frame or a matrix")
    }
  }
  
  # In case it is an array, make sure that there is a three-dimensional array
  if(is.array(models[[1]])){
    if(dim(models[[1]])[3]!=3){
      stop("The array must have three dimensions:
           rows for landmarks, columns for x, y, (z) coordinates,
           and a third dimension for observations")
    }
  }
}


models=list(NAM="NAM", NAM2="NAM2")
