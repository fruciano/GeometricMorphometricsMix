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
  # Check that all elements on the list are of the same type
  if(length(unique(sapply(models, class)))>1){
    stop("All elements in the list must be of the same type")
  }

  # In case the first model is an array, make sure that there is a three-dimensional array
  # and that all elements are three-dimensional arrays and have the same number of rows
  # and columns
  if(is.array(models[[1]])){
    model_format="landmark_array"
    if(dim(models[[1]])[3]!=3){
      stop("If arrays, models must be three-dimensional arrays")
    }
    for(i in 2:length(models)){
      if(dim(models[[i]])[3]!=3){
        stop("All models must be three-dimensional arrays")
      }
      if(dim(models[[i]])[1]!=dim(models[[1]])[1] | dim(models[[i]])[2]!=dim(models[[1]])[2]){
        stop("All models must have the same number of rows (landmarks) and columns (coordinates)")
      }
    }
    # If any of the models is a three dimensional array but has a number of columns
    # different than 2 or 3 stop returning an error
    if (any((sapply(models, function(x) ((dim(x)[2]) %in% c(2,3))==FALSE)))){
      stop("If models are threedimensional arrays, models must have 2 or 3 columns,
           corresponding to 2D or 3D landmark data")
    }
    # If the first model has two columns, return a message that the data is 2D
    if(dim(models[[1]])[2]==2){
      message("2D landmark data has been detected for the model")
      landmark_dims=2
    } else if (dim(models[[1]])[2]==3){
      message("3D landmark data has been detected for the model")
      landmark_dims=3
    }
    data_dims=dim(models[[1]])[1]*landmark_dims
  }
  # Otherwise if each model is a matrix or a data frame, make sure that they have the same number of columns
  # and that they have the same number of rows
  else if(is.matrix(models[[1]]) | is.data.frame(models[[1]])){
    model_format="general_data"
    data_dims=ncol(models[[1]])
    landmark_dims=0
    for(i in 2:length(models)){
      if(ncol(models[[i]])!=ncol(models[[1]])){
        stop("All models must have the same number of columns")
      }
    }
  }
  return(list(model_format=model_format, data_dims=data_dims, landmark_dims=landmark_dims))
}


models=list(NAM=array(NA, dim=c(4,3,1)), NAM2=array(NA, dim=c(4,3,1)))
