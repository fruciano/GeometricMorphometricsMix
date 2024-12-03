MorphoMender_projection=function(data_to_correct, models=NULL, user_defined_vec=NULL,
                                 model_align_method=c("Fruciano2020", "Valentin2008")){
  valid_model= model_input_validation(models)
  valid_data= data_input_validation(data_to_correct)
  model_align_method=model_align_method[1]
  # Validation for the data and models provided

  # In case both models and data to correct are landmarks,
  # perform alignment according to the choice set by the attribute model_align_method
  if (valid_data$data_format=="landmark_array" && valid_model$model_format=="landmark_array"){
    if (model_align_method=="Fruciano2020"){
      aligned_data_models=dt_model_align_OPA(data_to_correct, models)
    } else if (model_align_method=="Valentin2008"){
      aligned_data_models=dt_model_align_GPA(data_to_correct, models)
    }
  }


}

# Function to perform input validation for the models
# Returns messages and warnings to help the user or errors if the input is not correct
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
    if(length(dim(models[[1]]))!=3){
      stop("If arrays, models must be three-dimensional arrays")
    }
    for(i in 2:length(models)){
      if(length(dim(models[[1]]))!=3){
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

    # If the models are matrices or data frames, return a message that
    # the models are not recognised as landmark data and will be treated
    # as general models
    message("Models are not recognised as landmark data
            and will be treated as general models")

    # Check that the number of columns is the same for all models
    for(i in 2:length(models)){
      if(ncol(models[[i]])!=ncol(models[[1]])){
        stop("All models must have the same number of columns")
      }
    }
  }
  return(list(model_format=model_format, data_dims=data_dims, landmark_dims=landmark_dims))
}

# Function to perform input validation for the data to correct
# it returns messages and warnings to help the user or errors if the input is not correct
data_input_validation=function(data_to_correct){
  # Data should be provided as a three-dimensional array, a data frame or a matrix
  if(!is.array(data_to_correct) & !is.data.frame(data_to_correct) & !is.matrix(data_to_correct)){
    stop("data must be a three-dimensional array, a data frame or a matrix")
  }
  # If the data is a three-dimensional array, make sure that it has three dimensions
  if(is.array(data_to_correct)){
    data_format="landmark_array"
    if(length(dim(models[[1]]))!=3){
      stop("If arrays, data must be three-dimensional arrays")
    }
  # If the data
    # If the first model has two columns, return a message that the data is 2D
    if(dim(data_to_correct)[2]==2){
      message("2D landmark data has been detected for the data to be corrected")
      landmark_dims=2
    } else if (dim(data_to_correct)[2]==3){
      message("3D landmark data has been detected for the data to be corrected")
      landmark_dims=3
    }
    data_dims=dim(data_to_correct)[1]*landmark_dims
  } else if (is.matrix(data_to_correct) | is.data.frame(data_to_correct)){

    data_format="general_data"
    data_dims=ncol(data_to_correct)
    landmark_dims=0

    # If the data to correct are in a matrix or data frames, return a message that
    # the data are not recognised as landmark data and will be treated
    # as general models
    message("Data to correct not recognised as landmark data
            and will be treated as general data")

  }
  return(list(data_format=data_format, data_dims=data_dims, landmark_dims=landmark_dims))
  }



# Funtion to perform alignment using OPA
# Following Fruciano et al 2020 - ZJLS
# To align the models to the consensus of the data to correct
# So that models do not influence the alignment of data
#' @importFrom shapes procOPA procGPA
dt_model_align_OPA=function(Data_to_correct, Models_Landmarks){

  n_land=ncol(Data_to_correct[,,1])
  n_models=length(Models_Landmarks)
  # Number of models

  # First align each model and the dataset to correct separately
  Models_Landmarks_GPA=lapply(Models_Landmarks,procGPA)
  Data_to_correct_GPA=procGPA(Data_to_correct)

  Rotation_matrices_OPA_meanmodels_to_meandata=lapply(Models_Landmarks_GPA, function(X)
    shapes::procOPA(A=Data_to_correct_GPA$mshape,
                    B=X$mshape)$R )
  # Compute rotation matrices to align the consensus of each of the models
  # to the consensus of the data to correct

  # Now apply the rotation to the models to align them to the
  # consensus of the dataset to correct
  Models_Landmarks_GPA_Rotated=lapply(seq_len(n_models), function(i) {
    rotated_models=array(NA, dim = dim(Models_Landmarks_GPA[[i]]$rotated))
    for (sp in seq_len(dim(rotated_models)[3])) {
      rotated_models[,,sp]=Models_Landmarks_GPA[[i]]$rotated[,,sp]%*%
        Rotation_matrices_OPA_meanmodels_to_meandata[[i]]
    }
    return(rotated_models)
  })
return(list(rotated_data=Data_to_correct_GPA$rotated,
            rotated_models=Models_Landmarks_GPA_Rotated))
}


# Function to perform alignment using GPA following
# Valentin et al 2008 - Journal of Fish Biology
# This aligns both data and models with a single GPA
# Usually, not the preferred option as in principle the models
# will affect the GPA consensus and alignment of the dataset to correct
dt_model_align_GPA=function(Data_to_correct, Models_Landmarks){

  n_models=length(Models_Landmarks)
  # Number of models

  n_obs_data=dim(Data_to_correct)[3]
  n_obs_models=sapply(Models_Landmarks, function(X) dim(X)[3])
  # Number of observations in the data to correct and in the models

  combined_dt=array(NA, dim=c(dim(Data_to_correct)[seq(2)],n_obs_data+sum(n_obs_models)))
  # Pre-allocate array

  combined_idxs=data.frame(type=c(rep("data", n_obs_data),
                                  unlist(lapply(seq(n_obs_models), function(i){
                                    rep(paste0("model_",i), n_obs_models[i])
                                  }))))
  combined_idxs$idx=seq(nrow(combined_idxs))


  combined_dt[,,seq(n_obs_data)]=Data_to_correct
  # Add the data to correct to the first portion of the combined array

  for(i in seq(n_models)){
    tmp_idxs=combined_idxs[combined_idxs$type==paste0("model_",i),]$idx
    combined_dt[,,tmp_idxs]=Models_Landmarks[[i]]
  }
  # Add the models to the second portion of the combined array

  # Perform GPA on the combined dataset
  combined_dt_GPA=procGPA(combined_dt)
  rotated_data=combined_dt_GPA$rotated[,,seq(n_obs_data)]
  rotated_models=lapply(seq(n_models), function(i) {
    tmp_idxs=combined_idxs[combined_idxs$type==paste0("model_",i),]$idx
    return(combined_dt_GPA$rotated[,,tmp_idxs])
  })
  return(list(rotated_data=rotated_data,
              rotated_models=rotated_models))
}

