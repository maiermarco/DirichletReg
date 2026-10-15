update.DirichletRegModel <- function(
  object,         # an object, in this case of class "DirichletRegModel"
  formula.,       # a formula object that specifies the updates
  ...,            # extra arguments
  evaluate = TRUE # logical value that indicates whether to evaluate the updated call or not
){
  
  if(!inherits(object, "DirichletRegModel")) stop("\"object\" must be an instance of class \"DirichletRegModel\"") # error if object's class is not DirichletRegModel
  
  object_call <- getCall(object)
  if(is.null(object_call)) stop("\"object\" must have a \"call\" component") # error if object$call is NULL
  extras <- match.call(expand.dots = FALSE)$"..." # collect extra arguments
  
  if(!missing(formula.)){ # if formula. is specified
    if(!inherits(formula., c("formula", "Formula"))) stop("\"formula.\" must be an instance of class \"formula\" or \"Formula\"") # error if formula.'s class is neither "formula" nor "Formula"
    object_call$formula <- formula(update(Formula(formula(object)), formula.)) # replace the call's formula element with the updated formula
  }
  
  if(length(extras)){ # if any extra arguments were specified
    existing <- names(extras) %in% names(object_call) # elements in extras that already exists in  
    for(a in names(extras)[existing]) object_call[[a]] <- extras[[a]] # replace the already existing elements in object_call with new ones from extras
    if(!all(existing)){ # if there are any elements in extras that do not exists in object_call
      object_call <- as.call(c(as.list(object_call), extras[!existing])) # convert object_call to a list, append the remaining extras and convert everything in to a call again
    }
  }
  
  #
  # TO DO: Extract coefficients from the old model and use them as starting values for the updated one
  #
  
  # if evaluate is true, evaluate object_call (in the parent frame), else return the updated call
  if(evaluate){ eval(object_call, parent.frame()) } else { object_call }
}
