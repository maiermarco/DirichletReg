update.DirichletRegModel <- function(
  object,         # an object, in this case of class "DirichletRegModel"
  formula.,       # a formula object that specifies the updates
  ...,            # extra arguments
  evaluate = TRUE # logical value that indicates whether to evaluate the updated call or not
){
  if(!inherits(object, "DirichletRegModel")) stop('"object" must be an instance of class "DirichletRegModel"') # error if object's class is not DirichletRegModel
  if(is.null(call <- getCall(object))) stop('"object" must have a "call" component') # error if object$call is NULL
  extras <- match.call(expand.dots = FALSE)$"..." # collect extra arguments
  
  if(!missing(formula.)){ # if formula. is specified
    if(!inherits(formula., c("formula", "Formula"))) stop('"formula." must be an instance of class "formula" or "Formula"') # error if formula.'s class is neither "formula" nor "Formula"
    call$formula <- formula(update(Formula(formula(object)), formula.)) # replace the call's formula element with the updated formula
  }
  
  if(length(extras)){ # if any extra arguments were specified
    existing <- names(extras) %in% names(call) # elements in extras that already exists in  
    for(a in names(extras)[existing]) call[[a]] <- extras[[a]] # replace the already existing elements in call with new ones from extras
    if(any(!existing)){ # if there are any elements in extras that do not exists in call
      call <- as.call(c(as.list(call), extras[!existing])) # convert the call to a list, append the remaining extras and convert everything in to a call again
    }
  }
  
  #
  # TO DO: Extract coefficients from the old model and use them as starting values for the updated one
  #
  
  # if evaluate is true, evaluate the call (in the parent frame), else return the updated call
  if(evaluate){ eval(call, parent.frame()) } else { call }
}
