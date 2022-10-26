#' Bayes factors and posterior probabilities for linear regression models
#'
#' It Computes the Bayes factors and posterior probabilities of a list of linear
#' regression models proposed to explain a common response variable over the
#' same dataset
#'
#' The Bayes factors, BFi0, are expressed in relation with the simplest model
#' (the one nested in all the others). Then, the posterior probabilities of the
#' entertained models are obtained as
#'
#' Pr(Mi | \code{data})=Pr(Mi)*BFi0/C,
#'
#' where Pr(Mi) is the prior probability of model Mi and C is the normalizing
#' constant.
#'
#' The Bayes factor B_i depends on the prior assigned for the parameters in the regression
#'
#' Option is "eZellnerSiow" a prior which 
#'
#' With respect to the prior over the model space Pr(Mi) three possibilities
#' are implemented: "Constant", under which every model has the same prior
#' probability and "User". With this last option, the prior probabilities are
#' defined through the named list \code{priorprobs}. These probabilities can be
#' given unnormalized.
#'
#'
#' @aliases Btest
#' @export
#' @param models A named list with the entertained models defined with their
#' corresponding formulas. If the list is unnamed, default names are given by
#' the routine. One model must be nested in all the others.
#' @param data data frame containing the data.
#' @param prior.betas Prior distribution for regression parameters within each
#' model (to be literally specified). 
#' @param prior.models Type of prior probabilities of the models (to be literally specified). Possible choices are
#' "Constant" and "User" (see details).
#' @param priorprobs A named vector ir list (same length and names as in argument
#' \code{models}) with the prior probabilities of the models (used in combination
#' of \code{prior.models="User"}). If the provided object
#' is not named, then the order in the list of \code{models} is used to assign the prior
#' probabilities
#' @param null.model The name of the null model (eg. the one nested in all the others).
#' By default, the names of covariates in the different
#' models are used to identify the null model. An error is produced if such identification fails.
#' This identification
#' is not performed if the definition of the null model is provided, with this argument,
#' by the user.
#' Note that the (the \code{null.model} must coincide with that model with the largest
#' sum of squared errors and should be smaller in dimension to any other model).
#' @return \code{Btest} returns an object of type \code{Btest} which is a
#' \code{list} with the following elements: \item{BFio }{A vector with the
#' Bayes factor of each model to the simplest model.} \item{PostProbi }{A
#' vector with the posterior probabilities of each model.} \item{models }{A
#' list with the entertained models. } \item{nullmodel}{The position of the
#' simplest model. }
#' \item{prior.betas}{prior.betas}
#' \item{g}{The value of g}
#' \item{prior.models}{prior.models}
#' \item{priorprobs}{priorprobs}

#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: <anabel.forte@@uv.es>
#' @seealso \code{\link[BayesVarSel]{Bvs}} for variable selection within linear
#' regression models
#' @references Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.
#'
Btest <-
  function(models,
           data,
           prior.betas = "eZellnerSiow",
		   g = dim(data)[1],
           prior.models = "Constant",
           priorprobs = NULL,
		   null.model = NULL) {

	#N is the number of models:
  N <- length(models)

	if (!is.list(models)) stop("Argument models should be a list\n")

	#If competing models come wihtout a name, give one by default:
	if (is.null(names(models))){
		if (!is.null(null.model)) stop("Please provide a name for the competing models. The null model must be in that list\n")
		names(models) <- paste("model", 1:N, sep="")
	}

	#Check if the given null model is one of the competing models:
	if (!is.null(null.model)){
	  relax.nest = TRUE
		pos.user.null.model <- which(null.model == names(models))
		if (length(pos.user.null.model) == 0) stop("The null model provided is not in the list of competing models.\n")
	}
  else relax.nest = FALSE


    #n is the sample size
    n <- dim(data)[1]
    #SSE is a vector with SSE's for each model; Dim with the dimension (number of regressors in each)
    SSE <- rep(0, N)
    Dim <- rep(0, N)
    BFi0 <- rep(0, N)
    PostProbi <- rep(0, N)

    #prior for betas: (#pfb will contain the internal representation of the prior.betas argument)
	#pfb <- check.prior.betas(prior.betas) 
	
    #prior for model space:
    pfms <- substr(tolower(prior.models), 1, 1)
    if (pfms != "c" &&
        pfms != "u")
      stop("I am very sorry: prior for models not supported\n")
    if (pfms == "u" &&
        is.null(priorprobs)) {
      stop("A valid vector of prior probabilities must be provided\n")
    }
    if (pfms == "u" &&
        length(priorprobs) != N) {
      stop("Vector of prior probabilities with incorrect length\n")
    }
    if (pfms == "u" &&
        sum(priorprobs < 0) > 0) {
      stop("Prior probabilities must be positive\n")
    }
    #If priorprobs are unnamed, give this vector the name of the models
    if (pfms == "u" &&
        is.null(names(priorprobs))) {
          names(priorprobs)<- names(models)
    }

    #Prior probabilities of models:
    PriorModels <- rep(0, N)
    if (prior.models == "Constant") {
      PriorModels <- rep(1, N)
    }
    if (prior.models == "User") {
      #should coincide with the length of prior.models
      for (i in 1:N) {
        PriorModels[i] <- priorprobs[[names(models)[i]]]
      }
    }

    #list that contains the names of the covariates in each model
    covar.list <- list()
    for (i in 1:N) {
      temp <-
        lm(
          formula = as.formula(models[[i]]),
          data = data,
          y = TRUE,
          x = TRUE
        )
      SSE[i] <- sum(temp$residuals ^ 2)
      Dim[i] <- length(temp$coefficients)
      covar.list[[i]] <- dimnames(temp$x)[[2]]
    }
    ordered.SSE <- sort(SSE, index.return = TRUE, decreasing = TRUE)
    #Which acts as null model:
    nullmodel <- ordered.SSE$ix[1]

	#Check that the given null model and the model with largest SSE
	#(which in turns defines the null model)
	#coincides
	if (!is.null(null.model)){
		if (nullmodel != pos.user.null.model){
		stop("The given null model does not coincide with the one with the largest sum of squared error (and it should).\n")
		}
	}

    for (i in (1:N)[-nullmodel]) {
        #check if the "null" model is nested in all the others
        if (!relax.nest &
            sum(covar.list[[nullmodel]] %in% covar.list[[i]]) < Dim[nullmodel]) {
            stop("I suspect that perhaps the simplest (null) model is not nested in all the others. Define explicitly the simplest model if you are sure it is the case\n.\n")
        }
        Qi0 <- SSE[i] / SSE[nullmodel]
        #The .C to be used:
        if (prior.betas == "eZellnerSiow")
          BFi0[i] <-
          .C(
            "eZSBF",
			as.double(g),
            as.integer(n),
            as.integer(Dim[i]),
            as.integer(Dim[nullmodel]),
            as.double(Qi0),
            as.double(0.0)
          )[6][[1]]
		  else stop("Prior not supported")
        }

    BFi0[nullmodel] <- 1
    names(BFi0) <-
      paste(names(models), ".to.", names(models)[nullmodel], sep = "")
    PostProbi <- BFi0 * PriorModels / sum(BFi0 * PriorModels)
    names(PostProbi) <- names(models)
    result <- list()
    result$BFi0 <- BFi0
    result$PostProbi <- PostProbi
    result$models <- models
    result$nullmodel <- nullmodel
	
	result$prior.betas<- prior.betas	
	result$prior.models<- prior.models
	result$priorprobs<- priorprobs
	
    class(result) <- "Btest"
    result
  }



#' Print an object of class \code{Btest}
#'
#'Print an object of class \code{Btest}
#' @export
#' @param x Object of class Btest
#' @param ... Additional parameters to be passed
#'
#' @seealso See \code{\link[BayesVarSel]{Btest}} for creating objects of the class \code{Btest}.
#'
#'
#' @examples
#'  \dontrun{
#' #Analysis of Crime Data
#' #load data
#' data(UScrime)
#' #Model selection among the following models: (note model1 is nested in all the others)
#' model1<- y ~ 1 + Prob
#' model2<- y ~ 1 + Prob + Time
#' model3<- y ~ 1 + Prob + Po1 + Po2
#' model4<- y ~ 1 + Prob + So
#' model5<- y ~ .
#'
#' #Equal prior probabilities for models:
#' crime.BF<- Btest(models=list(basemodel=model1,
#' 	ProbTimemodel=model2, ProbPolmodel=model3,
#' 	ProbSomodel=model4, fullmodel=model5), data=UScrime)
#' 	crime.BF
#' 	}
print.Btest <- function(x, ...){
  cat("---------\n")
  cat("Models:\n")
  print(x$models)
  cat("---------\n")
  cat(paste("Bayes factors (expressed in relation to ",names(x$models)[x$nullmodel],")\n", sep=""))
  print(x$BFi0)
  cat("---------\n")
  cat("Posterior probabilities:\n")
  print(round(x$PostProbi,3))
}
