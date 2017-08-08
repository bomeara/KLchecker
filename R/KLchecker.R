#' Function to calculate expected KL
#'
#' Based on Jhwueng et al. 2014, especially equation 7
#' @param simulator.true.model Function to create a dataset given vector of true parameters
#' @param true.parameters Vector of parameters to use for simulations
#' @param param.estimator.candidate.model Function to estimate the parameters under a candidate model given an input dataset
#' @param simulator.candidate.model Function to simulate a dataset given vector of estimated parameters
#' @param likelihood.data.with.true.model Function to calculate the log likelihood (not neg log like) of the data under the true model
#' @param likelihoods.data.with.candidate.model Function to calculate the log likelihood (not neg log like) of the data under the candidate model
#' @param K Number of free parameters in candidate model
#' @param nrep.outer Number of outer loops
#' @param nrep.inner Number of inner loops
#' @return KL in same space (so multiplied by -2) as AIC
#' @examples
#' simulator.true.model <- function(parameters) {
#'   return(rnorm(n=1000, mean=parameters[1], sd=parameters[2]))
#' }
#'
#' true.parameters <- c(1000, 0.12)
#'
#' param.estimator.candidate.model <- function(x) {
#'   return(mean(1/x)) #here our candidate is exp
#' }
#'
#' simulator.candidate.model <- function(parameters) {
#'   return(rexp(n=1000, rate=parameters[1]))
#' }
#'
#' likelihood.data.with.true.model <- function(data, parameters) {
#'   return(sum(dnorm(data, mean=parameters[1], sd=parameters[2]), log=TRUE))
#' }
#'
#' likelihood.data.with.candidate.model <- function(data) {
#'   return(sum(dexp(data, rate=param.estimator.candidate.model(data)), log=TRUE))
#' }
#'
#' K <- 1
#'
#' NormVsExp <- EstimateKL(simulator.true.model, true.parameters, param.estimator.candidate.model, simulator.candidate.model, likelihood.data.with.true.model, likelihoods.data.with.candidate.model, K)
#'
#'
#' simulator.true.model <- function(parameters) {
#'   return(rnorm(n=1000, mean=parameters[1], sd=parameters[2]))
#' }
#'
#' true.parameters <- c(1000, 0.12)
#'
#' param.estimator.candidate.model <- function(x) {
#'   return(c(mean(x), sd(x))) #here our candidate is exp
#' }
#'
#' simulator.candidate.model <- function(parameters) {
#'   return(rnorm(n=1000, mean=parameters[1], sd=parameters[2]))
#' }
#'
#' likelihood.data.with.true.model <- function(data, parameters) {
#'   return(sum(dnorm(data, mean=parameters[1], sd=parameters[2]), log=TRUE))
#' }
#'
#' likelihood.data.with.candidate.model <- function(data) {
#'   return(sum(dnorm(data, mean=mean(data), sd=sd(data)), log=TRUE))
#' }
#'
#' K <- 2
#'
#' NormVsNorm <- EstimateKL(simulator.true.model, true.parameters, param.estimator.candidate.model, simulator.candidate.model, likelihood.data.with.true.model, likelihoods.data.with.candidate.model, K)
#' @export
EstimateKL <- function(simulator.true.model, true.parameters, param.estimator.candidate.model, simulator.candidate.model, likelihood.data.with.true.model, likelihoods.data.with.candidate.model, K, nrep.outer=100, nrep.inner=100) {
    KL.vector <- rep(NA, nrep.outer)
    AIC.vector <- rep(NA, nrep.outer)
    for (outer.index in sequence(nrep.outer)) {
        original.data <- simulator.true.model(true.parameters)
        parameters <- param.estimator.candidate.model(original.data)
        lnL.vector <- rep(NA, nrep.inner)
        for(inner.index in sequence(nrep.inner)) {
            simulated.data <- simulator.candidate.model(parameters)
            lnL.vector[inner.index] <- likelihood.data.with.true.model(simulated.data, true.parameters)
        }
        KL.vector[outer.index] <- mean(lnL.vector)
        AIC.vector[outer.index] <- -2*likelihood.data.with.candidate.model(original.data)+2*K
    }
    result <- c(KL=-2*mean(KL.vector), AIC=mean(AIC.vector))
    return(result)
}
