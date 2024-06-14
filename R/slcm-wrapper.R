# Construct the SLCM object
new_slcm_model = function(model_mcmc,
                          details,
                          estimates = NULL) {

    estimate_chain = lapply(model_mcmc, summarize_model)
    estimates = c(estimate_chain, estimates)
    chain = model_mcmc[!(names(model_mcmc) %in% estimates)]
    
    
    # Fix binary and order
    attribute_description = attribute_pattern_table_header(details$k, 2, details$k)
    colnames(estimates$beta$mean) <- paste0("B_", attribute_description)
    rownames(estimates$beta$mean) <- paste0("Item", seq_len(details$j))
    
    colnames(estimates$theta$mean) <- paste0("Theta_", attribute_description)
    rownames(estimates$theta$mean) <- paste0("Item", seq_len(details$j))
    
    colnames(estimates$delta) <- paste0("D_", attribute_description)
    rownames(estimates$delta) <- paste0("Item", seq_len(details$j))
    
    colnames(estimates$q) <- paste0("Q_", seq_len(details$k))
    rownames(estimates$q) <- paste0("Item", seq_len(details$j))
    
    names(estimates$pi$mean) <- attribute_description
    
    structure(
        list(
            estimates = estimates, # Iterates over each parameter and obtains summary information
            chain = chain,
            details = details
        ),
        class = c("slcm", "edm")
    )
}

#' Print the SLCM object
#' 
#' Custom printing class to reveal features of the fitted SLCM.
#' 
#' @param x      the `slcm` object.
#' @param digits the number of significant digits
#' @param ...    further arguments passed to or from other methods.
#' 
#' @returns 
#' Print details and estimates found within the fitted SLCM. 
#' Return the model invisibly (via [`invisible()`])
#' 
#' @export
print.slcm = function(x, digits = max(3L, getOption("digits") - 3L), ...) { 
  details = x$details
  
  
  cat("\nModel Details:\n",
      "- Observations (n): ", details$n, "\n",
      "- Items (j): ", details$j, "\n",
      "- Attributes (k): ", details$k, "\n",
      "- Runtime: ", details$runtime, "\n",
      "- Date: ", details$date_time, "\n",
      "- Package Version: ", details$package_version, "\n",
      "\n", sep = "")

  cat("Chain properties:\n",
      "- Burn in: ", details$burnin, "\n",
      "- Chain Length: ", details$chain_length, "\n",
      "- Total Iterations: ", details$chain_length + details$burnin, "\n",
      "\n", sep = "")  
  
  cat("Hyperparameter Details:\n",
      "- m0: ", details$n, "\n",
      "- bq: ", details$j, "\n",
      "- l1:\n", sep = "")
  print.default(format(t(details$l1), digits = digits), print.gap = 2L, 
                quote = FALSE)
  
  cat("\n\n")
  
  cat("Beta Coefficients:\n")
  
  print.default(format(x$estimates$beta$mean, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  
  cat("\nDelta activeness:\n")
  print.default(format(x$estimates$delta, digits = digits), print.gap = 2L, 
                quote = FALSE)
  
  cat("\nQ matrix:\n")
  print.default(format(x$estimates$q, digits = digits), print.gap = 2L, 
                quote = FALSE)
  
  cat("\n")
  invisible(x)
}


#' Sparse Latent Class Model for Cognitive Diagnosis (SLCM)
#'
#' Performs the Gibbs sampling routine for a sparse latent class model
#' as described in Chen et al. (2020) <doi: 10.1007/s11336-019-09693-2>
#'
#' @param y                 Item Matrix
#' @param k                 Dimension to estimate for Q matrix
#' @param burnin            Amount of Draws to Burn
#' @param chain_length      Number of Iterations for chain.
#' @param psi_invj,m0,bq    Additional tuning parameters.
#' 
#' @return
#'
#' An `slcm` object containing three named lists:
#'
#' - **`estimates`**
#'   - `beta`: Average beta coefficients
#'   - `theta`: Average theta coefficients
#'   - `delta`: Average activeness of coefficients
#'   - `class`: Average class membership
#'   - `pi`: Average attribute class probability.
#'   - `omega`: Average omega
#'   - `q`: Average activeness of Q matrix entries based on heuristic transformation.
#'   - `m2ll`: Average negative two times log-likelihood
#' - **`chain`**
#'   - `theta`: theta coefficients iterations
#'   - `beta`:  beta coefficients iterations
#'   - `class`:  class membership iterations
#'   - `pi`:  attribute class probability iterations
#'   - `omega`:  omega iterations
#'   - `m2ll`: Negative two times log-likelihood iterations
#' - **`details`**
#'   - `n`: Number of Subjects
#'   - `j`: Number of Items
#'   - `k`: Number of Traits
#'   - `l1`: Slab parameter
#'   - `m0`, `bq`: Additional tuning parameters
#'   - `burnin`: Number of Iterations to discard
#'   - `chain_length`: Number of Iterations to keep
#'   - `runtime`: Duration of model run inside of the C++ code. (Does not include summarization of MCMC chain.)
#'   - `package_version`: Version of the package the SLCM model was fit with.
#'   - `date_time`: Date and Time the model was fit.
#'
#' @details
#'
#' The **`estimates`** list contains the mean information from the sampling
#' procedure. Meanwhile, the **`chain`** list contains full MCMC values. Lastly,
#' the **`details`** list provides information regarding the estimation call.
#'
#' @rdname slcm
#' @export
#'
#' @examples
#' # Load Namespace
#' naming = requireNamespace("edmdata", quietly = TRUE)
#' 
#' # Use a demo data set from the paper
#' data("items_matrix_reasoning", package = "edmdata")
#'   
#' burnin = 50        # Set for demonstration purposes, increase to at least 1,000 in practice.
#' chain_length = 100 # Set for demonstration purposes, increase to at least 10,000 in practice.  
#'   
#' model_reasoning = slcm(items_matrix_reasoning, k = 4, 
#'                        burnin = burnin, chain_length = chain_length)
#'                          
#' print(model_reasoning)
slcm = function(
  y,
  k,
  burnin = 1000L,
  chain_length = 10000L,
  psi_invj = c(1, rep(2, 2 ^ k - 1)),
  m0 = 0,
  bq = 1) {

    # Perform some quality checks
    stopifnot(is.matrix(y))
    stopifnot(chain_length > burnin)
    
    # Disable order for now
    order = k
    
    # Restrict to binary approach
    m = 2
    
    # Check the correct number of psi_invj parameters are present.
    stopifnot(length(psi_invj) == 2^k)
    
    # Obtain information about when the model was fit.
    datetime_model_execution = Sys.time()

    # Launch routine and time it
    timing = system.time({
        model_mcmc <- slcm_cpp(y, k, 
                               M = m, order = order,
                               psi_invj = psi_invj,
                               m0 = m0, bq = bq,
                               burnin = burnin, chain_length = chain_length)
    })

    # Package object
    new_slcm_model(
        model_mcmc$chain,
        estimates = model_mcmc$estimates,
        details = list(
            n = nrow(y),
            j = ncol(y),
            k = k,
            l1 = psi_invj,
            m0 = m0,
            bq = bq,
            burnin = burnin,
            chain_length = chain_length,
            runtime = timing[["elapsed"]],
            package_version = as.character(utils::packageVersion("slcm")),
            date_time = as.character(datetime_model_execution)
        )
    )

}
