# Obtain summary statistics on each component of the MCMC chain
summarize_mcmc = function(x, margin, take_quantiles = TRUE, transpose_quantiles = FALSE) {
    if (take_quantiles) {
        quantiles = apply(
            x,
            MARGIN = margin,
            FUN = stats::quantile,
            probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
            na.rm = TRUE
        )

        if (transpose_quantiles) {
            quantiles = t(quantiles)
        }
    } else {
        quantiles = NULL
    }

    list(
        mean      = apply(x, MARGIN = margin, FUN = mean,
                          na.rm = TRUE),
        sd        = apply(x, MARGIN = margin, FUN = stats::sd,
                          na.rm = TRUE),
        quantiles = quantiles
    )
}

# Obtain an output related to the dimension of how the MCMC chain was stored.
summarize_1d_vector_output = function(x) {
    summarize_mcmc(x, margin = 2, transpose_quantiles = TRUE)
}

summarize_2d_matrix_output = function(x) {
    summarize_mcmc(x, margin = 1, transpose_quantiles = TRUE)
}

summarize_3d_array_output = function(x) {
    summarize_mcmc(x, margin = c(1, 2))
}

summarize_4d_array_output = function(x) {
    summarize_mcmc(x, margin = c(1, 2, 3), take_quantiles = FALSE)
}

# Route to the correct summary manipulation
summarize_model = function(x) {
    # Retrieve properties of object
    is_class_matrix = is.matrix(x)
    is_class_array = is.array(x)
    type_x = typeof(x)
    dim_x = dim(x)
    n_cols = dim_x[2]
    n_dim_x = length(dim(x))

    if (is_class_matrix && type_x == "double" && n_cols == 1L) {
        # Handles the 1D matrix case
        summarize_1d_vector_output(x)
    } else if (is_class_matrix && type_x == "double" && n_dim_x == 2) {
        # Handles the 2D matrix case
        summarize_2d_matrix_output(x)
    } else if(is_class_array && type_x == "double" && n_dim_x == 3) {
        # Handles the 3D cube case
        summarize_3d_array_output(x)
    } else if(is_class_array && type_x == "double" && n_dim_x == 4) {
        # Handles the 4D cube case
        summarize_4d_array_output(x)
    } else {
        stop("Unable to summarize the element model.")
    }
}
