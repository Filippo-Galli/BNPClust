library(Rcpp)
library(RcppEigen)

# Load the C++ module (linking prebuilt core library if available)
core_lib_dir <- normalizePath("build", mustWork = FALSE)

if (.Platform$OS.type == "windows") {
    core_lib_path <- file.path(core_lib_dir, "bnpclust_core.dll")
} else if (Sys.info()[["sysname"]] == "Darwin") {
    core_lib_path <- file.path(core_lib_dir, "libbnpclust_core.dylib")
} else {
    core_lib_path <- file.path(core_lib_dir, "libbnpclust_core.so")
}

if (!file.exists(core_lib_path)) {
    stop(
        "Core library not found. Build it with: cmake -S . -B build && cmake --build build"
    )
}

if (.Platform$OS.type == "windows") {
    Sys.setenv(PKG_LIBS = paste0("-L", core_lib_dir, " -lbnpclust_core"))
    Sys.setenv(
        PATH = paste(core_lib_dir, Sys.getenv("PATH"), sep = .Platform$path.sep)
    )
} else {
    Sys.setenv(
        PKG_LIBS = paste0(
            "-L",
            core_lib_dir,
            " -lbnpclust_core -Wl,-rpath,",
            core_lib_dir
        )
    )
}

rcpp_eigen_include <- system.file("include", package = "RcppEigen")
if (!nzchar(rcpp_eigen_include)) {
    stop("RcppEigen include directory not found. Is RcppEigen installed?")
}
existing_cppflags <- Sys.getenv("PKG_CPPFLAGS")
eigen_cppflags <- paste0("-I", rcpp_eigen_include)
Sys.setenv(PKG_CPPFLAGS = paste(existing_cppflags, eigen_cppflags))

object_files <- list.files(
    "src",
    pattern = "\\.o$",
    recursive = TRUE,
    full.names = TRUE
)
if (length(object_files) > 0) {
    unlink(object_files)
}

Rcpp::sourceCpp("src/r_bindings.cpp", cacheDir = "tmp/rcpp_cache")

run_mcmc <- function(
    likelihood_param,
    process_param,
    utils_param,
    initial_allocations = integer(0),
    W,
    continuos_covariates = NULL,
    binary_covariates = NULL,
    categorical_covariates = NULL
) {
    # Ensure types are correct for C++
    initial_allocations <- as.integer(initial_allocations)

    continuos_cache <- create_Continuos_cache(
        initial_allocations,
        continuos_covariates
    )
    binary_cache <- create_Binary_cache(initial_allocations, binary_covariates)
    # spatial_cache <- create_Spatial_cache(initial_allocations, W)

    # Instantiate Data using factory function
    # data <- create_Datax(
    #     utils_param,
    #     list(binary_cache, continuos_cache),
    #     initial_allocations
    # )
    data <- create_Data(utils_param, initial_allocations)

    # Instantiate Likelihood using factory function
    likelihood <- create_GaussianMixtureModel_likelihood(data, likelihood_param)
    # likelihood <- create_Natarajan_likelihood(data, likelihood_param, utils_param)
    # likelihood <- create_Natarajan_likelihood_summaryStats(data, params)
    # likelihood <- create_Null_likelihood(data, params) # Placeholder likelihood
    # likelihood <- create_Gamma_likelihood(data, params)

    # Instantiate U_sampler (RWMH) using factory function
    # Constructor: Params&, Data&, bool use_V, double proposal_sd, bool tuning_enabled
    # u_sampler <- create_RWMH(process_param, data, TRUE, 2.0, TRUE)
    u_sampler <- NULL

    # Instantiate Process (NGGPx) using modules
    # 1. Spatial module
    #mod_spatial <- create_SpatialModule(data, W, spatial_coefficient = 0.1)
    # mod_spatial <- create_SpatialModuleContinuous(data, W, spatial_coefficient = 1.0)
    # mod_spatial <- create_SpatialModuleCache(data, cache = spatial_cache, spatial_coefficient = 1.0)

    # 2. Covariate module (cached)
    fixed_v <- TRUE
    B <- 10 * var(continuos_covariates) # prior variance
    m <- 0 # prior mean
    v <- 0.5 * var(continuos_covariates) # known variance
    nu <- 1
    S0 <- 1.0

    mod_cont <- create_ContinuosCovariatesModuleCache(
        data,
        continuos_cache,
        fixed_v,
        m,
        B,
        v,
        nu,
        S0
    )
    # mod_cont <- create_ContinuosCovariatesModule(data, continuos_covariates, fixed_v = TRUE, m = m, B = B, v = v)

    # 3. Binary covariate module
    # mod_binary <- create_BinaryCovariatesModule(data, binary_covariates, 0.1, 0.1)
    mod_binary <- create_BinaryCovariatesModuleCache(
        data,
        binary_cache,
        0.1,
        0.1
    )

    # 4. Categorical covariate module
    # alphas <- rep(1.0, length(unique(categorical_covariates)))
    # mod_categorical <- create_CategoricalCovariatesModule(data, categorical_covariates, alphas)

    # Combine modules into NGGPx process
    # process <- create_NGGPx(data, process_param, u_sampler, list(mod_spatial, mod_cont, mod_binary))
    # process <- create_NGGP(data, params, u_sampler)
    process <- create_DP(
        data,
        process_param
    )

    # Instantiate Sampler (SplitMerge_LSS_SDDS) using factory function
    sm <- create_SplitMerge(
        data,
        likelihood,
        process,
        shuffle = TRUE
    )

    neal3 <- create_Neal3(data, likelihood, process)

    # Get parameters for loop using getter functions
    BI <- params_get_BI(utils_param)
    NI <- params_get_NI(utils_param)
    total_iters <- BI + NI

    # Results storage
    allocations_out <- vector("list", total_iters)
    K_out <- integer(total_iters)
    U_out <- rep(NA_real_, total_iters)

    cat("Starting MCMC with", NI, "iterations after", BI, "burn-in...\n")

    start_time <- Sys.time()

    for (i in 1:total_iters) {
        # Update process parameters (U)
        process_update_params(process)

        # MCMC Step
        sampler_step(sm)

        # Neal3 Step
        if (i %% 25 == 0) {
            sampler_step(neal3)
        }
        # Store results
        allocations_out[[i]] <- data_get_allocations(data)
        K_out[i] <- data_get_K(data)
        if (!is.null(u_sampler)) {
            U_out[i] <- u_sampler_get_U(u_sampler)
        }

        # Progress
        if (i %% max(1, floor(total_iters / 20)) == 0) {
            elapsed <- as.numeric(difftime(
                Sys.time(),
                start_time,
                units = "secs"
            ))
            iter_per_sec <- i / elapsed
            eta <- (total_iters - i) / iter_per_sec
            cat(sprintf(
                "Iteration %d: Clusters: %d - iter/s: %.2f eta: %.2f\n ",
                i,
                data_get_K(data),
                iter_per_sec,
                eta
            ))
        }
    }

    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat("MCMC completed.\n")
    cat("Total time (secs):", elapsed_time, "\n")
    if (!is.null(u_sampler)) {
        cat(
            "U acceptance rate:",
            u_sampler_get_acceptance_rate(u_sampler) * 100,
            "%\n"
        )
    }
    if (exists("sampler")) {
        lss_sdds_accepted_moves(sampler)
    }

    return(list(
        allocations = allocations_out,
        K = K_out,
        U = U_out,
        BI = BI,
        NI = NI,
        elapsed_time = elapsed_time
    ))
}
