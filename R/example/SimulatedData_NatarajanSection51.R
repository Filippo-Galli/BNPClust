##############################################################################
# BNPClust — Natarajan Section 5.1 Simulation Study
#
# Re-runs all 6 configurations of the simulation study under the updated
# MCMC algorithm and updates all results and plots in plot/simulated_data/
##############################################################################

source("R/utils.R")
source("R/utils_plot.R")
library(Rcpp)
library(coda)
library(salso)
library(aricode)

dyn.load("build/libbnpclust_r.so")
bnp_mod <- Rcpp::Module("bnpclust_module", "libbnpclust_r")

scenarios <- list(
    list(tag = "018_10", sigma = 0.18, d = 10),
    list(tag = "025_10", sigma = 0.25, d = 10),
    list(tag = "020_50", sigma = 0.20, d = 50)
)

models <- c("NGGP", "NGGPW")

results_summary <- list()

for (sc in scenarios) {
    cat("\n=======================================================\n")
    cat(sprintf("Running Scenario: sigma = %.2f, d = %d\n", sc$sigma, sc$d))
    cat("=======================================================\n")

    set.seed(44)
    data_gen <- generate_mixture_data(N = 100, sigma = sc$sigma, dim = sc$d)
    all_data <- data_gen$points
    ground_truth <- data_gen$clusts
    dist_mat <- as.matrix(dist(all_data))

    W <- retrieve_W(dist_mat, neighbours = 8)
    W <- matrix(as.integer(W), nrow = nrow(W), ncol = ncol(W))

    hyperparams <- set_hyperparameters(
        dist_mat,
        k_elbow = 3,
        plot_clustering = FALSE,
        plot_distribution = FALSE
    )
    init_clusters <- as.integer(hyperparams$initial_clusters - 1)

    for (m in models) {
        config_name <- paste0(m, "_", sc$tag)
        folder <- paste0("plot/simulated_data/", config_name, "/")
        dir.create(folder, recursive = TRUE, showWarnings = FALSE)

        cat(sprintf("\n--- Model: %s -> Output folder: %s ---\n", config_name, folder))

        BI <- 5000
        NI <- 15000
        total_iters <- BI + NI

        proc_param <- bnp_mod$create_NGGP_params(1.0, 0.1, 1.0)
        ut_param <- bnp_mod$create_utils_params(BI, NI, dist_mat)
        lik_param <- bnp_mod$create_Natarajan_params(
            hyperparams$delta1,
            hyperparams$alpha,
            hyperparams$beta,
            hyperparams$delta2,
            hyperparams$gamma,
            hyperparams$zeta
        )

        data_obj <- bnp_mod$create_Data(ut_param, init_clusters)
        lik_obj <- bnp_mod$create_Natarajan_likelihood(data_obj, lik_param, ut_param)
        u_samp <- bnp_mod$create_RWMH(proc_param, data_obj, TRUE, 2.0, TRUE)

        if (m == "NGGP") {
            proc_obj <- bnp_mod$create_NGGP(data_obj, proc_param, u_samp)
        } else {
            mod_spatial <- bnp_mod$create_SpatialModule(data_obj, W, spatial_coefficient = 1.0)
            proc_obj <- bnp_mod$create_NGGPx(data_obj, proc_param, u_samp, list(mod_spatial))
        }

        sm <- bnp_mod$create_SplitMerge_LSS_SDDS(data_obj, ut_param, lik_obj, proc_obj, TRUE)
        neal3 <- bnp_mod$create_Neal3(data_obj, lik_obj, proc_obj)

        allocations_out <- vector("list", total_iters)
        K_out <- integer(total_iters)
        U_out <- rep(NA_real_, total_iters)

        cat("Sampling 20,000 iterations...\n")
        t0 <- Sys.time()
        for (i in 1:total_iters) {
            bnp_mod$process_update_params(proc_obj)
            bnp_mod$sampler_step(sm)
            if (i %% 25 == 0) {
                bnp_mod$sampler_step(neal3)
            }
            allocations_out[[i]] <- bnp_mod$data_get_allocations(data_obj)
            K_out[i] <- bnp_mod$data_get_K(data_obj)
            U_out[i] <- bnp_mod$u_sampler_get_U(u_samp)
        }
        t1 <- Sys.time()
        elapsed_time <- as.numeric(difftime(t1, t0, units = "secs"))
        cat(sprintf("MCMC completed in %.2f seconds (%.1f iter/s)\n", elapsed_time, total_iters / elapsed_time))

        mcmc_result <- list(
            allocations = allocations_out,
            K = K_out,
            U = U_out,
            elapsed_time = elapsed_time,
            BI = BI,
            NI = NI
        )

        cat("Saving plots and stats...\n")
        pe <- plot_cls_est(mcmc_result, BI = BI, save = TRUE, folder = folder)
        plot_stats(mcmc_result, ground_truth, BI = BI, save = TRUE, folder = folder)
        plot_post_distr(mcmc_result, BI = BI, save = TRUE, folder = folder)
        plot_trace_cls(mcmc_result, BI = BI, save = TRUE, folder = folder)
        plot_post_sim_matrix(mcmc_result, BI = BI, save = TRUE, folder = folder)
        plot_trace_U(mcmc_result, BI = BI, save = TRUE, folder = folder)
        plot_acf_U(mcmc_result, BI = BI, save = TRUE, folder = folder)

        pe_vec <- as.vector(pe)
        ari_val <- arandi(pe, ground_truth)
        nmi_val <- NMI(pe_vec, ground_truth)
        vi_val <- NVI(pe_vec, ground_truth)
        k_est <- length(unique(pe))

        results_summary[[config_name]] <- list(
            config = config_name,
            scenario = sc$tag,
            model = m,
            K_est = k_est,
            ARI = ari_val,
            NMI = nmi_val,
            NVI = vi_val,
            time = elapsed_time
        )
    }
}

cat("\n=======================================================\n")
cat("FINAL COMPARISON TABLE ACROSS ALL 6 RUNS\n")
cat("=======================================================\n")
df_res <- do.call(rbind, lapply(results_summary, as.data.frame))
print(df_res)
write.csv(df_res, file = "plot/simulated_data/summary_comparison.csv", row.names = FALSE)
cat("\nResults saved to plot/simulated_data/summary_comparison.csv\n")
