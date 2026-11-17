"""BNPClust — Natarajan Section 5.1 Simulation Study in Python.

Replicates the 6 configurations of the simulation study in example/R/SimulatedData_NatarajanSection51.R
using the Python bindings (bnpclust) with direct class construction and MCMC sampling.
"""

import argparse
import os
import sys
import time
from typing import Any, Dict, List

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# Ensure build/ directory is in python module search path
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "build"))
)
import bnpclust as bnp
from utils import (
    generate_mixture_data,
    retrieve_W,
    set_hyperparameters,
)
from utils_plot import (
    plot_acf_U,
    plot_cls_est,
    plot_post_distr,
    plot_post_sim_matrix,
    plot_stats,
    plot_trace_cls,
    plot_trace_U,
    variation_of_information,
)


def run_simulation(
    scenarios: List[Dict[str, Any]],
    models: List[str],
    burn_in: int,
    num_samples: int,
    output_base_dir: str = "plot/simulated_data",
) -> List[Dict[str, Any]]:
    """Run simulation study across configurations."""
    total_iters = burn_in + num_samples
    results_summary: List[Dict[str, Any]] = []

    for sc in scenarios:
        print("\n=======================================================")
        print(f"Running Scenario: sigma = {sc['sigma']:.2f}, d = {sc['d']}")
        print("=======================================================")

        # Generate synthetic mixture data
        all_data, ground_truth = generate_mixture_data(
            N=100, sigma=sc["sigma"], dim=sc["d"], random_state=44
        )
        dist_mat = cdist(all_data, all_data)

        # Retrieve symmetric 8-NN spatial adjacency matrix
        W = retrieve_W(dist_mat, neighbours=8)

        # Estimate empirical hyperparameters from k-means initialization
        hyperparams = set_hyperparameters(
            dist_mat,
            k_elbow=3,
            random_state=44,
        )
        init_clusters = hyperparams["initial_clusters"]

        for m in models:
            config_name = f"{m}_{sc['tag']}"
            folder = os.path.join(output_base_dir, config_name)
            os.makedirs(folder, exist_ok=True)

            print(f"\n--- Model: {config_name} -> Output folder: {folder} ---")

            # 1. Parameter specification
            proc_param = bnp.NGGP_params(1.0, 0.1, 1.0)
            ut_param = bnp.utils_params(burn_in, num_samples, dist_mat)
            lik_param = bnp.Natarajan_params(
                hyperparams["delta1"],
                hyperparams["alpha"],
                hyperparams["beta"],
                hyperparams["delta2"],
                hyperparams["gamma"],
                hyperparams["zeta"],
            )

            # 2. Data and Likelihood
            data_obj = bnp.Data(ut_param, init_clusters)
            lik_obj = bnp.Natarajan_likelihood(data_obj, lik_param, ut_param)
            u_samp = bnp.RWMH(proc_param, data_obj, True, 2.0, True)

            # 3. Process prior (NGGP or NGGPx with spatial module)
            if m == "NGGP":
                proc_obj = bnp.NGGP(data_obj, proc_param, u_samp)
            else:
                mod_spatial = bnp.SpatialModule(data_obj, W, spatial_coefficient=1.0)
                proc_obj = bnp.NGGPx(data_obj, proc_param, u_samp, [mod_spatial])

            # 4. MCMC Samplers
            sm = bnp.SplitMerge_LSS_SDDS(data_obj, ut_param, lik_obj, proc_obj, True)
            neal3 = bnp.Neal3(data_obj, lik_obj, proc_obj)

            allocations_out: List[np.ndarray] = []
            K_out: List[int] = []
            U_out: List[float] = []

            print(
                f"Sampling {total_iters:,} iterations (Burn-in: {burn_in:,}, Samples: {num_samples:,})..."
            )
            t0 = time.time()
            for i in range(1, total_iters + 1):
                proc_obj.update_params()
                sm.step()
                if i % 25 == 0:
                    neal3.step()

                allocations_out.append(data_obj.get_allocations().copy())
                K_out.append(data_obj.get_K())
                U_out.append(u_samp.get_U())

                if i % 1000 == 0 or i == total_iters:
                    elapsed = time.time() - t0
                    print(
                        f"  Iteration {i:6d}/{total_iters:6d} | "
                        f"K: {K_out[-1]:2d} | "
                        f"U: {U_out[-1]:.3f} | "
                        f"Speed: {i / max(elapsed, 1e-6):.1f} iter/s"
                    )

            elapsed_time = time.time() - t0
            print(
                f"MCMC completed in {elapsed_time:.2f} seconds "
                f"({total_iters / elapsed_time:.1f} iter/s)"
            )

            mcmc_result = {
                "allocations": allocations_out,
                "K": K_out,
                "U": U_out,
                "elapsed_time": elapsed_time,
                "BI": burn_in,
                "NI": num_samples,
            }

            print("Saving plots and statistics...")
            pe = plot_cls_est(mcmc_result, BI=burn_in, save=True, folder=folder)
            plot_stats(mcmc_result, ground_truth, BI=burn_in, save=True, folder=folder)
            plot_post_distr(mcmc_result, BI=burn_in, save=True, folder=folder)
            plot_trace_cls(mcmc_result, BI=burn_in, save=True, folder=folder)
            plot_post_sim_matrix(mcmc_result, BI=burn_in, save=True, folder=folder)
            plot_trace_U(mcmc_result, BI=burn_in, save=True, folder=folder)
            plot_acf_U(mcmc_result, BI=burn_in, save=True, folder=folder)

            ari_val = float(adjusted_rand_score(ground_truth, pe))
            nmi_val = float(normalized_mutual_info_score(ground_truth, pe))
            _, nvi_val = variation_of_information(ground_truth, pe)
            k_est = int(len(np.unique(pe)))

            results_summary.append(
                {
                    "config": config_name,
                    "scenario": sc["tag"],
                    "model": m,
                    "K_est": k_est,
                    "ARI": ari_val,
                    "NMI": nmi_val,
                    "NVI": nvi_val,
                    "time_sec": elapsed_time,
                }
            )

    return results_summary


def main():
    parser = argparse.ArgumentParser(
        description="BNPClust Natarajan Section 5.1 Simulation Study (Python)"
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run a quick test with 500 burn-in and 1500 samples on a single scenario.",
    )
    parser.add_argument(
        "--burn-in",
        type=int,
        default=None,
        help="Number of burn-in iterations (default: 5000, or 500 if --quick).",
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        default=None,
        help="Number of post-burn-in iterations (default: 15000, or 1500 if --quick).",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        default=None,
        help="Scenario tags to run (choices: '018_10', '025_10', '020_50').",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=["NGGP", "NGGPW"],
        help="Models to run (choices: 'NGGP', 'NGGPW').",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="plot/simulated_data",
        help="Base output directory for plots and summary.",
    )
    args = parser.parse_args()

    all_scenarios = [
        {"tag": "018_10", "sigma": 0.18, "d": 10},
        {"tag": "025_10", "sigma": 0.25, "d": 10},
        {"tag": "020_50", "sigma": 0.20, "d": 50},
    ]

    if args.quick:
        burn_in = args.burn_in if args.burn_in is not None else 500
        num_samples = args.num_samples if args.num_samples is not None else 1500
        selected_scenarios = (
            [all_scenarios[0]]
            if args.scenarios is None
            else [s for s in all_scenarios if s["tag"] in args.scenarios]
        )
    else:
        burn_in = args.burn_in if args.burn_in is not None else 5000
        num_samples = args.num_samples if args.num_samples is not None else 15000
        if args.scenarios is None:
            selected_scenarios = all_scenarios
        else:
            selected_scenarios = [
                s for s in all_scenarios if s["tag"] in args.scenarios
            ]

    print("=======================================================")
    print("BNPClust Python Simulation Study")
    print(f"Scenarios: {[s['tag'] for s in selected_scenarios]}")
    print(f"Models: {args.models}")
    print(f"Burn-in: {burn_in:,} | Samples: {num_samples:,}")
    print(f"Output directory: {args.output_dir}")
    print("=======================================================")

    results = run_simulation(
        scenarios=selected_scenarios,
        models=args.models,
        burn_in=burn_in,
        num_samples=num_samples,
        output_base_dir=args.output_dir,
    )

    print("\n=======================================================")
    print("FINAL COMPARISON TABLE")
    print("=======================================================")
    header = f"{'Config':<15} {'Scenario':<10} {'Model':<8} {'K_est':<6} {'ARI':<8} {'NMI':<8} {'NVI':<8} {'Time (s)':<10}"
    print(header)
    print("-" * len(header))
    for r in results:
        print(
            f"{r['config']:<15} {r['scenario']:<10} {r['model']:<8} "
            f"{r['K_est']:<6d} {r['ARI']:<8.4f} {r['NMI']:<8.4f} "
            f"{r['NVI']:<8.4f} {r['time_sec']:<10.2f}"
        )

    # Save summary CSV
    csv_file = os.path.join(args.output_dir, "summary_comparison_python.csv")
    with open(csv_file, "w") as f:
        f.write("config,scenario,model,K_est,ARI,NMI,NVI,time_sec\n")
        for r in results:
            f.write(
                f"{r['config']},{r['scenario']},{r['model']},{r['K_est']},"
                f"{r['ARI']:.4f},{r['NMI']:.4f},{r['NVI']:.4f},{r['time_sec']:.2f}\n"
            )
    print(f"\nResults saved to {csv_file}\n")


if __name__ == "__main__":
    main()
