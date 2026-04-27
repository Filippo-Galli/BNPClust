# BNPClust: A C++ Framework for Bayesian Nonparametric Clustering with MCMC

Inspired by the need for flexible and efficient tools in Bayesian nonparametric clustering, this repository presents a C++ framework that leverages Markov Chain Monte Carlo (MCMC) methods. The framework is designed to be modular and extensible, allowing users to implement various stochastic processes, likelihood models, and sampling strategies.

**Documentation**: [https://filippo-galli.github.io/BNPClust/](https://filippo-galli.github.io/BNPClust/)

**License**: GPL-3.0 — See [LICENSE](LICENSE) file for details

## 🚀 Quick Start

To get started quickly, clone the repository and open it in R/RStudio. Ensure you have the required R packages installed (`Rcpp`, `RcppEigen`).

```bash
git clone https://github.com/filippo-galli/BNPClust.git
cd BNPClust
```

Open `R/mcmc_loop.R` to create your own preferred MCMC scheme and see how to create and pass data to each component.

An example of a complete MCMC scheme is provided in `R/launcher.R`, which uses `R/mcmc_loop.R` as the underlying MCMC loop.

### 📊 Data Retrieval

Check the [guide](R/data_retrieval/README.md)

## 📁 Directory Structure

- `R/`: R interface functions and scripts for MCMC orchestration and analysis
- `src/`: C++ source code implementing the core framework
- `docs/`: Generated documentation (Doxygen)
- `doxygen-theme/`: Doxygen theme configuration files

### R Scripts

Each C++ class in the framework has corresponding R bindings to facilitate interaction with R users. The typical workflow is:

- **`R/launcher.R`**: Entry point to load data, set parameters, invoke MCMC functions, and save results
- **`R/mcmc_loop.R`**: The MCMC iteration loop where you select the combination of Process, Likelihood, and Sampler to use
- **`R/mcmc_analysis.R`**: Functions to visualize and summarize MCMC output
- Other scripts provide utility functions for data loading, visualization, result plotting, and data fetching/cleaning

### src: C++ Core Framework

The C++ source code is organized into the following subdirectories:

- **`src/processes/`**: Implementations of stochastic processes (e.g., Dirichlet Process, Normalized Generalized Gamma Process) and their modular extensions (continuous covariate, spatial, and binary covariate modules)
- **`src/likelihoods/`**: Likelihood model implementations, including distance-based clustering, gamma likelihood, and null likelihood
- **`src/samplers/`**: MCMC sampling algorithms (Neal's Algorithm 3, Split-Merge variants, SAMS, etc.)
- **`src/utils/`**: Utility functions, base classes, and shared infrastructure
- **`src/bindings.cpp`**: Rcpp bindings exposing C++ classes and functions to R

### docs and Documentation

The `docs/` folder contains documentation generated with Doxygen, including class diagrams and detailed code documentation. To regenerate the documentation, ensure Doxygen is installed and run:

```bash
doxygen Doxyfile
```

## 🏗️ Architecture

The framework is built around five main logical components:

1. **`Params`**: Manages model hyperparameters and MCMC configuration
2. **`Likelihood`**: Defines the data observation model
3. **`Process`**: Handles the Bayesian nonparametric prior
4. **`Data`**: Manages cluster assignments and data handling
5. **`Sampler`**: Implements the MCMC inference engine

## 📦 Implemented Methods

### Stochastic Processes

Currently implemented processes:

- **Dirichlet Process (DP)**: The foundational nonparametric prior
- **Normalized Generalized Gamma Process (NGGP)**: A flexible family of which the DP is a special case

The framework supports modular extensions to incorporate domain-specific structure:

- **Continuous Covariate Module**: Incorporates continuous covariates into the clustering process
- **Spatial Module**: Accounts for spatial dependencies given a neighbor adjacency matrix
- **Binary Covariate Module**: Incorporates discrete covariates into the clustering process

Cached versions of these modules are available for improved computational performance.

### Likelihood Models

All likelihood components are implemented as log-likelihood functions for numerical stability. Currently available:

- **Distance-based Clustering**: Based on "Cohesion and Repulsion in Bayesian Distance Clustering" (Natarajan et al., 2023)
- **Gamma Likelihood**: Variant of distance-based clustering without the repulsion term
- **Null Likelihood**: A placeholder that does not contribute to the posterior, useful for prior inspection

### MCMC Samplers

The framework includes several MCMC sampling strategies suited to different problem structures:

- **Neal's Algorithm 3**: Conjugate Gibbs sampler for efficient cluster updates
- **ZDNAM**: Gibbs sampling with Zero-self Downward Nested Antithetic Modification (Neal, 2024) — *experimental*
- **Split-Merge (Jain & Neal, 2004)**: Standard Split-Merge MCMC procedure
- **SAMS**: Sequentially-Allocated Merge-Split sampler (Dahl, 2021)
- **LSS Split-Merge**: Locality Sensitive Sampling for scalable Split-Merge (Luo et al., 2018)
- **LSS-SDDS Split-Merge**: Split-Merge with Smart-Dumb/Dumb-Smart moves, LSS, and SAMS enhancements

**Status Legend**: Core methods (Neal's Algorithm 3, standard Split-Merge) are production-ready; samplers marked experimental are under development.

## 🛠️ Installation & Usage

### Prerequisites

- **R** (version 3.5.0 or later)
- **C++ Compiler** supporting C++11 or later (e.g., g++, clang)
- **R Packages**: `Rcpp`, `RcppEigen`

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/filippo-galli/BNPClust.git
   cd BNPClust
   ```

2. Open the project in R/RStudio

3. Install dependencies:
   ```R
   install.packages(c("Rcpp", "RcppEigen"))
   ```

### Running a Basic Example

The Rcpp bindings in `src/bindings.cpp` expose all C++ classes and functions to R, allowing you to compose custom MCMC schemes directly from R. Use `R/launcher.R` as a template:

1. Configure your Process, Likelihood, and Sampler in `R/mcmc_loop.R`
2. Prepare your data and hyperparameters in `R/launcher.R`
3. Execute `R/launcher.R` to run the MCMC chain
4. Analyze results using `R/mcmc_analysis.R`

## 📚 References

- Natarajan, L., et al. (2023). *Cohesion and Repulsion in Bayesian Distance Clustering*
- Neal, R. M. (2000). *Markov Chain Sampling Methods for Dirichlet Process Mixture Models*. Journal of Computational and Graphical Statistics, 9(2), 249–265
- Neal, R. M. (2024). *Modifying Gibbs Sampling to Avoid Self Transitions*
- Jain, S., & Neal, R. M. (2004). *A Split-Merge Markov Chain Monte Carlo Procedure for the Dirichlet Process Mixture Model*. Journal of Computational and Graphical Statistics, 13(1), 158–182
- Dahl, D. B. (2021). *Sequentially-allocated merge-split samplers for conjugate Bayesian nonparametric models*
- Luo, L., et al. (2018). *Scaling-up Split-Merge MCMC with Locality Sensitive Sampling (LSS)*

## 💡 Usage & Licensing

BNPClust is released under the **GPL-3.0 license**, which means:

- ✅ **Use freely**: You can use BNPClust for any purpose, including commercial applications
- ✅ **Modify and extend**: You can modify the code to suit your needs
- ✅ **Study the source**: Full source code access for learning and verification
- 📋 **Share improvements**: If you distribute modified versions, they must also be released under GPL-3.0 with source code provided

For details on GPL-3.0 compliance and your obligations, see the [LICENSE](LICENSE) file.

---

## 💡 Contributing

See [CONTRIBUTING.md](CONTRIBUTING.MD) for guidelines on how to contribute to BNPClust.


If you find this project useful, please leave a star! ⭐
