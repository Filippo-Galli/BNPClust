# BNPClust: A C++ Framework for Bayesian Nonparametric Clustering with MCMC

Inspired by the need for flexible and efficient tools in Bayesian nonparametric clustering, this repository presents a C++ framework that leverages Markov Chain Monte Carlo (MCMC) methods. The framework is designed to be modular and extensible, allowing users to implement various stochastic processes, likelihood models, and sampling strategies.

**Documentation**: [https://filippo-galli.github.io/BNPClust/](https://filippo-galli.github.io/BNPClust/)

**License**: GPL-3.0 — See [LICENSE](LICENSE) file for details

## Quick Start

To get started quickly, clone the repository and open it in your preferred R/Python (bindings in progress) environment.

```bash
git clone https://github.com/filippo-galli/BNPClust.git
cd BNPClust
```

The core library is implemented in C++ so you need to compile it before using it in R/Python. You can do this by running the following command in the terminal:

```bash
cmake -S . -B build
cmake --build build --parallel
```

This will compile the C++ library and place the executable in `build/`. Once done, you can use the library in your R/Python environment.

As example of how it works in R, you can take a look at the [Natarajan example](R/example/Natarajan_NGGP.R) script.

### API Documentation
For R users, I suggest to take a look at the [R interface](src/r_bindings.cpp). For Python users, the API is documented in the [Python interface](src/python_bindings.cpp). For the C++ interface the docs is available at [Doxygen documentation](https://filippo-galli.github.io/BNPClust/).

## 📁 Directory Structure

- `R/`: R interface functions and scripts for MCMC orchestration and analysis
- `src/`: C++ source code implementing the core framework
- `docs/`: Generated documentation (Doxygen)
- `doxygen-theme/`: Doxygen theme configuration files

### src: C++ Core Framework

The C++ source code is organized into the following subdirectories:

- **`src/processes/`**: Implementations of stochastic processes (e.g., Dirichlet Process, Normalized Generalized Gamma Process) and their modular covariates extensions.
- **`src/likelihoods/`**: Likelihood model implementations, including distance-based clustering, gamma likelihood, null likelihood and a Mixture of Gaussian Likelihood.
- **`src/samplers/`**: MCMC sampling algorithms (Neal's Algorithm 3, Split-Merge variants, SAMS, etc.)
- **`src/utils/`**: Utility functions, base classes, and shared infrastructure
- **`src/r_bindings.cpp`**: Rcpp bindings exposing C++ classes and functions to R

### docs and Documentation

The `docs/` folder contains documentation generated with Doxygen, including class diagrams and detailed code documentation. To regenerate the documentation, ensure Doxygen is installed and run:

```bash
doxygen Doxyfile
```

## 🏗️ Architecture

The framework is built around five main logical components:

1. **`likelihood/`**: Defines the data observation model
2. **`process/`**: Handles the Bayesian nonparametric prior
3. **`params/`**: Manages process, likelihood and sampler hyperparameters
4. **`sampler/`**: Implements the MCMC inference engine
5. **`utils/`**: Utility objects for data handling

Each part of the framework is developed to work inside the library but also give the opportunity to extract a specific class and use it as a standalone component.

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
