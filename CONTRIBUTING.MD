# Contributing to BNPClust

I welcome contributions to this project! Whether it's bug fixes, new features, or documentation improvements, your help is appreciated.

## Getting Started

### 1. Fork the Repository

Click the **"Fork"** button at the top right of the [BNPClust repository](https://github.com/Filippo-Galli/BNPClust) to create your own copy.

### 2. Clone Your Fork

Clone your forked repository to your local machine:

```bash
git clone https://github.com/yourusername/BNPClust.git
cd BNPClust
```

### 3. Create a New Branch

Create a new branch for your feature or bug fix:

```bash
git checkout -b feature-or-bugfix-name
```

Use descriptive branch names, e.g., `add-split-merge-sampler`, `fix-spatial-module-bug`.

### 4. Make Your Changes

Implement your changes in the codebase. Ensure that your code follows the existing style and conventions.

**Code Style Guidelines:**
- Follow the existing C++ naming conventions (CamelCase for classes, snake_case for functions)
- Use consistent indentation (same as surrounding code)
- Write clear, descriptive variable and function names
- Comment complex logic, especially in MCMC algorithms

### 5. Test Your Changes

Run existing tests and add new tests if necessary:

```bash
# Compile and run your changes
# Add your test cases to verify functionality
```

Ensure your changes:
- Compile without warnings
- Don't break existing functionality
- Pass new and existing tests

### 6. Commit Your Changes

Commit with a clear, descriptive message:

```bash
git add .
git commit -m "Description of your changes"
```

**Commit Message Guidelines:**
- Start with a verb: "Add", "Fix", "Improve", "Refactor"
- Be specific about what changed: `git commit -m "Add SAMS sampler with caching support"`
- Reference issues if applicable: `git commit -m "Fix spatial module bug (closes #42)"`

### 7. Push to Your Fork

Push your changes to your forked repository:

```bash
git push origin feature-or-bugfix-name
```

### 8. Create a Pull Request

Go to the [original repository](https://github.com/Filippo-Galli/BNPClust) and create a pull request from your forked branch. Provide:

- A clear title summarizing your changes
- A detailed description of what you changed and why
- Reference to any related issues
- Examples or test cases if applicable

### 9. Address Feedback

Be responsive to feedback and requests for changes from the maintainers. This is a collaborative process, and iteration is normal.

---

## Design Philosophy

BNPClust is built on a **modular architecture** to support flexible composition of components. When contributing, keep these principles in mind:

### KISS Principle

Follow the **Keep It Simple, Stupid (KISS)** principle:
- Avoid unnecessary complexity
- Focus on clear, maintainable code
- Solve problems directly without over-engineering
- Write code that others can understand easily

### Modularity

The framework is organized around four main logical components:

1. **Processes**: Bayesian nonparametric priors (DP, NGGP, with modular extensions)
2. **Likelihoods**: Data observation models
3. **Samplers**: MCMC inference algorithms
4. **Data**: Cluster assignments and data management

When adding new features:
- **New likelihood model?** → Inherit from the `Likelihood` abstract class so it can be swapped with other likelihoods
- **New sampler?** → Follow the `Sampler` interface for seamless integration
- **New process?** → Ensure it respects the `Process` abstract interface and supports modular extensions

### Code Organization

- **`BNPClust/processes/`**: Place new stochastic process implementations here
- **`BNPClust/likelihoods/`**: Add new likelihood models here
- **`BNPClust/samplers/`**: Implement new MCMC samplers here
- **`BNPClust/utils/`**: Add utility functions and base classes for shared infrastructure
- **`example/R/`**: Add R wrapper functions and examples

### Documentation

Code documentation is critical:

- **Header comments**: Explain the purpose of classes and methods
- **Inline comments**: Document complex logic, especially in MCMC algorithms
- **Parameter documentation**: Clearly describe function parameters and their constraints
- **References**: Cite papers when implementing algorithms (e.g., Neal, 2000; Jain & Neal, 2004)

Example:

```cpp
/**
 * Split-Merge MCMC sampler based on Jain & Neal (2004)
 * 
 * This sampler proposes splitting one cluster into two or merging two clusters
 * into one, with acceptance probability computed via Metropolis-Hastings.
 * 
 * References:
 *   Jain, S., & Neal, R. M. (2004). A Split-Merge Markov Chain Monte Carlo 
 *   Procedure for the Dirichlet Process Mixture Model.
 *   Journal of Computational and Graphical Statistics, 13(1), 158-182.
 */
class SplitMergeSampler : public Sampler { ... };
```

### R Bindings

BNPClust uses **Rcpp external pointers (XPtr)** pattern for seamless C++ object management in R. When adding new C++ features, expose them via `BNPClust/r_bindings.cpp`:

**Factory Functions:**
Every C++ class needs a factory function to create instances from R:

```cpp
// [[Rcpp::export]]
Rcpp::XPtr<MyNewClass> create_MyNewClass(Rcpp::XPtr<Params> params, 
                                         SEXP data_sexp) {
    Data *data = get_data_ptr(data_sexp);  // Extract C++ pointer from R
    return Rcpp::XPtr<MyNewClass>(new MyNewClass(*data, *params), true);
}
```

**Wrapper Functions:**
Expose methods via simple wrapper functions:

```cpp
// [[Rcpp::export]]
void mynewclass_step(Rcpp::XPtr<MyNewClass> obj) { 
    obj->step(); 
}

// [[Rcpp::export]]
Eigen::VectorXi mynewclass_get_result(Rcpp::XPtr<MyNewClass> obj) { 
    return obj->get_result(); 
}
```

**Usage in R:**

```R
# Create objects using factory functions
params <- create_Params(delta1=1.0, alpha=1.0, ...)
data <- create_Data(params, initial_allocations)
sampler <- create_MyNewClass(params, data)

# Call methods via wrappers
mynewclass_step(sampler)
result <- mynewclass_get_result(sampler)
```

**Key Points:**
- Use `get_data_ptr()` helper to extract Data pointers (handles both `Data` and `Datax` types)
- Add `// [[Rcpp::export]]` above every function meant to be accessible from R
- Memory is managed by Rcpp (set second parameter to `true` for automatic cleanup)
- Add examples in `R/launcher.R` or `R/mcmc_loop.R` showing usage

---

## Types of Contributions

### 🐛 Bug Reports

- Open an issue with a clear title and description
- Include minimal reproducible example
- Specify your OS, R version, C++ compiler version

### ✨ Feature Requests

- Check if the feature already exists or is planned
- Describe the use case and why it would be valuable
- Discuss implementation approach if possible

### 📖 Documentation

- Fix typos, clarify explanations, add examples
- Improve Doxygen documentation
- Update README or contributing guidelines
- Add academic references or method citations

### 🔧 Code Improvements

- Refactor for clarity and efficiency
- Optimize MCMC sampling performance
- Add missing error handling
- Improve test coverage

---

## Questions?

If you have questions about contributing, feel free to:
- Open a discussion in the GitHub repository
- Check the documentation at [https://filippo-galli.github.io/BNPClust/](https://filippo-galli.github.io/BNPClust/)

Thank you for contributing to BNPClust! 🙏
