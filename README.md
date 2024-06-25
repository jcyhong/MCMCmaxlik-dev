# Fast maximum likelihood estimation for general hierarchical models

This project provides efficient sampling-based methods for maximum likelihood estimation in the context of general hierarchical models. The key algorithms revolve around stochastic gradient methods.

## Getting Started

### Prerequisites

* [R](https://www.r-project.org/)
* R packages: `devtools`, `nimble` with automatic differentiation (install the package from the branch [ADoak](https://github.com/nimble-dev/nimble/tree/ADoak))

### Installing

Set the working directory to the folder containing MCMCmaxlik. Then use `install()` in `devtools` to install the package.

```
library("devtools")
install("MCMCmaxlik")
```

## Built With

* [R](https://www.r-project.org/) - the software environment
* [nimble](https://r-nimble.org/) - NIMBLE (NIMBLE: Numerical Inference for Hierarchical Models Using Bayesian and Likelihood Estimation)
* Requires ADoak branch of NIMBLE

## Authors

* **Johnny Hong**  - [Website](https://jcyhong.github.io/)
* **Sara Stoudt** - [Website](https://sastoudt.github.io/)
* **Perry de Valpine**  - [Website](https://nature.berkeley.edu/~pdevalpine/)
