# Sampling-based Maximum Likelihood Estimation for Latent Variable Models

This project provides efficient sampling-based methods for maximum likelihood estimation in the context of latent variable modeling. The key algorithms revolve around stochastic gradient methods.

## Getting Started

### Prerequisites

* [R](https://www.r-project.org/)
* R packages: `devtools`, `nimble` with automatic differentiaion

The folder nimbleAD is the package `nimble` with automatic differentiaion.

### Installing

Set the working directory to the folder containing MCMCmaxlik. Then use `install()` in `devtools` to install the package.

```
library("devtools")
install("MCMCmaxlik")
```

## Built With

* [R](https://www.r-project.org/) - the software environment
* [nimble](https://r-nimble.org/) - NIMBLE (NIMBLE: Numerical Inference for Hierarchical Models Using Bayesian and Likelihood Estimation)

## Authors

* **Johnny Hong**  - [Website](https://jcyhong.github.io/)
* **Sara Stoudt** - [Website](https://www.stat.berkeley.edu/~sstoudt/)
* **Perry de Valpine**  - [Website](https://nature.berkeley.edu/~pdevalpine/)
