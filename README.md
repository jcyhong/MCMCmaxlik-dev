# Sampling-based Maximum Likelihood Estimation for Latent Variable Models

This project provides efficient sampling-based methods for maximum likelihood estimation in the context of latent variable modeling. The key algorithms revolve around stochastic gradient methods.

## Getting Started

### Prerequisites

* [R](https://www.r-project.org/)
* R packages: `devtools`, `nimble`


### Installing

First, set the working directory to the folder containing MCMCmaxlik. Then use `install()` in `devtools` to install the package.

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

## Acknowledgments

We thank Chris Paciorek for his numerous suggestions in the project, in particular the idea of warm start in the MCMC sampling. We thank Nick Michaud for NIMBLE's implementation of MCEM. The second author is supported by a  National Physical Sciences Consortium fellowship. 
