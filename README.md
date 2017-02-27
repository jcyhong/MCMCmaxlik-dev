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

* **Perry de Valpine**  - [Website](https://nature.berkeley.edu/~pdevalpine/)
* **Johnny Hong**  - [Website](https://jcyhong.github.io/)
* **Sara Stoudt** - [Website](https://www.stat.berkeley.edu/~sstoudt/)

## License

## Acknowledgments

This project started as a class project for Chris Paciorek's Bayesian Statistics class at UC Berkeley and was inspired by Perry de Valpine's idea to combine Fisher's identity with finite element approximation to estimate the gradient, as well as apply a one-dimensional sampler to estimate the step size. Perry has helped us work through some of the details coding in NIMBLE. We would like to thank Chris for his numerous suggestions along the way as well.
