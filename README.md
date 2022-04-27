# Double Corrected Variance Estimator

[![status: static](https://opensource.twitter.dev/status/static.svg)](https://opensource.twitter.dev/status/#static)

This is a repo for the code used for reproducing our [De-biasing "bias" measurement paper]() (will add link when it's out). 

If you plan to use this code please cite our paper as follows:

```
@ARTICLE{debiasing-bias-measurement,
       author = {Lum, Kristian and Zhang, Yunfeng and Bower, Amanda},
        title = {De-biasing "bias" measurement},
         year = 2022
}
```

# Instructions
- clone the repo.
- Run `RScript bootstrap_variance_estimates_fast.R` to reproduce the results in the paper.
- Run `adult income.ipynb` to reproduce the results in Section 6.
