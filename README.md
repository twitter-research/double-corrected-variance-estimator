# Double Corrected Variance Estimator

[![status: static](https://opensource.twitter.dev/status/static.svg)](https://opensource.twitter.dev/status/#static)

This is a repo for the code used for reproducing our [De-biasing "bias" measurement paper](https://doi.org/10.1145/3531146.3533105).

If you plan to use this code please cite our paper as follows:

```
@inproceedings{
	address = {Seoul, Republic of Korea},
	title = {De-biasing "bias" measurement},
	doi = {10.1145/3531146.3533105},
	booktitle = {2022 {ACM} {Conference} on {Fairness}, {Accountability}, and {Transparency} ({FAccT} ’22)},
	publisher = {ACM},
	author = {Lum, Kristian and Zhang, Yunfeng and Bower, Amanda},
	month = jun,
	year = {2022},
}
```

# Instructions
- clone the repo.
- Run `RScript bootstrap_variance_estimates_fast.R` to reproduce the results in the paper.
- Run `adult income.ipynb` to reproduce the results in Section 6.
