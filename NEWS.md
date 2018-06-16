# fuser 1.0.1

- Introduced option to scale sum-squared loss per subgroup by the subgroup sample size. This improves performance in situations with unbalanced groups.
- Fixed warning caused by deprecated function in the Matrix package.

# fuser 1.0.0

First full release of fuser package for subgroup regression with fusion 
penalties.

## Features

- L1 and L2 fusion penalized regression across heterogeneous subgroups.
- Ability to deal with high-dimensional datasets (tens of thousands of covariates).
- C implementation for additional speed gains.
