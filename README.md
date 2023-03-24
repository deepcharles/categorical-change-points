# Change-point detection in categorical time series

A categorical signal $(y_{t})_t$ is a discrete-valued sequence, meaning that $y_{t}$ can take only $D$ different values.
Computationnally, categorical signals are one-hot encoded: $y_t$ is represented by a multivariate binary signal $ (z_{t})_{t} $ where $z_{t} \in \{0,1 \}^D $ and $z_{t,d} = \mathbb{1}(y_t \text{ is equal to the d-th symbol})$.

## Python

To install the Python package, run in the repository

```
cd python/
python -m pip install .
```
 

## R

To install the R package...
