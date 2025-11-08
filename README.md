# Square root QR covariance steering

## TODO
- Implement comparison with other methods (e.g., standard covariance steering, LQG)

# Observations from experiments
- When the QR method converges, it always matches the optimal cost of the full covariance steering method.
- Sometimes, the QR method stalls even when the full covariance method is solved. --- Why?
- Updating the SCvx* algorithm penalty might improve convergence rates
- Comparison (solution time, optimality, condition number) with varying state sizes
- Comparison with varying horizon lengths

- Come up with a bit more complex chance constrained problem where the full covariance might converge to local optimum.
- Convexity in different coordinate systems. Is there a paper out there?
- SCvxStar with inexact linearization? proof? 
Random chance constraints, discard ones where constraints are not active
Ill-scaled problem (half of the state has small magnitude)

## Debugging numerical issues in MOSEK

In order to debug numerical issues in MOSEK, follow these steps:
1. In the `sdp_settings` passed to YALMIP, set the file name for MOSEK's log to be dumped to.
`sdp_settings.mosektaskfile = 'data/mosek_dump.ptf';`
2. Use MOSEK's Python console to load the dumped data and analyze it.
 - `pip install -r requirements.txt`
 - `python mosekconsole.py`

3. This will open MOSEK's interactive console. From there, you can use the following commands to analyze the problem:
 - `read data/mosek_dump.ptf`
 - `anapro`

 4. It is especially useful to look at the `|A|` output to see the ratio of coefficients of the lienar constraint matrix: [source](https://docs.mosek.com/11.0/toolbox/debugging-numerical.html)