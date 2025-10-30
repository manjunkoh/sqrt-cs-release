# Square root QR covariance steering

# TODO
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