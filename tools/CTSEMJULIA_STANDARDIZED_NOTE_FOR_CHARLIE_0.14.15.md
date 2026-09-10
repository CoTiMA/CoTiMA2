# Note for Charlie: stationary standardized multilevel validation

Hi Charlie,

I built a complementary validation for the native Julia/Laplace path under an
idealized covariance-stationary latent model. Marginal latent variances are
fixed at one, study effects act on DRIFT, and DIFFUSION is derived through

`Q = -(A P + P A')`.

After correcting two issues in my own validation model (R matrix orientation
and an insufficiently constrained drift parameterization), the results are
encouraging. In a bivariate full design with 16 studies, 560 subjects, and
11,200 rows, six start configurations all reached the same solution. Two
simulation seeds recovered generating raw study SDs of 0.20 and 0.16 as
0.159/0.169 and 0.214/0.185. All four 95% intervals covered their generating
values. The Laplace-versus-five-node-quadrature gaps were 0.000278 and
0.000347 per subject. Both optimizations converged with gradient norms below
8e-9.

One API limitation may be useful for you to know about. The native Julia
backend currently rejects a parameter-dependent expression in T0VAR. This
prevents us from allowing a study-varying T0 correlation while keeping both
marginal T0 variances exactly at one through the Cholesky identity. I therefore
fixed the T0 correlation at zero in the final validation rather than fitting a
free Cholesky element that changes a marginal variance.

Expressions in DRIFT and DIFFUSION through PARS worked correctly. Is support
for parameter-dependent T0VAR entries planned, or is there another supported
parameterization you would recommend for a varying correlation with fixed unit
marginal variances?

The earlier non-positive-definite Q failures came from my provisional model
allowing incompatible combinations of A and P; I do not regard those as a
ctsem bug.
