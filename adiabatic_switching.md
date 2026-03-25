# Why does causality require adiabatic switching?

## The problem

The linearized EoM is a first-order inhomogeneous ODE:

$$i\hbar\,\partial_t \rho^{(1)} = \mathcal{L}[\rho^{(1)}] + o\,f(t)$$

Its general solution = particular solution + arbitrary homogeneous solution ($\rho^{(1)}_\text{hom} \propto e^{-i\omega_{cv}t}$, free particle-hole oscillations). We need a **boundary condition** to fix the solution uniquely.

The physical (causal) requirement is:

$$\rho^{(1)}(t \to -\infty) = 0$$

i.e., no response before the perturbation. But for a monochromatic field $f(t) = e^{-i\omega t}$, the perturbation extends to $t = -\infty$ and never vanishes — the boundary condition is **ill-defined**.

## The resolution

Replace $f(t) \to f(t)\,e^{\eta t}$ with $\eta \to 0^+$. Now $f(t) \to 0$ as $t \to -\infty$, so the system is genuinely unperturbed in the remote past and $\rho^{(1)}(-\infty) = 0$ is unambiguous.

Solving the ODE with this prescription:

- The particular solution acquires the denominator $\hbar(\omega + i\eta) - \Delta\xi$.
- The boundary condition $\rho^{(1)}(-\infty) = 0$ kills the homogeneous part entirely.

Taking $\eta \to 0^+$, the $i\eta$ survives as a **pole prescription**: all poles are pushed below the real $\omega$-axis, which is the defining property of the **retarded** Green's function.

## Summary

Adiabatic switching is not a physical assumption about how the field is turned on. It is a **mathematical device** to implement the causal boundary condition $\rho^{(1)}(-\infty) = 0$ in frequency space, where it would otherwise be ambiguous.

| Choice | Pole location | Response type |
|--------|--------------|---------------|
| $\eta \to 0^+$ | Below real axis | Retarded |
| $\eta \to 0^-$ | Above real axis | Advanced |
