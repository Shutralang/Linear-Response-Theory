# Time Dependent Hartree-Fock and Linear Response Theory

**Matthew Shu Liang** (Department of Physics, HKUST) & **Claude Opus 4.6** (Anthropic)

> This note reviews time dependent Hartree-Fock theory and how it can be used to calculate the linear response of a many-body system to an external perturbation. We specifically demonstrate that the TDHF response function has corrections from exchange interactions compared to the random phase approximation (RPA) response function.

---

## 1. Density Matrix Formalism of EoM

We start with a general many-body Hamiltonian:

$$\hat{H} = \sum_{ij}[t_{ij} + o_{ij} f(t)] \hat{c}^\dagger_i \hat{c}_j + \frac{1}{2} \sum_{ijkl} V_{ijkl} \hat{c}^\dagger_i \hat{c}^\dagger_j \hat{c}_k \hat{c}_l$$

where $o_{ij} f(t)$ is the external field perturbation. Hermiticity requires:

$$t_{ij} = t_{ji}^*, \quad o_{ij} = o_{ji}^*, \quad V_{ijkl} = V_{klij}^*$$

Fermionic anti-commutation requires:

$$V_{ijkl} = -V_{jikl} = -V_{ijlk} = V_{jilk}$$

The single-particle density matrix (time-dependent) is defined as:

$$\rho_{mn} = \langle \hat{c}^\dagger_n \hat{c}_m \rangle = \langle \Psi_0 | e^{iHt/\hbar} \hat{c}^\dagger_n \hat{c}_m e^{-iHt/\hbar} | \Psi_0 \rangle$$

The equation of motion is derived as:

$$i\hbar \partial_t \rho_{mn} = \langle [\hat{c}^\dagger_n \hat{c}_m, \hat{H}] \rangle$$
$$= \langle \hat{c}^\dagger_n [\hat{c}_m, \hat{H}] \rangle - \langle [\hat{H}, \hat{c}^\dagger_n] \hat{c}_m \rangle$$
$$= (t_{mj} + o_{mj} f(t)) \rho_{jn} - (t_{in} + o_{in} f(t)) \rho_{mi} + V_{mikl} \langle \hat{c}^\dagger_n \hat{c}^\dagger_i \hat{c}_k \hat{c}_l \rangle - V_{klni} \langle \hat{c}^\dagger_l \hat{c}^\dagger_k \hat{c}_i \hat{c}_m \rangle$$

where Einstein summation is implied. The commutator involving the interaction term:

$$[\hat{c}_m, \hat{V}] = \sum_{ikl} V_{mikl} \hat{c}^\dagger_i \hat{c}_k \hat{c}_l, \qquad [\hat{H}, \hat{c}^\dagger_n] = ([\hat{c}_n, \hat{H}])^\dagger$$

---

## 2. Time Dependent Hartree-Fock

Under Hartree-Fock approximation, the ground state is assumed to be a Slater determinant, so the EoM can be written in a self-consistent mean-field of $\rho$:

$$i\hbar \partial_t \rho_{mn} = [h(\rho) + o f(t),\; \rho]_{mn}$$

$$h(\rho)_{mn} \coloneqq t_{mn} + \sum_{ij}(V_{mijn} - V_{minj})\rho_{ji}$$

Assuming uniqueness of $\rho$, we can prove **idempotency** $\rho^2 = \rho$:

$$i\hbar \partial_t(\rho^2) = \rho[h+of,\rho] + [h+of,\rho]\rho = [h+of,\rho^2] \;\Rightarrow\; \rho^2 = \rho$$

When the external field is absent ($f(t)=0$), the system is in its static ground state and the mean-field EoM reduces to:

$$0 = [h(\rho^0), \rho^0]$$

The mean-field eigenequation ($f(t) = 0$):

$$\langle\hat{H}\rangle = t_{ij}\rho^0_{ji} + \frac{1}{2}\sum_{ijkl} V_{ijkl}(\rho^0_{li}\rho^0_{kj} - \rho^0_{ki}\rho^0_{lj})$$

$$\frac{\delta\langle\hat{H}\rangle}{\delta\langle\Psi_0|} = E^0|\Psi_0\rangle \;\Rightarrow\; \sum_{ij} h(\rho^0)_{ij}\hat{c}^\dagger_i \hat{c}_j |\Psi_0\rangle = E^0|\Psi_0\rangle$$

which solves $\rho^0$ as:

$$\rho^0_{\alpha\beta} = \sum_v \langle\alpha|v\rangle\langle v|\beta\rangle, \quad \text{where } |v\rangle \text{ is the occupied single-particle eigenstate}$$

$\rho^0$ is idempotent and acts as a projector onto the occupied subspace, with $\sigma^0 \coloneqq \mathbf{1} - \rho^0$. We expand $\rho$ as a perturbation series in $f(t)$:

$$\rho = \sum_{n=0} \rho^{(n)}$$

The linear response EoM:

$$i\hbar\partial_t \rho^{(1)} = [h(\rho^0), \rho^{(1)}] + [h'(\rho^0)\cdot\rho^{(1)}, \rho^0] + [of(t), \rho^0]$$

The idempotent property at first order implies:

$$\rho^0\rho^{(1)} + \rho^{(1)}\rho^0 = \rho^{(1)}$$
$$\Rightarrow\; \rho^{(1)}_{VV} = 0, \quad \rho^{(1)}_{CC} = 0$$

So $\rho^{(1)}$ only has matrix elements between valence and conduction subspaces. Projecting onto CV and VC blocks (with $v,c$ denoting valence/conduction eigenstates of $h(\rho^0)$):

$$i\hbar\partial_t \rho^{(1)}_{cv} = \left[(\xi_c - \xi_v)\delta_{cc'}\delta_{vv'} + V_{cv'c'v} - V_{cv'vc'}\right]\rho^{(1)}_{c'v'} + (V_{cc'v'v} - V_{cc'vv'})\rho^{(1)}_{v'c'} + o_{cv}f(t)$$

where

$$[h'(\rho^0)\cdot\rho^{(1)}]_{mn} = \sum_{ij} h'_{mn,ij}\rho^{(1)}_{ij}, \qquad \frac{\partial h_{mn}}{\partial \rho_{ij}} = V_{mjin} - V_{mjni}$$

We define the **dynamical matrix elements**:

$$\mathcal{E}_{cv,c'v'} = (\xi_c - \xi_v)\delta_{cc'}\delta_{vv'} + V_{cv'c'v} - V_{cv'vc'}, \qquad \Gamma_{cv,v'c'} = V_{cc'v'v} - V_{cc'vv'}$$

$$\mathcal{E}_{vc,v'c'} = -\mathcal{E}^*_{cv,c'v'}, \qquad \Gamma_{vc,c'v'} = -\Gamma^*_{cv,v'c'}$$

To select the causal (retarded) solution, we switch on the perturbation adiabatically from $t=-\infty$: $f(t) \to f(t)e^{\eta t}$ with $\eta \to 0^+$. Upon Fourier transform this amounts to $\omega \to \omega + i\eta$ in every energy denominator. The two EoMs can be written in frequency space as:

$$\begin{pmatrix} \hbar\omega & \\ & -\hbar\omega \end{pmatrix} \begin{pmatrix} \rho^{(1)}_{CV} \\ \rho^{(1)}_{VC} \end{pmatrix} = \begin{pmatrix} \mathcal{E} & \Gamma \\ \Gamma^* & \mathcal{E}^* \end{pmatrix} \begin{pmatrix} \rho^{(1)}_{CV} \\ \rho^{(1)}_{VC} \end{pmatrix} + \begin{pmatrix} o_{CV} & \\ & o_{VC} \end{pmatrix} f(\omega)$$

---

## 3. Application to Homogeneous Electron Gas

The HEG Hamiltonian is $\hat{H} = \hat{H}^0 + \hat{V}$, with:

$$\hat{V} = \frac{1}{2\Omega}\sum_{\mathbf{q}\mathbf{k}\mathbf{p}}\sum_{\alpha\beta} V(\mathbf{q})\,\hat{c}^\dagger_{\mathbf{k}\alpha}\hat{c}^\dagger_{\mathbf{p}\beta}\hat{c}_{\mathbf{p}+\mathbf{q}\,\beta}\hat{c}_{\mathbf{k}-\mathbf{q}\,\alpha}$$

$$V_{ijkl} = \frac{1}{2\Omega}\left[\underbrace{\delta_{s_is_l}\delta_{s_js_k}V(\mathbf{p}_i - \mathbf{p}_l)}_{\text{direct}} - \underbrace{\delta_{s_is_k}\delta_{s_js_l}V(\mathbf{p}_i - \mathbf{p}_k)}_{\text{exchange}}\right]\delta_{\mathbf{p}_i+\mathbf{p}_j,\,\mathbf{p}_k+\mathbf{p}_l}$$

$$\hat{H}^0 = \sum_{\alpha\mathbf{k}} \varepsilon_\mathbf{k}\,\hat{c}^\dagger_{\mathbf{k}\alpha}\hat{c}_{\mathbf{k}\alpha}, \qquad t_{ij} = \varepsilon_{\mathbf{p}_i}\delta_{s_is_j}\delta_{\mathbf{p}_i,\mathbf{p}_j}$$

The linear response of the density operator:

$$\delta\langle\hat{n}(\mathbf{q})\rangle(\omega) = \sum_{\mathbf{k}\alpha}\left([\rho^{(1)}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}]_{VC} + [\rho^{(1)}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}]_{CV}\right)$$

with $[\rho^{(1)}]_{VC} \propto \Theta(k_F - |\mathbf{k}+\mathbf{q}|)\Theta(|\mathbf{k}|-k_F)$ and $[\rho^{(1)}]_{CV} \propto \Theta(|\mathbf{k}+\mathbf{q}|-k_F)\Theta(k_F-|\mathbf{k}|)$.

The explicit EoM from the matrix block form is:

$$\hbar\omega[\rho^{(1)}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}]_{CV} = \mathcal{E}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha),(\mathbf{p}\,\beta)(\mathbf{l}\,\gamma)}[\rho^{(1)}_{(\mathbf{p}\,\beta)(\mathbf{l}\,\gamma)}]_{CV} + \Gamma_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha),(\mathbf{p}'\beta')(\mathbf{l}'\gamma')}[\rho^{(1)}_{(\mathbf{p}'\beta')(\mathbf{l}'\gamma')}]_{VC} + o_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}f(\omega)$$

### Evaluating $\mathcal{E}$ for HEG

Mapping: $c = (\mathbf{k}+\mathbf{q},\alpha)$, $v = (\mathbf{k},\alpha)$, $c' = (\mathbf{p},\beta)$, $v' = (\mathbf{l},\gamma)$, with $|\mathbf{k}+\mathbf{q}|,|\mathbf{p}| > k_F$ and $|\mathbf{k}|,|\mathbf{l}| < k_F$.

**Term $V_{cv'c'v}$:** Setting $i=(\mathbf{k}+\mathbf{q},\alpha)$, $j=(\mathbf{l},\gamma)$, $k=(\mathbf{p},\beta)$, $l=(\mathbf{k},\alpha)$:
- Momentum conservation: $\mathbf{l} = \mathbf{p} - \mathbf{q}$
- Direct: $\delta_{\gamma\beta}V(\mathbf{q})$
- Exchange: $\delta_{\alpha\beta}\delta_{\alpha\gamma}V(\mathbf{k}+\mathbf{q}-\mathbf{p})$

$$V_{cv'c'v} = \frac{1}{2\Omega}[\delta_{\beta\gamma}V(\mathbf{q}) - \delta_{\alpha\beta}\delta_{\alpha\gamma}V(\mathbf{k}+\mathbf{q}-\mathbf{p})]\,\delta_{\mathbf{l},\mathbf{p}-\mathbf{q}}$$

**Term $V_{cv'vc'}$:** Setting $i=(\mathbf{k}+\mathbf{q},\alpha)$, $j=(\mathbf{l},\gamma)$, $k=(\mathbf{k},\alpha)$, $l=(\mathbf{p},\beta)$:
- Same momentum conservation
- Direct: $\delta_{\alpha\beta}\delta_{\gamma\alpha}V(\mathbf{k}+\mathbf{q}-\mathbf{p})$
- Exchange: $\delta_{\gamma\beta}V(\mathbf{q})$

$$V_{cv'vc'} = \frac{1}{2\Omega}[\delta_{\alpha\beta}\delta_{\alpha\gamma}V(\mathbf{k}+\mathbf{q}-\mathbf{p}) - \delta_{\beta\gamma}V(\mathbf{q})]\,\delta_{\mathbf{l},\mathbf{p}-\mathbf{q}}$$

**Result:** The two terms are exact negatives, so they double upon subtraction:

$$\boxed{\mathcal{E}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha),(\mathbf{p}\,\beta)(\mathbf{l}\,\gamma)} = (\xi_{\mathbf{k}+\mathbf{q}} - \xi_\mathbf{k})\delta_{\mathbf{p},\mathbf{k}+\mathbf{q}}\delta_{\mathbf{l},\mathbf{k}}\delta_{\alpha\beta}\delta_{\alpha\gamma} + \frac{1}{\Omega}\left[\delta_{\beta\gamma}V(\mathbf{q}) - \delta_{\alpha\beta}\delta_{\alpha\gamma}V(\mathbf{k}+\mathbf{q}-\mathbf{p})\right]\delta_{\mathbf{l},\mathbf{p}-\mathbf{q}}}$$

where $\xi_\mathbf{k} = \varepsilon_\mathbf{k} + \Sigma^{\text{HF}}(\mathbf{k})$ is the HF single-particle energy.

**Notes:**
1. Momentum conservation forces $\mathbf{l} = \mathbf{p} - \mathbf{q}$, so only particle-hole pairs with the same $\mathbf{q}$ couple; different $\mathbf{q}$-sectors decouple.
2. The direct term $\frac{1}{\Omega}\delta_{\beta\gamma}V(\mathbf{q})$ depends only on $\mathbf{q}$. Retaining only this term recovers **RPA**.
3. The exchange term $-\frac{1}{\Omega}\delta_{\alpha\beta}\delta_{\alpha\gamma}V(\mathbf{k}+\mathbf{q}-\mathbf{p})$ is the **TDHF correction beyond RPA**.

### Evaluating $\Gamma$ for HEG

Mapping: $c=(\mathbf{k}+\mathbf{q},\alpha)$, $v=(\mathbf{k},\alpha)$, $v'=(\mathbf{p}',\beta')$, $c'=(\mathbf{l}',\gamma')$, with $|\mathbf{p}'|<k_F$ and $|\mathbf{l}'|>k_F$.

By the same procedure (the two terms $V_{cc'v'v}$ and $V_{cc'vv'}$ are again exact negatives):

$$\boxed{\Gamma_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha),(\mathbf{p}'\beta')(\mathbf{l}'\gamma')} = \frac{1}{\Omega}\left[\delta_{\beta'\gamma'}V(\mathbf{q}) - \delta_{\alpha\beta'}\delta_{\alpha\gamma'}V(\mathbf{k}+\mathbf{q}-\mathbf{p}')\right]\delta_{\mathbf{l}',\mathbf{p}'-\mathbf{q}}}$$

$\Gamma$ has the same interaction structure as the off-diagonal part of $\mathcal{E}$ (with $\mathbf{p}\to\mathbf{p}'$, $\beta\to\beta'$, $\gamma\to\gamma'$), but without the diagonal quasiparticle energy term.

### RPA Equation of Motion

Dropping exchange terms:

$$\mathcal{E}^{\text{RPA}} = (\xi_{\mathbf{k}+\mathbf{q}} - \xi_\mathbf{k})\delta_{\mathbf{p},\mathbf{k}+\mathbf{q}}\delta_{\mathbf{l},\mathbf{k}}\delta_{\alpha\beta}\delta_{\alpha\gamma} + \frac{1}{\Omega}\delta_{\beta\gamma}V(\mathbf{q})\delta_{\mathbf{l},\mathbf{p}-\mathbf{q}}$$

$$\Gamma^{\text{RPA}} = \frac{1}{\Omega}\delta_{\beta'\gamma'}V(\mathbf{q})\delta_{\mathbf{l}',\mathbf{p}'-\mathbf{q}}$$

**$\mathcal{E}$ contraction:** The diagonal part collapses via Kronecker deltas. The direct part, after eliminating $\mathbf{l}=\mathbf{p}-\mathbf{q}$ and $\beta=\gamma$, gives $\frac{V(\mathbf{q})}{\Omega}\sum_{\mathbf{k}'\beta}[\rho^{(1)}_{(\mathbf{k}'+\mathbf{q}\,\beta)(\mathbf{k}'\beta)}]_{CV}$.

**$\Gamma$ contraction:** Similarly gives $\frac{V(\mathbf{q})}{\Omega}\sum_{\mathbf{k}'\beta'}[\rho^{(1)}_{(\mathbf{k}'+\mathbf{q}\,\beta')(\mathbf{k}'\beta')}]_{VC}$.

The two sums reconstruct the density response:

$$\delta\langle\hat{n}(\mathbf{q})\rangle(\omega) = \sum_{\mathbf{k}'\beta}\left([\rho^{(1)}_{(\mathbf{k}'+\mathbf{q}\,\beta)(\mathbf{k}'\beta)}]_{CV} + [\rho^{(1)}_{(\mathbf{k}'+\mathbf{q}\,\beta)(\mathbf{k}'\beta)}]_{VC}\right)$$

giving the **RPA EoM**:

$$\boxed{\hbar\omega\;[\rho^{(1)}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}]_{CV} = (\xi_{\mathbf{k}+\mathbf{q}} - \xi_\mathbf{k})[\rho^{(1)}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}]_{CV} + \frac{V(\mathbf{q})}{\Omega}\delta\langle\hat{n}(\mathbf{q})\rangle(\omega) + o_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}f(\omega)}$$

The coupling to all other particle-hole pairs is entirely through $\delta\langle\hat{n}(\mathbf{q})\rangle(\omega)$ — the hallmark of RPA screening.

### Solving the RPA EoM

Using Hermiticity $\rho^{(1)}_{vc}(\omega) = (\rho^{(1)}_{cv}(-\omega))^*$ to obtain the VC block, the solutions are:

$$\boxed{[\rho^{(1)}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}]_{CV} = \frac{\frac{V(\mathbf{q})}{\Omega}\delta\langle\hat{n}(\mathbf{q})\rangle(\omega) + o_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}f(\omega)}{\hbar(\omega+i\eta) - (\xi_{\mathbf{k}+\mathbf{q}} - \xi_\mathbf{k})}(1-f_{\mathbf{k}+\mathbf{q}})f_\mathbf{k}}$$

$$\boxed{[\rho^{(1)}_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}]_{VC} = \frac{\frac{V(\mathbf{q})}{\Omega}\delta\langle\hat{n}(\mathbf{q})\rangle(\omega) + o_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)}f(\omega)}{-\hbar(\omega+i\eta) - (\xi_\mathbf{k} - \xi_{\mathbf{k}+\mathbf{q}})}f_{\mathbf{k}+\mathbf{q}}(1-f_\mathbf{k})}$$

where $f_\mathbf{k} = \Theta(\varepsilon_F - \xi_\mathbf{k})$. Taking $o_{(\mathbf{k}+\mathbf{q}\,\alpha)(\mathbf{k}\,\alpha)} = o_\mathbf{q}$ and self-consistently solving:

$$\Pi^0(\mathbf{q},\omega) = \sum_\mathbf{k} \frac{(1-f_{\mathbf{k}+\mathbf{q}})f_\mathbf{k} - f_{\mathbf{k}+\mathbf{q}}(1-f_\mathbf{k})}{\hbar(\omega+i\eta) - (\xi_{\mathbf{k}+\mathbf{q}} - \xi_\mathbf{k})}$$

$$\delta\langle\hat{n}(\mathbf{q})\rangle(\omega) = \frac{\Pi^0(\mathbf{q},\omega)}{1 - \frac{V(\mathbf{q})}{\Omega}\Pi^0(\mathbf{q},\omega)}\;o_\mathbf{q}\,f(\omega)$$

The $i\eta$ prescription makes $\Pi^0$ the **retarded** response function, whose poles lie below the real $\omega$-axis, ensuring causality. $\Pi^0(\mathbf{q},\omega)$ has the form of the Lindhard function, except the energies are from the HF spectrum:

$$\xi_{\mathbf{k}\alpha} = \varepsilon_\mathbf{k} + \frac{1}{\Omega}\sum_{\mathbf{q}\gamma}\left[V(0) - \delta_{\alpha\gamma}V(\mathbf{k}-\mathbf{q})\right]\rho^0_{\mathbf{q}\gamma}$$

The diverging $V(0)$ is cancelled by the positive background charge, and ignoring the remaining exchange term recovers the original RPA density-density response function.

---

## References

- Fetter & Walecka, *Quantum Theory of Many-Particle Systems*
- Shao (2026)
