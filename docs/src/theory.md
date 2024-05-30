# A brief example - The Rabi-model

To demonstrate the use of the method, we start with a simple but illustrative example -- the Rabi model. The Rabi model describes the interaction between a single linear cavity mode and a two-level (spin) system. However, unlike the Jaynes-Cummings model, it does not assume the rotating-wave approximation (RWA) and retains the counter-rotating terms. In the strong-coupling regime, the Rabi model is of special interest since the RWA breaks down as the spin-cavity coupling strength becomes comparable or greater than the mode frequencies\cite{Braumüller_USC_Rabi,Casanova_DSC_JC,Ashhab_USC_Rabi}. In order to showcase the time-coarse graining method, we focus on this strong-coupling regime in order to investigate the effects of the counter-rotating terms on the system dynamics when all measurements have limited time resolution.

More concretely, the model we consider here is described by the following Hamiltonian terms in the interaction picture,
The Hamiltonian terms in the interaction picture are described by the following equations:

\[
\begin{align}
\label{eq:Rabi Hamiltonian}
    \hat{H} &= \frac{g}{2} \left ( \hat{a}^\dagger \hat{\sigma}_+ e^{-i(\omega_c + \omega_a) t} + \hat{a}\hat{\sigma}_-e^{+i(\omega_c + \omega_a) t} \right ) \\
    &+ \frac{g}{2} \left ( \hat{a}^\dagger \hat{\sigma_-} e^{-i(\omega_c - \omega_a) t} + \hat{a}\hat{\sigma}_+e^{i(\omega_c - \omega_a) t} \right )
\end{align}
\]
where we $\omega_c$  $(\omega_a)$  is the cavity (atom) resonance.

In particular, we consider the situation where the coupling strength $g$ is not much smaller than the spin and cavity frequencies $\omega_{a}$ and $\omega_{c}$. i.e. $g \approx \omega_c \ (\omega_a)$. In that regime, the counter-rotating terms assume significance, and the induced dynamics that depend on the time resolution of the measurement apparatus.

When observing the dynamics of a coherent-state cavity mode interacting with a single atom, the Jaynes-Cummings model (which is simply the Rabi model with the RWA employed), show collapse-revival dynamics. Due to the photon-number dispersion of a coherent state, the Rabi oscillation decoheres and revives. Interestingly, these collapse-revivla cycles are completely absent in the full Rabi-model, and as we will see, only appear under finite-time resolution (or equivalently, time-coarse graining).

For example, the numerical simulations in this subsection assume the following set of parameters:
\[
\begin{align}
\frac{\omega_{c}}{2\pi} &= \frac{\omega_{a}}{2\pi} = 2 \textrm{GHz}; \\
\frac{g}{4\pi} &= 0.4 \textrm{GHz}.
\end{align}
\]

With the TCG method, we obtain an effective description that gives us \textbf{directly} the time-averaged observables that would be obtained from a bandwidth-limited measurement apparatus.  The TCG description produces a set of operators, comprised of products of the original Hamiltonian operators (i.e. multi-body transitions) $h^{(k)}_{\vec{\mu}}$, and their corresponding coupling strengths $g_{\vec{\mu}}^{(k)}$. In general, the produced TCG evolution does not need to be unitary, and produces also a set of pairs of pseudo-dissipators $(\hat{L}_{\vec\mu}, \hat{J}_{\vec\nu})$ with complex coupling rate $i\gamma_{\vec{\mu}, \vec{\nu}}^{(k)}$. These coupling strengths are encoded in a set of frequency-dependent scalars we call "contraction coefficients" $C_{l,r}(\vec{\mu}, \vec{\nu})$.
$$
    \begin{subequations}
        \begin{align}
        g^{(k)}_{\vec\mu}
        &
        =
        \sum_{\vec{\mu} \in \mathcal{P}_{k,0}[\Omega]} \frac{1}{2} \big( C_{k,0}(\vec{\mu}) + C_{k,0}(-\vec{\mu}^{\,\textrm{rev}}) \big) \\
        h^{(k)}_{\vec{\mu}}
        &
        =
        \prod_{\vec{\mu} \in \mathcal{P}_{k,0}[\Omega]} \hat{h}_{\mu_{k}} \hat{h}_{\mu_{k-1}} \cdots \hat{h}_{\mu_{1}} \\
        i \gamma^{(k)}_{\vec{\mu}, \vec{\nu}}
    &=
    \sum_{(\vec{\mu}, \vec{\nu}) \in \mathcal{P}_{l,r}[\Omega]} C_{l,r}(\vec{\mu}, \vec{\nu}) - C_{r,l}(-\vec{\nu}^{\,\textrm{rev}}, -\vec{\mu}^{\,\textrm{rev}}) \\
    \big( \hat{L}_{\vec{\mu}}, \hat{J}_{\vec{\nu}} \big)
    &
    =
    \Big\{ \big( \hat{h}_{\mu_{l}} \hat{h}_{\mu_{l-1}} \cdots \hat{h}_{\mu_{1}}, \hat{h}_{\nu_{1}} \hat{h}_{\nu_{2}} \cdots \hat{h}_{\nu_{r}} \big) \quad \vert \quad (\vec{\mu}, \vec{\nu}) \in \mathcal{P}_{l,r}[\Omega]  \Big\} 
        \end{align}
    \end{subequations}
$$
where $k = l + r$ is the order of up approximation in the original coupling strengths. We explain the details of the calculation in the following sections, but for now we will just assume these are given to us based on the original Hamiltonian, as illustrated in the figure, for example by the symbolic software package we developed \textbf{QuantumGraining.jl}, as illustrated in the figure. The full explicit calculation up to second-order is shown in the appendix \ref{app:rabi-model_contraction_coefficients}. 

In order to effectively capture the coarse-grained dynamics of the Rabi model, we apply the TCG perturbation theory up to the third order and derive the corresponding master equation. In particular, we will see that the TCG procedure reproduces the RWA Hamiltonian at the first-order, and goes beyond it starting at the second-order. 

The TCG prduces many different terms, but the high-frequency contributions would be exponentially suppressed by the filter function $f(\omega)$. Ignoring these exponentially suppressed terms, we find the most significant effective Hamiltonian terms to be,
$$
    \begin{equation}
    \hat{H}_{\rm TCG}^{(2)}
    \approx
    \hat{H}_{\rm{RWA}}
    +
    \frac{g^{2}}{8} \Big[
    \frac{1}{2\omega_{a}}
    -
    \big(
    \tau^{2}
    +
    \frac{1}{4\omega_{a}^{2}}
    \big) (\omega_{c} - \omega_{a})
    \Big]
    (1 + 2 \hat{a}^{\dagger} \hat{a}) \cdot \hat{\sigma}_z
    \end{equation}
$$
if we work in the limit where $\frac{1}{\omega_{a}} \ll \tau \ll \frac{1}{\abs{\omega_{c}-\omega_{a}}}, \frac{2}{g}$. 

The exact form of these corrections can be obtained using the QuantumGraining.jl package. Note that up to first-order, we get exactly the RWA Hamiltonian,
$$
\begin{equation}
\begin{split}
    \hat{H}_{\textrm{RWA}}
    =&
    \frac{g}{2}
    e^{-\frac{(\omega_{a}-\omega_{c})^{2}\tau^{2}}{2}}
    \Big(
    e^{i(\omega_{a}-\omega_{c})t} \hat{a} \hat{\sigma}_{+}
    +
    e^{-i(\omega_{a}-\omega_{c})t} \hat{a}^{\dagger} \hat{\sigma}_{-}
    \Big)\\
    \approx&
    \frac{g}{2}
    \Big(
    e^{i(\omega_{a}-\omega_{c})t} \hat{a} \hat{\sigma}_{+}
    +
    e^{-i(\omega_{a}-\omega_{c})t} \hat{a}^{\dagger} \hat{\sigma}_{-}
    \Big).
\end{split}
\end{equation}
$$
Note the appearance of a filter-dependent correction already at this order, which is an important difference from the standard RWA. Unlike the RWA, here the coarse-graining time-scale is a tunable parameter, to be determined by the experimental apparatus. Generally, the filter functions would be determined by the measurement apparatus, but for simplicity here we chose a Gaussian filter with a coarse-graining time-scale $\tau$.

A unique property of our TCG method is that it also produces non-unitary contributions to the dynamics. In second-order, we also get the following set of pseudo-dissipators in addition to the dispersive correction in the Hamiltonian,
$$
\begin{subequations}
    \begin{align}
    \hat{L}_1 &= \hat{J}_1 = \hat{a} \hat{\sigma}_+ &
    \gamma_1 &= - i \frac{g^2 \tau^2}{2} (\omega_c - \omega_a) e^{-2i(\omega_c - \omega_a) t} \\
    \hat{L}_2 &= \hat{J}_2 = \hat{a}^\dagger \hat{\sigma}_- &
    \gamma_2 &= i \frac{g^2 \tau^2}{2} (\omega_c - \omega_a) e^{2i(\omega_c - \omega_a) t}.
\end{align}
\end{subequations}
$$
These contributions are unique to the TCG method, and are not included in effective static Hamiltonian methods. Notice that the coefficients $\gamma_{1}$ and $\gamma_{2}$ are purely imaginary, so the corresponding pseudo-dissipators do not break the time reversal symmetry (hence the prefix ``pseudo-''). In addition, we do not expect significant secular effects from the second-order pseudo-dissipators since their coefficients are oscillatory at frequency $\pm 2(\omega_{c} - \omega_{a})$ and vanish in the resonant limit when $\omega_{c} \rightarrow \omega_{a}$. In particular, they do not induce any secular gain or loss of the system energy, but rather only add small fluctuating corrections (micro-oscillations) to the entropy and energy of the system.

However, this is not the case for higher-order corrections in general. For example, looking at the third order contributions,
$$
\begin{equation}
\begin{split}
    \lim_{\omega_{c} \rightarrow \omega_{a}}
    \hat{H}_{\textrm{TCG}}^{(3)}
    \approx
    \lim_{\omega_{c} \rightarrow \omega_{a}}
    \hat{H}_{\textrm{TCG}}^{(2)}
    -
    \frac{g^{3}}{32 \omega_{a}^{2}} \hat{a}^{\dagger} \hat{a} \hat{a}^{\dagger} \hat{\sigma}_{-}
    +
    h.c.
    \end{split}
    \end{equation}
    and
    \begin{equation}
    \begin{split}
    \lim_{\omega_{c} \rightarrow \omega_{a}}
    \hat{D}_{\textrm{TCG}}^{(3)}
    \approx
    \frac{i g^{3}}{32 \omega_{a}^{2}}
    \Big(
    \mathcal{D}[\hat{a}^{2} \hat{\sigma}_{z}, \hat{a}^{\dagger} \hat{\sigma}_{+}]
    -
    \mathcal{D}[\hat{a}^{\dagger 2} \hat{\sigma}_{z}, \hat{a} \hat{\sigma}_{-}]
    \Big)
    +
    h.c.
\end{split}
\end{equation}
$$
we see that both the unitary and non-unitary contributions both give rise to time-independent corrections to the Jaynes-Cummings model in the resonant limit. From the numerical simulation in Fig.\ref{fig:rabi-models} (obtained with a coarse-graining time scale of $\tau = 0.2 \textrm{ns}$), we see that they both have secular effects on the collapse and revival of spin-state population. In particular, although $\hat{D}_{\textrm{TCG}}^{(3)}$ does not cause any dissipation of energy over long periods of time due to its purely imaginary pre-factor, it does have observable secular effects on the \textbf{coherence }of the spin, which affects the collapse-revival pattern of the spin-state population. For example, both the RWA (first-order TCG) and the third-order TCG Hamiltonians make the false prediction of a double-revival pattern between $t=15\textrm{ns}$ and $t=35\textrm{ns}$, which is removed by including the third-order pseudo-dissipators in $\hat{D}_{\textrm{TCG}}^{(3)}$; in addition, adding those pseudo-dissipators also makes the resulting TCG dynamics closer to the one obtained by directly coarse-graining the exact dynamics, as can be observed in \ref{fig:rabi-models}.  It is also important to emphasize that the TCG master equations can be numerically simulated with much larger time steps without encountering any stiffness problems compared to the exact von-Neumann equation which contains fast-oscillating counter-rotating terms. Therefore, the TCG master equation not only offers analytical insights into the physics of an interacting system, but also allows much more efficient numerical simulation of the dynamics. 

In this example, we have treated the method as black-box, producing a new effective description that approximates the time-coarse grained dynamics. In the following sections, we will delve into how the method works in more detail.