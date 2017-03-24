% LLG Equation
% Oscar David Arbelaez Echeverri

# LLG Equation

The LLG equation reads,
$$
\frac{d\vec s}{dt} = - \frac{\gamma}{1 + \lambda^2}\left\lbrack
    \vec s \times \vec h_{eff} +
    \lambda \vec s \times \left( \vec s \times \vec h_{eff} \right)
\right\rbrack
$$

where  $\vec s$ is a unit vector representing a spin, $\gamma$ is the gyromagnetic
ratio, $\vec h_{eff}$ is the effective magnetic field acting upon the spin and
$\lambda$ is the Gilbert damping parameter.

The effective field can be computed out of a given hamiltonian by,
$$
\vec h_{eff} = - \frac{1}{\mu _s} \frac{\partial \mathcal{H}}{\partial \vec s}
$$

# Heun method for the LLG Equation

Now, a popular method to solve this equation is the Heun half step method, which
is a predictor corrector scheme, first you compute,
$$
\vec s ^ {\,\prime} = \vec s + \Delta \vec s \Delta t
$$
where,
$$
\Delta \vec s = - \frac{\gamma}{1 + \lambda^2}\left\lbrack
    \vec s \times \vec h_{eff} +
    \lambda \vec s \times \left( \vec s \times \vec h_{eff} \right)
\right\rbrack
$$

Here is the time to re normalize the new vector, since this scheme does not
preserve the norm, then recompute the field given the new spin prima, then apply
the corrector step,
$$
\vec s^{\,(t+\Delta t)} = \vec s^{\,(t)} + \frac{1}{2}\left\lbrack
    \Delta \vec s + \Delta \vec s ^ {\,\prime}
\right\rbrack
$$
where,
$$
\Delta \vec s^{\,\prime} = - \frac{\gamma}{1 + \lambda^2}\left\lbrack
    \vec s^{\,\prime} \times \vec h_{eff}^{\,\prime} +
    \lambda \vec s^{\,\prime} \times \left(
        \vec s^{\,\prime} \times \vec h_{eff}^{\,\prime}
    \right)
\right\rbrack
$$

And it's worthwhile to remember to normalize the spin after this step as well.
