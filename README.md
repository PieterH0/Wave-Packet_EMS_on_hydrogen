# Wave-Packet_EMS_on_hydrogen
Computation of wave-packet EMS scattering probability in the symmetric noncoplanar geometry

This program calculates the exact (within the first Born and PWIA approximations) scattering probability for an electron wave-packet on a fixed hydrogen atom. The analytic expression is the following
$$ \frac{\partial^2\mathcal{P}}{\partial \Omega_a \partial\Omega_b}(\phi) = \int d k_f \  k_f^3 \Big( \frac{1}{4} |\mathcal{A}(k_f, \Omega_a, \Omega_b, s=0)|^2 + \frac{3}{4} |\mathcal{A}(k_f, \Omega_a, \Omega_b, s=1)|^2 \Big) $$
with the scattering amplitude given by
$$ \mathcal{A}_{\bm k_a \bm k_b i} &= 2\pi i \sum_{n,l,m}  c_{n,l,m} k_n \int d\Omega_i \ a_e(k_n \bm{\hat{k}}_i) \frac{1}{(2\pi)^3} \frac{4\pi}{\Delta^2} \phi_{n,l,m}(\bm q) e^{-i k_n\bm{\hat{k}}_i \cdot \bm b}$$
