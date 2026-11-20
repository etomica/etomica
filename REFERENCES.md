# References

Publications describing algorithms implemented in Etomica, mapped to the code that implements them.

- Entries are keyed by DOI, so they work for anyone. `PDF:` is an optional local path.
  David's papers are in `~/Documents/Papers Library/`.
- `Code:` lists the main implementing classes, not every class that uses them.
- `Notes:` records how the code relates to the paper, especially where it departs from or extends it.

## Virial coefficients (`etomica.virial`)

### Recursive evaluation of biconnected Mayer clusters

- R. J. Wheatley, "Calculation of High-Order Virial Coefficients with Applications to Hard and Soft
  Spheres," Phys. Rev. Lett. **110**, 200601 (2013). doi:[10.1103/PhysRevLett.110.200601](https://doi.org/10.1103/PhysRevLett.110.200601)
  - PDF: `~/Documents/Papers Library/2013-Wheatley-Wheatley-Physical Review Letters-Calculation of High-Order Virial Coefficients with Applications to Hard and Soft Spheres.pdf`
  - Code: `etomica.virial.wheatley` package. `ClusterWheatleySoft` (general pair potentials),
    `ClusterWheatleyHS` (hard spheres, with clique screening), `ClusterWheatleyPartitionScreening`,
    plus the `*Mix`, `*BD`, `*Derivatives`, `*Multibody`, `ExtendSW` and `PT` variants.
  - Notes: computes the sum over all biconnected graphs (f_B) in O(3^N) time using recursions
    f_Q → f_C → f_B (Eqs. 4–8), rather than enumerating graphs. For hard spheres, the paper screens
    out zero-f_B configurations with a depth-first-search biconnectivity test and simple clique
    separators (m = 2–5). The paper handles round-off by zeroing |f_B| below ~1e-14. Instead,
    `ClusterWheatleySoft(n, f, tol)` recomputes small values (|f_B| < 100·tol) with the BigDecimal
    `*BD` variant. The paper's own sampling scheme (chain placement,
    permutation-averaged P) is not what Etomica uses: Etomica uses Mayer-sampling/overlap MC instead.
    The derivative, multibody and temperature-perturbation extensions go beyond this paper; see
    Wheatley et al. (2020) below.

- R. J. Wheatley, A. J. Schultz, H. Do, N. Gokul, and D. A. Kofke, "Cluster integrals and virial
  coefficients for realistic molecular models," Phys. Rev. E **101**, 051301(R) (2020).
  doi:[10.1103/PhysRevE.101.051301](https://doi.org/10.1103/PhysRevE.101.051301)
  - PDF: `~/Documents/Papers Library/2020-Wheatley-Kofke-Physical Review B-Cluster integrals and virial coefficients for realistic molecular models_1.pdf`
    (the filename says "Physical Review B", but the paper is in PRE)
  - Code: `etomica.virial.wheatley.ClusterWheatleyMultibody` (non-additive potentials),
    `ClusterWheatleySoftDerivatives` (β-derivatives), `ClusterWheatleyMultibodyDerivatives`
    (both), plus their `Mix` and `BD` variants. Drivers: `etomica.virial.simulations.VirialH2OGCPMD`
    (GCPM water, B_N with temperature derivatives, the paper's application),
    `simulations.helium.VirialHeNonAdditiveD`, `VirialCO2H2OGCPM`, `VirialN2PI`.
  - Notes: generalizes the 2013 recursion to arbitrary (non-pairwise) energies. It starts from
    f_Q(N) = exp[-βE(N)] and uses f_Q → f_C (Eq. 4) → f_B via the cluster-star recursion (Eqs. 7a–d),
    with no reference to pair bonds or biconnected graphs. In `ClusterWheatleyMultibody.calcFullFQ`,
    f_Q for each subset is built as the pairwise product times exp(-βΔU) of the non-additive part,
    taken from each sub-subset of 3 or more molecules. Temperature derivatives propagate through
    the same recursion with binomial coefficients (Eqs. 8–10); `ClusterWheatleySoftDerivatives`
    keeps a `binomial` table for this. The paper's Supplemental Material has the proofs, timings
    for N = 4–12 and a link to an implementation.
