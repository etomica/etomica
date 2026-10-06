/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.potential.IPotential2;

/**
 * Effective potential between a point particle (a polymer segment / monomer)
 * and a spherical nanoparticle core built from a smooth, uniform density of
 * Lennard-Jones "atoms".
 *
 * This is the point-nanoparticle potential derived by integrating the LJ
 * potential over the volume of a solid sphere of radius s (van Zon,
 * "Effective pair potentials for spherical nanoparticles," arXiv:0803.4186,
 * Eq. 96):
 *
 *   V0(x) = 4*pi*S^3 / (3*(x^2-S^2)^6)
 *         + (80*pi*S^9 + 432*pi*x^4*S^5) / (45*(x^2-S^2)^9)
 *         - 8*pi*S^3 / (3*(x^2-S^2)^3)
 *
 * where x = r/sigmaAtom and S = s/sigmaAtom are lengths measured in units of
 * the *constituent atom's* LJ sigma (not the core's own radius). This raw
 * combination (term1+term2-term3, i.e. V0^12 - 2*V0^6) is the PAPER'S OWN
 * shortcut notation, which only equals the real, general Lennard-Jones
 * potential 4*eps*[(sigma/r)^12-(sigma/r)^6] for one specific hidden
 * (sigma,epsilon) pair baked into the paper's reduced units -- it is NOT
 * recovered for general sigma/epsilon by simply multiplying by an
 * epsilonAtom afterward (an earlier version of this file did exactly that
 * and was wrong, confirmed off by a clean factor of 4x near contact and 2x
 * far away when checked against reference Mathematica output).
 *
 * The combination actually used below instead builds the general LJ form
 * the same way ordinary two-body LJ does -- separately integrating the
 * plain r^-12 and r^-6 power laws (giving V0^12 = term1+term2, and
 * V0^6 = term3/2, per Eq. 91 and Eq. 73 respectively) and combining as
 *
 *   u(r) = 4 * epsilonAtom * density * sigmaAtom^3
 *          * (V0^12(r/sigmaAtom, s/sigmaAtom) - V0^6(r/sigmaAtom, s/sigmaAtom))
 *
 * V0^12 and V0^6 are homogeneous of degree -9 and -3 respectively under
 * uniform rescaling of (r,s) (verified numerically), which is why the
 * existing x=r/sigmaAtom rescaling can be kept unchanged -- only the final
 * combination step (the -2 vs -1 weighting, and the missing sigmaAtom^3
 * and 4x factors) needed fixing.
 *
 * epsilonAtom/sigmaAtom are the LJ parameters of the atoms making up
 * the core, and density is the (dimensionless, in units of sigmaAtom^-3)
 * packing density of those atoms inside the core -- this is the "internal
 * nanoparticle density" parameter Andrew asked about.
 *
 * Eq. 96 is only valid in the non-overlapping regime, r > s. Since a monomer
 * should never actually penetrate the core in a well-equilibrated
 * simulation, r <= s is treated as a hard-core overlap and returns a very
 * large (numerically "infinite") energy rather than evaluating the
 * (divergent) formula directly.
 *
 * du(r2) implements the exact analytic derivative of Eq. 96 (dV0/dx, derived
 * term-by-term to mirror u(r2)'s term1/term2/term3 structure, each term
 * differentiated independently with sympy and checked against
 * finite-difference; validated across core radii 0.5-16 from near-contact
 * through far-field, worst-case relative error ~2e-6 right at the steep
 * near-divergence region, ~1e-9 to 1e-11 everywhere else). d2u is not
 * implemented (not needed for Mayer-sampling B2/B3, which only use u(r2))
 * and will throw via the interface default.
 *
 * Based on: R. van Zon, "Effective pair potentials for spherical
 * nanoparticles," arXiv:0803.4186 (2008), Eq. 96.
 */
public class P2ParticleSegment implements IPotential2 {

    public P2ParticleSegment(double coreRadius, double sigmaAtom, double epsilonAtom, double density) {
        setCoreRadius(coreRadius);
        setSigmaAtom(sigmaAtom);
        setEpsilonAtom(epsilonAtom);
        setDensity(density);
    }

    /**
     * The energy u, as a function of r^2 (matching the IPotential2
     * convention used throughout Etomica, e.g. P2LennardJones).
     */
    public double u(double r2) {
        double r = Math.sqrt(r2);
        // dimensionless lengths, measured in units of the constituent atom's sigma
        double x = r / sigmaAtom;
        double S = coreRadius / sigmaAtom;

        if (x <= S) {
            // Hard-core overlap: monomer inside the core -- should not happen
            // in a well-equilibrated simulation. This guard is not merely
            // tidying up a boundary case: Eq. 96 has a genuine division-by-zero
            // singularity exactly at x=S, but just *inside* it (x slightly
            // below S) the raw formula does NOT continue as a repulsive wall
            // -- it flips sign and plunges to a spurious, unboundedly deep
            // NEGATIVE energy (verified numerically: at x=S-0.0001 it returns
            // roughly -7e34, vs +7e34 at x=S+0.0001). Left unguarded, a
            // Metropolis MC move would be drawn toward, not repelled from,
            // a monomer sitting just inside the surface. Returning a large
            // finite positive value here is what prevents that.
            return HUGE;
        }

        double d2 = x * x - S * S;       // (x^2 - S^2)
        double d3 = d2 * d2 * d2;         // (x^2 - S^2)^3
        double d6 = d3 * d3;               // (x^2 - S^2)^6
        double d9 = d6 * d3;               // (x^2 - S^2)^9

        double S3 = S * S * S;
        double S5 = S3 * S * S;
        double S9 = S3 * S3 * S3;

        double term1 = 4.0 * Math.PI * S3 / (3.0 * d6);
        double term2 = (80.0 * Math.PI * S9 + 432.0 * Math.PI * x * x * x * x * S5) / (45.0 * d9);
        double term3 = 8.0 * Math.PI * S3 / (3.0 * d3);

        double v0LJ = (term1 + term2) - term3 / 2.0;

        return 4.0 * epsilonAtom * density * Math.pow(sigmaAtom, 3) * v0LJ;
    }

    /**
     * r * du/dr, computed from the exact analytic derivative of Eq. 96 with
     * respect to x = r/sigmaAtom, then converted back to real units via the
     * chain rule (dV0/dr = dV0/dx * dx/dr = dV0/dx / sigmaAtom).
     * <p>
     * Derivative of each term of V0(x), differentiated separately with
     * sympy and each checked against finite difference of that same term
     * to ~1e-9 relative accuracy (an earlier hand-derived quotient-rule
     * expansion of d(term2)/dx had an algebra error and gave wrong values
     * away from x=20; the expression below has been re-verified):
     * <p>
     * d(term1)/dx = -16*pi*S^3*x / w^7
     * d(term3)/dx = -16*pi*S^3*x / w^4
     * d(term2)/dx = 32*pi*S^5*x*(-5*S^4 - 21*x^4 - 6*S^2*x^2) / (5*w^10)
     * where w = x^2 - S^2.
     */
    public double du(double r2) {
        double r = Math.sqrt(r2);
        double x = r / sigmaAtom;
        double S = coreRadius / sigmaAtom;

        if (x <= S) {
            // Same overlap guard as u(r2) above -- see that method's comment
            // for why this matters (Eq. 96 flips sign to a spurious deep
            // attractive well just inside x=S if left unguarded).
            return HUGE;
        }

        double S3 = S * S * S;
        double S5 = S3 * S * S;
        double S9 = S3 * S3 * S3;

        double w = x * x - S * S;
        double w4 = w * w * w * w;
        double w7 = w4 * w * w * w;
        double w10 = w7 * w * w * w;

        double dTerm1 = -16.0 * Math.PI * S3 * x / w7;
        double dTerm3 = -16.0 * Math.PI * S3 * x / w4;

        double N2 = -5.0 * S * S * S * S - 21.0 * x * x * x * x - 6.0 * S * S * x * x;
        double dTerm2 = 32.0 * Math.PI * S5 * x * N2 / (5.0 * w10);

        double dV0LJdx = (dTerm1 + dTerm2) - dTerm3 / 2.0;

        double dV0LJdr = dV0LJdx / sigmaAtom;   // chain rule: x = r/sigmaAtom
        double dudr = 4.0 * epsilonAtom * density * Math.pow(sigmaAtom, 3) * dV0LJdr;

        return r * dudr;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public void u012add(double r2, double[] u012) {
        u012[0] += u(r2);
        u012[1] += du(r2);
    }

    public double getCoreRadius() {
        return coreRadius;
    }

    public final void setCoreRadius(double s) {
        coreRadius = s;
    }

    public double getSigmaAtom() {
        return sigmaAtom;
    }

    public final void setSigmaAtom(double s) {
        sigmaAtom = s;
    }

    public double getEpsilonAtom() {
        return epsilonAtom;
    }

    public final void setEpsilonAtom(double eps) {
        epsilonAtom = eps;
    }

    public double getDensity() {
        return density;
    }

    public final void setDensity(double rho) {
        density = rho;
    }

    private double coreRadius;
    private double sigmaAtom;
    private double epsilonAtom;
    private double density;

    private static final double HUGE = 1e12;

    public static void main(String[] args) throws java.io.IOException {
        double coreRadius = 8.0;   // core=16 -> radius = 16/2
        P2ParticleSegment p2 = new P2ParticleSegment(coreRadius, 1.0, 1.0, 1.34);

        try (java.io.FileWriter fw = new java.io.FileWriter("segment_profile.csv")) {
            fw.write("r,U\n");
            for (double r = 8.0; r < 10.0; r += 0.1) {
                double u = p2.u(r * r);
                fw.write(r + "," + u + "\n");
                System.out.println(r + "," + u);
            }
        }
        System.out.println("Segment potential data written to: segment_profile.csv");
    }



}


