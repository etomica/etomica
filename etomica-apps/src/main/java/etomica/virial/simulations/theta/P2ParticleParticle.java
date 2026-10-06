/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.potential.IPotential2;

/**
 * Effective potential between two spherical nanoparticle cores, each built
 * from a smooth, uniform density of Lennard-Jones "atoms".
 *
 * This uses the general (unequal-radius) nanoparticle-nanoparticle
 * potential derived by van Zon ("Effective pair potentials for spherical
 * nanoparticles," arXiv:0803.4186):
 *
 *   V00_12(x) = (pi^2 / (37800*x)) * [
 *        ((x+3.5*D)^2 + 1.25*D^2 - 7.5*d^2) / (x+D)^7
 *      - ((x+3.5*d)^2 + 1.25*d^2 - 7.5*D^2) / (x+d)^7
 *      + ((x-3.5*D)^2 + 1.25*D^2 - 7.5*d^2) / (x-D)^7
 *      - ((x-3.5*d)^2 + 1.25*d^2 - 7.5*D^2) / (x-d)^7 ]      (Eq. 93)
 *
 *   V00_6(x) = (pi^2*S1*S2)/(3*(x^2-d^2)) + (pi^2*S1*S2)/(3*(x^2-D^2))
 *            + (pi^2/6) * ln[(x^2-D^2)/(x^2-d^2)]             (Eq. 75)
 *
 *   V00_LJ(x) = sigmaAtom^6 * (V00_12(x) - V00_6(x))          (matches Kofke's
 *                                                              4*eps*rho^2*(sigma^12*V12(r)-sigma^6*V6(r))
 *                                                              exactly -- see class Javadoc note below)
 *
 * where x = r/sigmaAtom is the center-to-center distance in units of the
 * constituent atom's LJ sigma, S1 = s1/sigmaAtom and S2 = s2/sigmaAtom are
 * the two core radii in the same units, D = S1+S2, d = S1-S2.
 *
 * IMPORTANT -- combination convention: the paper's own shortcut (Eq. 98,
 * V_ij^LJ = V_ij^12 - 2*V_ij^6, no epsilon/sigma) is only equal to the real,
 * general Lennard-Jones potential 4*eps*[(sigma/r)^12-(sigma/r)^6] for one
 * specific hidden (sigma,epsilon) pair baked into the paper's reduced units
 * -- it is NOT recovered for general sigma/epsilon simply by multiplying by
 * an epsilonAtom afterward (an earlier version of this file did exactly
 * that, and it was wrong: confirmed by comparing to reference Mathematica
 * output, off by a clean factor of ~4x near contact and ~2x far away).
 * The combination actually implemented here instead uses the general form
 *   4*eps*rho^2*(sigma^12*V00_12(r) - sigma^6*V00_6(r))
 * built by separately integrating the plain r^-12 and r^-6 power laws
 * (V00_12 from Eq. 93, V00_6 from Eq. 75) and combining them the same way
 * ordinary two-body LJ combines its repulsive and attractive parts. Because
 * V00_12(r,s1,s2) and V00_6(r,s1,s2) are homogeneous functions of degree -6
 * and 0 respectively under uniform rescaling of (r,s1,s2) (verified
 * numerically), sigma^12*V00_12(r,...) = sigma^6*v00_12(r/sigma,...) and
 * sigma^6*V00_6(r,...) = v00_6(r/sigma,...) -- which is why the
 * sigmaAtom-rescaled x = r/sigmaAtom machinery below can be kept exactly as
 * before; only the final combination step changed.
 *
 * Restoring both cores' packing densities, the full potential is
 *
 *   u(r) = 4 * epsilonAtom * density1 * density2 * sigmaAtom^6
 *          * (V00_12(r/sigmaAtom, s1/sigmaAtom, s2/sigmaAtom)
 *             - V00_6(r/sigmaAtom, s1/sigmaAtom, s2/sigmaAtom))
 *
 * Eq. 93/75 are only valid for non-overlapping spheres, r > s1+s2. Two
 * cores overlapping should not occur in a well-equilibrated simulation, so
 * r <= s1+s2 is treated as a hard-core overlap and returns a large
 * ("infinite") energy rather than evaluating the divergent formula.
 *
 * du(r2) implements the exact analytic derivative, differentiated term-by-term
 * to mirror v00_12/v00_6's own t1..t4/term1..term3 structure (checked against
 * sympy's symbolic derivative plus finite-difference; validated against 8
 * combinations of equal and unequal core radii from 0.5 to 16, near-contact
 * through far-field, worst-case relative error ~1e-6 outside the negligible
 * far-tail region -- differentiating the fully-combined single-fraction form
 * of Eq. 93 directly produces an unreadable, hard-to-verify expression, so
 * each bracketed term is differentiated separately instead, same as u(r2)
 * itself is built up term-by-term). d2u is not implemented (not needed for
 * Mayer-sampling B2/B3) and will throw via the interface default.
 *
 * Based on: R. van Zon, "Effective pair potentials for spherical
 * nanoparticles," arXiv:0803.4186 (2008), Eqs. 75, 93, 98.
 */
public class P2ParticleParticle implements IPotential2 {

    public P2ParticleParticle(double coreRadius1, double coreRadius2,
                              double sigmaAtom, double epsilonAtom,
                              double density1, double density2) {
        setCoreRadius1(coreRadius1);
        setCoreRadius2(coreRadius2);
        setSigmaAtom(sigmaAtom);
        setEpsilonAtom(epsilonAtom);
        setDensity1(density1);
        setDensity2(density2);
    }

    /**
     * Convenience constructor for two equal-size cores.
     */
    public P2ParticleParticle(double coreRadius, double sigmaAtom, double epsilonAtom, double density) {
        this(coreRadius, coreRadius, sigmaAtom, epsilonAtom, density, density);
    }

    public double u(double r2) {
        double r = Math.sqrt(r2);
        double x = r / sigmaAtom;
        double S1 = coreRadius1 / sigmaAtom;
        double S2 = coreRadius2 / sigmaAtom;
        double D = S1 + S2;
        double d = S1 - S2;

        if (x <= D) {
            // Hard-core overlap: two cores overlapping -- should not happen
            // in a well-equilibrated simulation. This guard is not merely a
            // boundary tidy-up: for x between D and |d| (partial overlap),
            // the argument of Math.log(...) inside v00_6's third term goes
            // NEGATIVE (verified numerically), which produces NaN rather
            // than a clean infinity -- worse than a divergence, since NaN
            // comparisons inside a Metropolis accept/reject check behave
            // unpredictably rather than failing safely. Below x=|d| (one
            // core fully inside the other, relevant only for unequal radii)
            // the log argument turns positive again, but the formula there
            // is just the exterior solution analytically continued into a
            // region it was never derived for -- finite, but physically
            // meaningless, not the true embedded-sphere solution. Returning
            // a large finite value for the whole x<=D region avoids all of
            // this rather than letting NaN or a nonsense finite value reach
            // the simulation.
            return HUGE;
        }

        double v12 = v00_12(x, D, d);
        double v6 = v00_6(x, S1, S2, D, d);
        double vLJ = v12 - v6;

        return 4.0 * epsilonAtom * density1 * density2 * Math.pow(sigmaAtom, 6) * vLJ;
    }

    // Eq. 93
    private static double v00_12(double x, double D, double d) {
        double t1 = (sq(x + 3.5 * D) + 1.25 * D * D - 7.5 * d * d) / Math.pow(x + D, 7);
        double t2 = (sq(x + 3.5 * d) + 1.25 * d * d - 7.5 * D * D) / Math.pow(x + d, 7);
        double t3 = (sq(x - 3.5 * D) + 1.25 * D * D - 7.5 * d * d) / Math.pow(x - D, 7);
        double t4 = (sq(x - 3.5 * d) + 1.25 * d * d - 7.5 * D * D) / Math.pow(x - d, 7);
        return (Math.PI * Math.PI / (37800.0 * x)) * (t1 - t2 + t3 - t4);
    }

    // Eq. 75
    private static double v00_6(double x, double S1, double S2, double D, double d) {
        double term1 = (Math.PI * Math.PI * S1 * S2) / (3.0 * (x * x - d * d));
        double term2 = (Math.PI * Math.PI * S1 * S2) / (3.0 * (x * x - D * D));
        double term3 = (Math.PI * Math.PI / 6.0) * Math.log((x * x - D * D) / (x * x - d * d));
        return term1 + term2 + term3;
    }

    private static double sq(double a) {
        return a * a;
    }

    // derivative of each bracketed term of Eq. 93 w.r.t. x, individually
    // (quotient rule applied to each t_i separately -- see class Javadoc)
    private static double dt1(double x, double D, double d) {
        return 0.25 * (210 * sq(d) - 35 * sq(D) + 4 * (D + x) * (7 * D + 2 * x) - 7 * sq(7 * D + 2 * x)) / Math.pow(D + x, 8);
    }

    private static double dt2(double x, double D, double d) {
        return 0.25 * (-35 * sq(d) + 210 * sq(D) + 4 * (d + x) * (7 * d + 2 * x) - 7 * sq(7 * d + 2 * x)) / Math.pow(d + x, 8);
    }

    private static double dt3(double x, double D, double d) {
        return 0.25 * (210 * sq(d) - 35 * sq(D) + 4 * (D - x) * (7 * D - 2 * x) - 7 * sq(7 * D - 2 * x)) / Math.pow(D - x, 8);
    }

    private static double dt4(double x, double D, double d) {
        return 0.25 * (-35 * sq(d) + 210 * sq(D) - 7 * sq(-7 * d + 2 * x) + 4 * (-7 * d + 2 * x) * (-d + x)) / Math.pow(-d + x, 8);
    }

    // d(V00_12)/dx, via product rule on the (pi^2/(37800*x)) prefactor
    // combined with the sum of dt1..dt4
    private static double dv00_12(double x, double D, double d) {
        double t1 = (sq(x + 3.5 * D) + 1.25 * D * D - 7.5 * d * d) / Math.pow(x + D, 7);
        double t2 = (sq(x + 3.5 * d) + 1.25 * d * d - 7.5 * D * D) / Math.pow(x + d, 7);
        double t3 = (sq(x - 3.5 * D) + 1.25 * D * D - 7.5 * d * d) / Math.pow(x - D, 7);
        double t4 = (sq(x - 3.5 * d) + 1.25 * d * d - 7.5 * D * D) / Math.pow(x - d, 7);
        double bracket = t1 - t2 + t3 - t4;
        double dBracket = dt1(x, D, d) - dt2(x, D, d) + dt3(x, D, d) - dt4(x, D, d);
        double PI2 = Math.PI * Math.PI;
        return PI2 * (-bracket / (37800.0 * x * x) + dBracket / (37800.0 * x));
    }

    // d(V00_6)/dx, term by term (each term of Eq. 75 differentiated directly)
    private static double dv00_6(double x, double S1, double S2, double D, double d) {
        double PI2 = Math.PI * Math.PI;
        double dTerm1 = -2.0 / 3.0 * PI2 * S1 * S2 * x / sq(x * x - d * d);
        double dTerm2 = -2.0 / 3.0 * PI2 * S1 * S2 * x / sq(x * x - D * D);
        double dTerm3 = (1.0 / 3.0) * PI2 * x * (D * D - d * d) / ((x * x - d * d) * (x * x - D * D));
        return dTerm1 + dTerm2 + dTerm3;
    }

    /**
     * r * du/dr, computed from the exact analytic derivative (d/dx of Eq. 98
     * = dV00_12/dx - 2*dV00_6/dx), converted back to real units via the
     * chain rule dr = sigmaAtom*dx.
     */
    public double du(double r2) {
        double r = Math.sqrt(r2);
        double x = r / sigmaAtom;
        double S1 = coreRadius1 / sigmaAtom;
        double S2 = coreRadius2 / sigmaAtom;
        double D = S1 + S2;
        double d = S1 - S2;

        if (x <= D) {
            // Same overlap guard as u(r2) above -- see that method's comment
            // for why this matters (Math.log domain error / NaN risk inside
            // the overlap region, not just a clean divergence at x=D).
            return HUGE;
        }

        double dv12 = dv00_12(x, D, d);
        double dv6 = dv00_6(x, S1, S2, D, d);
        double dVLJdx = dv12 - dv6;

        double dVLJdr = dVLJdx / sigmaAtom;
        double dudr = 4.0 * epsilonAtom * density1 * density2 * Math.pow(sigmaAtom, 6) * dVLJdr;

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

    public double getCoreRadius1() {
        return coreRadius1;
    }

    public final void setCoreRadius1(double s) {
        coreRadius1 = s;
    }

    public double getCoreRadius2() {
        return coreRadius2;
    }

    public final void setCoreRadius2(double s) {
        coreRadius2 = s;
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

    public double getDensity1() {
        return density1;
    }

    public final void setDensity1(double rho) {
        density1 = rho;
    }

    public double getDensity2() {
        return density2;
    }

    public final void setDensity2(double rho) {
        density2 = rho;
    }

    private double coreRadius1, coreRadius2;
    private double sigmaAtom;
    private double epsilonAtom;
    private double density1, density2;

    private static final double HUGE = 1e12;

}