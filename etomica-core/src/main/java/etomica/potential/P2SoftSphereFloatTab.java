/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;

/**
 * Soft-sphere interatomic potential.
 * Characterized by the pairwise-additive, spherically symmetric,
 * inverse-power potential of the form:
 * <p>
 * u(r) = epsilon*[sigma/ r]^n
 * <p>
 * where epsilon:  describes the strength of the pair interaction,
 * sigma  :  is the atom size parameter
 * n    :  is the degree of softness, s=1/n,
 * eg: hard sphere when s=0 or n->infinity
 * <p>
 * sigma/r  : sig_r
 * <p>
 * This class uses interpolation to handle the not-even-integer part of the
 * exponent.  The P2SoftSphere class will handle the even-integer-exponent
 * contribution (using the rounded-down even integer).  Interpolation of the
 * residual floating-point part [(sigma/r)^x] is handled as a function of r^2
 * such that sqrt is unneeded during the simulation.  Exponents close to even
 * are interpolated most efficiently because the interpolated function is
 * either very flat [x close to 0] or nearly linear [x close to 1].  Nearly odd
 * exponents are most difficult as x is nearly 1/2.
 *
 * This class will crash if asked about distances beyond the cutoff specified
 * at construction.
 */
public class P2SoftSphereFloatTab extends P2SoftSphere {

    public static Potential2Soft makeTruncated(Space space, double sigma, double epsilon, double nn, double rc, int ntab, TruncationFactory tf) {
        return tf.make(new P2SoftSphereFloatTab(space, sigma, epsilon, nn, rc, ntab));
    }

    public P2SoftSphereFloatTab(Space space, double sigma, double epsilon, double nn, double rc, int ntab) {
        super(space, sigma, epsilon, ((int) nn / 2) * 2);
        if (ntab > 100000) {
            System.err.println("Allocating for " + ntab + " values of r for tabulated potential is probably overkill");
        }
        nFloat = nn - n;
        xFac = ntab / (rc * rc);
        // push out 1 extra bin so we can return a value at rc
        rpTab = new double[ntab + 2][4];
        setSigma(sigma);
        setEpsilon(epsilon);
    }

    private void populateTabulatedValues() {
        if (rpTab == null) return;
        // compute extra bits beyond rc so we can accurately compute
        // derivatives up to rc
        for (int i = 0; i < rpTab.length; i++) {
            double r2 = i / xFac;
            // value
            rpTab[i][0] = Math.pow(r2, 0.5 * nFloat);
            // first order TODO
            rpTab[i][1] = 0.5 * nFloat * rpTab[i][0] / r2 / xFac;
        }

        // now set quadratic, cubic terms to enforce continuity
        for (int i = 0; i < rpTab.length - 1; i++) {
            rpTab[i][2] = 3 * (rpTab[i + 1][0] - rpTab[i][0]) - 2 * rpTab[i][1] - rpTab[i + 1][1];
            rpTab[i][3] = 2 * (rpTab[i][0] - rpTab[i + 1][0]) + rpTab[i][1] + rpTab[i + 1][1];
        }

    }

    private double rpInterp(double r2) {
        double x = r2 * xFac;
        int idx = (int) x;
        x -= idx;
        double[] irp = rpTab[idx];
        return irp[0] + x * (irp[1] + x * (irp[2] + x * irp[3]));
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        if (nFloat > 0) return super.u(r2) / rpInterp(r2);
        return super.u(r2) * rpInterp(r2);
    }

    public void u012add(double r2, double[] u012) {
        double u = u(r2);
        u012[0] += u;
        u012[1] += -(n + nFloat) * u;
        u012[2] += (n + nFloat) * (n + nFloat + 1) * u;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        return -(n + nFloat) * u(r2);
    }

    /**
     * The second derivative of the pair energy, times the square of the
     * separation:  r^2 d^2u/dr^2.
     */
    public double d2u(double r2) {
        return (n + nFloat) * (n + nFloat + 1) * u(r2);
    }

    /**
     * Mutator method for soft sphere size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public void setSigma(double sig) {
        super.setSigma(sig);
        if (epsilon > 0) populateTabulatedValues();
    }

    /**
     * Mutator method for soft-sphere energy parameter
     */
    public void setEpsilon(double eps) {
        super.setEpsilon(eps);
        if (sigma > 0) populateTabulatedValues();
    }

    public double integral(double rC) {
        double A = space.sphereArea(1.0);  //multiplier for differential surface element
        int D = space.D();                 //spatial dimension
        double rCD = space.powerD(rC);

        double sig_rCn = sigma / rC;
        return epsilon * A * rCD * Math.pow(sig_rCn, n+nFloat) / (n+nFloat - D);
    }

    private final double nFloat;
    private final double xFac;
    private final double[][] rpTab;
}
