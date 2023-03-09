/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  The energy
 * is switched from fully-on to 0 over a short range (i.e. from 0.95*rC to rC).
 */
public class P2SoftSphericalSumTruncatedSwitched extends P2SoftSphericalSum {

    protected double rCutoff, r2Cutoff;
    protected int taperOrder;
    protected double switchFac, r2Switch;
    protected double a, b, c;

    public P2SoftSphericalSumTruncatedSwitched(double truncationRadius, IPotential2... potential) {
        super(potential);
        setTruncationRadius(truncationRadius);
        setTaperOrder(3);
        setSwitchFac(0.95);
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut * rCut;
        update();
    }

    /**
     * Accessor method for the radial cutoff distance.
     */
    public double getTruncationRadius() {
        return rCutoff;
    }

    public double getSwitchFac() {
        return switchFac;
    }

    public void setSwitchFac(double newSwitchFac) {
        switchFac = newSwitchFac;
        update();
    }

    public int getTaperOrder() {
        return taperOrder;
    }

    public void setTaperOrder(int newOrder) {
        taperOrder = newOrder;
        switch (taperOrder) {
            case 1:
                // linear from 0(rc) to 1(rswitch)
                // energy is continuous
                // force is discontinuous at switch and rc
                a = 1;
                b = c = 0;
                break;
            case 2:
                // parabola with max (0 slope) at rc, 0 to 1
                // force is continuous at rc, but not at switch
                a = 0;
                b = 1;
                c = 0;
                break;
            case 3:
                // cubic, 1 to 0; 0 slope at both start and end
                // force is continuous everywhere
                a = 0;
                b = 3;
                c = -2;
        }
    }

    /**
     * Returns the truncation radius.
     */
    public double getRange() {
        return rCutoff;
    }


    protected void update() {
        r2Switch = r2Cutoff * switchFac * switchFac;
    }

    protected double getF(double r) {
        double x = (rCutoff - r) / (rCutoff * (1 - switchFac));
        return a * x + b * x * x + c * x * x * x;
    }

    protected double getdFdr(double r) {
        // dF/dr = dF/dx dx/dr
        double x = (rCutoff - r) / (rCutoff * (1 - switchFac));
        double dFdx = a + 2 * b * x + 3 * c * x * x;
        double dxdr = -1 / (rCutoff * (1 - switchFac));
        return dFdx * dxdr;
    }

    protected double getd2Fdr2(double r) {
        // d2F/dr2 = d2F/dx2 (dx/dr)^2 + dF/dx d2x/dr2
        //         = d2F/dx2 (dx/dr)^2
        double x = (rCutoff - r) / (rCutoff * (1 - switchFac));
        double d2Fdx = 2 * b + 6 * c * x;
        double dxdr = -1 / (rCutoff * (1 - switchFac));
        return d2Fdx * dxdr * dxdr;
    }

    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public Dimension getTruncationRadiusDimension() {
        return Length.DIMENSION;
    }

    @Override
    public double u(double r2) {
        if (r2 > r2Cutoff) return 0;
        double u = uWrapped(r2);
        if (r2 < r2Switch) return u;
        return u * getF(Math.sqrt(r2));
    }

    @Override
    public double du(double r2) {
        if (r2 > r2Cutoff) return 0;

        double du = duWrapped(r2);

        if (r2 < r2Switch) return du;

        // U = u F
        // r (dU/dr) = r F (du/dr) + r u dF/dr
        //           = F du + r u dF/dr
        double u = uWrapped(r2);
        double r = Math.sqrt(r2);
        return getF(r) * du + r * u * getdFdr(r);
    }

    @Override
    public double d2u(double r2) {
        if (r2 > r2Cutoff) return 0;

        double d2u = d2uWrapped(r2);

        if (r2 < r2Switch) return d2u;

        // U = u F
        // r (dU/dr) = r F (du/dr) + r u dF/dr
        //           = F du + r u dF/dr
        // r2 (d2U/dr2) = r2 (dF/dr) (du/dr) + r2 F (d2u/dr2) + r2 (du/dr) (dF/dr) + r2 u d2F/dr2
        //           = r2 d2F/dr2 u + 2r (dF/dr) du + F d2u
        double u = uWrapped(r2);
        double du = duWrapped(r2);
        double r = Math.sqrt(r2);
        return r2 * getd2Fdr2(r) * u + 2 * r * getdFdr(r) * du + getF(r) * d2u;
    }

    @Override
    public double integral(Space space, double rC) {
        throw new RuntimeException("nope");
    }

    @Override
    public void u012add(double r2, double[] u012) {
        if (r2 > r2Cutoff) {
            return;
        }

        double[] thisU012 = new double[3];
        uduWrapped(r2, thisU012);
        if (r2 < r2Switch) {
            u012[0] += thisU012[0];
            u012[1] += thisU012[1];
            u012[2] += thisU012[2];
            return;
        }

        // U = u F
        // r (dU/dr) = r F (du/dr) + r u dF/dr
        //           = F du + r u dF/dr
        double r = Math.sqrt(r2);
        double F = getF(r);
        u012[0] += F * thisU012[0];
        u012[1] += F * thisU012[1] + r * thisU012[0] * getdFdr(r);
        u012[2] += r2 * getd2Fdr2(r) * thisU012[0] + 2 * r * getdFdr(r) * thisU012[1] + getF(r) * thisU012[2];
    }

    public static void main(String[] args) {
        P2LennardJones p = new P2LennardJones();
        P2SoftSphericalSumTruncatedSwitched pt = new P2SoftSphericalSumTruncatedSwitched(2, p);
        for (double x = 1.900001; x < 1.999999; x += 0.001) {
            System.out.println(x + " " + pt.getF(x) + " " + pt.getdFdr(x));
        }
        for (double x = 1.900001; x < 1.999999; x += 0.001) {
            System.out.println(x + " " + p.du(x * x) + " " + pt.du(x * x));
        }
    }
}
