/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

/**
 * Torsion potential expressed as power series in cos(theta)
 * https://doi.org/10.1063/1.474882
 */
public class P4BondTorsionPow implements IPotentialBondTorsion {

    protected double[] a;

    public P4BondTorsionPow(double[] a) {
        this.a = a;
    }

    @Override
    public double u(double costheta) {
        double rv = a[0];
        double ctp = costheta;
        for (int i=1; i<a.length; i++) {
            rv += a[i]*ctp;
            ctp *= costheta;
        }
        return rv;
    }

    @Override
    public void udu(double costheta, double[] u, double[] du) {
        u[0] = a[0];
        du[0] = 0;
        double ctp = 1;
        for (int i=1; i<a.length; i++) {
            du[0] += i*a[i]*ctp;
            ctp *= costheta;
            u[0] += a[i]*ctp;
        }
    }
}
