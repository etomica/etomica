/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.amoeba;

import etomica.potential.IPotential2;

public class P2AmoebaChargeTransfer implements IPotential2 {

    protected final double a, b;

    public P2AmoebaChargeTransfer(double a, double b) {
        this.a = a;
        this.b = b;
    }

    @Override
    public double u(double r2) {
        double r = Math.sqrt(r2);
        return -a*Math.exp(-b*r);
    }

    @Override
    public double du(double r2) {
        double r = Math.sqrt(r2);
        return a*b*r*Math.exp(-b*r);
    }

    public double d2u(double r2) {
        double r = Math.sqrt(r2);
        return -a*b*b*r2*Math.exp(-b*r);
    }

    @Override
    public void u012add(double r2, double[] u012) {
        double r = Math.sqrt(r2);
        double aexpbr = a*Math.exp(-b*r);
        u012[0] = -aexpbr;
        u012[1] = b*aexpbr;
        u012[2] = -b*b*aexpbr;
    }
}
