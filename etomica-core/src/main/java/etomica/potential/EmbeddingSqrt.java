/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

public class EmbeddingSqrt implements IPotentialEmbedding {

    private final double epsilon;

    public EmbeddingSqrt(double epsilon) {
        this.epsilon = epsilon;
    }

    @Override
    public double u(double rho) {
        return -epsilon * Math.sqrt(rho);
    }

    @Override
    public void udud2u(double rho, double[] f, double[] df, double[] d2f) {
        double s = Math.sqrt(rho);
//        f[0]   = -epsilon * s;
//        df[0]  = -0.5 * epsilon / s;  //df/drho
//        d2f[0] = 0.25*epsilon/(s*s*s);//d2f/drho2
        f[0]   = 0;
        df[0]  = 0;
        d2f[0] = 0;
    }
}
