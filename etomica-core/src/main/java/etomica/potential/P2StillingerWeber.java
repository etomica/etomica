/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.exception.MethodNotImplementedException;
import etomica.units.Erg;

/**
 * 2-body contribution to Stillinger-Weber potential.
 *
 * <a href="https://doi.org/10.1103/PhysRevB.31.5262">DOI:10.1103/PhysRevB.31.5262</a>
 */
public class P2StillingerWeber implements IPotential2 {

    public double sigma, epsilon, A, B, a, sigma2;
    public int p, q;

    public static IPotential2 makeTruncated(double sigma, double epsilon, double A, double B, int p, int q, double a, TruncationFactory tf) {
        return tf.make(new P2StillingerWeber(sigma, epsilon, A, B, p, q, a));
    }

    public P2StillingerWeber(double sigma, double epsilon, double A, double B, int p, int q, double a) {
        super();
        setParameters(sigma, epsilon, A, B, p, q, a);
    }

    public void setParameters(double sigma, double epsilon, double A, double B, int p, int q, double a) {
        this.sigma = sigma;
        this.epsilon = epsilon;
        this.A = A;
        this.B = B;
        this.p = p;
        this.q = q;
        this.a = a;
        sigma2 = sigma*sigma;
    }

    protected double intpow(double b2, double b, int e) {
        switch (e) {
            case 0:
                return 1;
            case 1:
                return b;
            case 2:
                return b2;
            case 3:
                return b2 * b;
            case 4:
                return b2 * b2;
            case 5:
                return b2 * b2 * b;
            case 6:
                return b2 * b2 * b2;
            case 7:
                return b2 * b2 * b2 * b;
            default:
                throw new RuntimeException("Implement high powers!");
        }
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        double sr2 = sigma2/r2;
        double sr = Math.sqrt(sr2);
        double rs = 1/sr;
        if (rs > a) return 0;
        return epsilon*A*(B*intpow(sr2, sr, p) - intpow(sr2, sr, q) * Math.exp(1/(rs - a)));
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double sr2 = sigma2/r2;
        double sr = Math.sqrt(sr2);
        double rs = 1/sr;
        if (rs > a) return 0;
        double srp = intpow(sr2, sr, p), srq = intpow(sr2, sr, q);
        double exp = Math.exp(1/(rs - a));
        return -epsilon*A*exp*(B*p*srp - q*srq + (B*srp - srq) * rs / ((rs-a)*(rs-a)));
    }

    public void u012add(double r2, double[] u012) {
        double sr2 = sigma2/r2;
        double sr = Math.sqrt(sr2);
        double rs = 1/sr;
        if (rs > a) return;
        double srp = intpow(sr2, sr, p), srq = intpow(sr2, sr, q);
        double exp = Math.exp(1/(rs - a));
        u012[0] += epsilon*A*(B*srp - srq) * exp;
        u012[1] += -epsilon*A*exp*(B*p*srp - q*srq + (B*srp - srq) * rs / ((rs-a)*(rs-a)));
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        throw new MethodNotImplementedException("nope");
    }

    public static void main(String[] args) {
        P2StillingerWeber p2 = new P2StillingerWeber(2.0951, Erg.UNIT.toSim(3.4723e-12), 7.049556277, 0.6022245584, 4, 0, 1.80);
        for (int i=0; i<100; i++) {
            double r = 2.0951*0.80 + i*0.02;
            double[] u012 = new double[3];
            p2.u012add(r*r, u012);
            System.out.println(r+" "+u012[0]+" "+u012[1]);
        }
    }
}
