package etomica.potential.mbnrg;

import etomica.space.Vector;

public class Definitions {
    public static final double EMAX1B = 100000.0;
    public static final double EPSILON = 1E-50;
    static double v_intra(double k, double r0, Vector v1, Vector v2) {
        double d = Math.sqrt(v1.Mv1Squared(v2));

        return Math.exp(-k * (d - r0));
    }
    static void g_intra(double g, double k, double r0, double[] a1, double[] a2, double[] g1,
                        double[] g2) {
        double[] dx = {a1[0] - a2[0], a1[1] - a2[1], a1[2] - a2[2]};
        double dsq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
        double d = Math.sqrt(dsq);

        double gg = -k * g * Math.exp(-k * (d - r0)) / d;

        for (int i = 0; i < 3; ++i) {
            g1[i] += gg * dx[i];
            g2[i] -= gg * dx[i];
        }
    }
    static double v_exp(double r0, double k, Vector v1, Vector v2){

        double r = Math.sqrt(v1.Mv1Squared(v2));

        return Math.exp(k * (r0 - r));
//            double gg = -k * exp1 / r;
//
//            g[0] *= gg;
//            g[1] *= gg;
//            g[2] *= gg;


    };
    static double v_exp(double k, Vector v1, Vector v2){

        double r = Math.sqrt(v1.Mv1Squared(v2));

        return Math.exp(k * (- r));
//            double gg = -k * exp1 / r;
//
//            g[0] *= gg;
//            g[1] *= gg;
//            g[2] *= gg;


    };

}
