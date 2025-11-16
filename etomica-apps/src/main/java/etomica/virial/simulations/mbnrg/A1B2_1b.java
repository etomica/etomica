package etomica.virial.simulations.mbnrg;

import java.util.ArrayList;

public class A1B2_1b {
    String mon1;
    Double[] coefficients;
    double d_intra_AB;
    double d_intra_BB;
    double k_intra_AB;
    double k_intra_BB;

    public A1B2_1b(String mon1){
        this.mon1 = mon1;
        A1B2_init_1b.init(this);
    }
    ArrayList<Double> eval(Double[] xyz, int nmon) {
        ArrayList<Double> energies = new ArrayList<>(nmon);
        for (int i = 0; i < nmon; i++) {
            energies.add(0.0);
        }

        for (int j = 0; j < nmon; j++) {

//            std::copy(xyz + j * 9, xyz + (j + 1) * 9, xcrd);
//        const double* A0 = xcrd + 0;
//        const double* B0 = xcrd + 3;
//        const double* B1 = xcrd + 6;

            int srcPos = j * 9;
            double[] A0 = new double[3];
            double[] B0 = new double[3];
            double[] B1 = new double[3];

            System.arraycopy(xyz, srcPos, A0, 0, 3);
            System.arraycopy(xyz, srcPos + 3, B0, 0, 3);
            System.arraycopy(xyz, srcPos + 6, B1, 0, 3);

            double[] x = new double[3];
            double[] g = new double[3];
            x[0] = Definitions.v_intra(k_intra_AB, d_intra_AB, A0, B0);
            x[1] = Definitions.v_intra(k_intra_AB, d_intra_AB, A0, B1);
            x[2] = Definitions.v_intra(k_intra_BB, d_intra_BB, B0, B1);
            energies.set(j, poly_eval(coefficients, x, g));
        }
        return energies;
    }
//    static ArrayList<Double> eval(ArrayList<Double> xyz, ArrayList<Double> grad, Integer nmon, ArrayList<Double> virial) {
//
//        ArrayList<Double> energies = new ArrayList<>(nmon);
//        for (int i = 0; i < nmon; i++) {
//            energies.add(0.0);
//        }
//
//        for (int j = 0; j < nmon; j++) {
////            double[] xcrd = new double[9];
////            std::copy(xyz + j * 9, xyz + (j + 1) * 9, xcrd);
////            Double A0 = xcrd + 0;
////            Double B0 = xcrd + 3;
////            Double B1 = xcrd + 6;
//            int base = j * 9;
//            double[] A0 = { xyz.get(base), xyz.get(base+1), xyz.get(base+2) };
//            double[] B0 = { xyz.get(base+3), xyz.get(base+4), xyz.get(base+5) };
//            double[] B1 = { xyz.get(base+6), xyz.get(base+7), xyz.get(base+8) };
//
//            double[] x = new double[3];
//            double[] g = new double[3];
//            //Convert coordinates to non-linear form
//            x[0] = Definitions.v_intra(k_intra_AB, d_intra_AB, A0, B0);
//            x[1] = Definitions.v_intra(k_intra_AB, d_intra_AB, A0, B1);
//            x[2] = Definitions.v_intra(k_intra_BB, d_intra_BB, B0, B1);
//            energies.set(j, poly_eval(coefficients, x, g));
//            double[] gA0 = {0.0, 0.0, 0.0};
//            double[] gB0 = {0.0, 0.0, 0.0};
//            double[] gB1 = {0.0, 0.0, 0.0};
//            Definitions.g_intra(g[0], k_intra_AB, d_intra_AB, A0, B0, gA0, gB0);
//            Definitions.g_intra(g[1], k_intra_AB, d_intra_AB, A0, B1, gA0, gB1);
//            Definitions.g_intra(g[2], k_intra_BB, d_intra_BB, B0, B1, gB0, gB1);
//            for (int i = 0; i < 9; i++) grad.set(i + j * 9, xgrd[i]);
//
//            if (virial != null) {
//                virial.set(0, virial.get(0)
//                        + (-A0[0] * gA0[0] - B0[0] * gB0[0] - B1[0] * gB1[0]));
//
//                virial.set(1, virial.get(1)
//                        + (-A0[0] * gA0[1] - B0[0] * gB0[1] - B1[0] * gB1[1]));
//
//                virial.set(2, virial.get(2)
//                        + (-A0[0] * gA0[2] - B0[0] * gB0[2] - B1[0] * gB1[2]));
//
//                virial.set(4, virial.get(4)
//                        + (-A0[1] * gA0[1] - B0[1] * gB0[1] - B1[1] * gB1[1]));
//
//                virial.set(5, virial.get(5)
//                        + (-A0[1] * gA0[2] - B0[1] * gB0[2] - B1[1] * gB1[2]));
//
//                virial.set(8, virial.get(8)
//                        + (-A0[2] * gA0[2] - B0[2] * gB0[2] - B1[2] * gB1[2]));
//
//                virial.set(3, virial.get(1));
//                virial.set(6, virial.get(2));
//                virial.set(7, virial.get(5));
//            }
//        }
//
//
//        return energies;
//    }
    static double poly_eval(Double[] a, double[] x, double[] g) {
        double t1 = a[1];
        double t2 = a[4];
        double t6 = x[2];
        double t4 = a[20] * t6;
        double t5 = a[11];
        double t7 = (t4 + t5) * t6;
        double t9 = (t2 + t7) * t6;
        double t12 = a[0];
        double t13 = a[3];
        double t15 = a[14] * t6;
        double t16 = a[8];
        double t18 = (t15 + t16) * t6;
        double t20 = (t13 + t18) * t6;
        double t21 = a[5];
        double t23 = a[18] * t6;
        double t24 = a[7];
        double t26 = (t23 + t24) * t6;
        double t27 = a[13];
        double t25 = x[1];
        double t28 = t27 * t25;
        double t29 = a[17];
        double t30 = t29 * t6;
        double t31 = a[6];
        double t33 = (t28 + t30 + t31) * t25;
        double t35 = (t21 + t26 + t33) * t25;
        double t38 = a[2];
        double t40 = a[16] * t6;
        double t41 = a[9];
        double t43 = (t40 + t41) * t6;
        double t44 = a[15];
        double t45 = t44 * t25;
        double t46 = a[12];
        double t47 = t46 * t6;
        double t48 = a[10];
        double t50 = (t45 + t47 + t48) * t25;
        double t52 = (t38 + t43 + t50) * t25;
        double t54 = a[19] * t25;
        double t56 = (t54 + t47 + t48) * t25;
        double t53 = x[0];
        double t57 = t27 * t53;
        double t59 = (t57 + t45 + t30 + t31) * t53;
        double t61 = (t21 + t26 + t56 + t59) * t53;
        double t93 = (2.0 * t15 + t16) * t6;
        double t95 = 2.0 * t23;
        double t100 = t46 * t25;
        g[0] = ((2.0 * t57 + t45 + t30 + t31) * t53 + t21 + t26 + t56 + t59) * t53 + t12 + t20 + t52 + t61;
        g[1] = ((2.0 * t28 + t30 + t31) * t25 + t21 + t26 + t33) * t25 + t12 + t20 + t35 +
                ((2.0 * t45 + t47 + t48) * t25 + t38 + t43 + t50 + (t44 * t53 + t47 + t48 + 2.0 * t54) * t53) * t53;
        g[2] = ((2.0 * t4 + t5) * t6 + t2 + t7) * t6 + t1 + t9 + (t93 + t13 + t18 + (t29 * t25 + t24 + t95) * t25) * t25 +
                (t93 + t13 + t18 + (t100 + 2.0 * t40 + t41) * t25 + (t29 * t53 + t100 + t24 + t95) * t53) * t53;
        double e = (t1 + t9) * t6 + (t12 + t20 + t35) * t25 + (t12 + t20 + t52 + t61) * t53;


        return e;
    }


}
