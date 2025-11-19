package etomica.potential.mbnrg;

import java.util.ArrayList;

public class A1B2_2b {
    double[] coefficients;
    String mon1;
    String mon2;
    double m_d_intra_AB;
    double m_k_intra_AB;
    double m_d_intra_BB;
    double m_k_intra_BB;
    double m_d_AA;
    double m_k_AA;
    double m_d_AB;
    double m_k_AB;
    double m_d_BB;
    double m_k_BB;
    // Inner cutoff
    double m_r2i;
    // Outer cutoff
    double m_r2f;
    // constructor, initialization of coefficients
    public A1B2_2b(String mon1, String mon2) {
        this.mon1 = mon1;
        this.mon2 = mon2;
        //initialize PIP coefficients and parameters
        A1B2_init_2b.init(this);

    }
    double f_switch(double r){
//        double g;
        if (r > m_r2f) {
//            g = 0.0;
            return 0.0;
        } else if (r > m_r2i) {
            double t1 = Math.PI / (m_r2f - m_r2i);
            double x = (r - m_r2i) * t1;
//            g = -Math.sin(x) * t1 / 2.0;
            return (1.0 + Math.cos(x)) / 2.0;
        } else {
//            g = 0.0;
            return 1.0;
        }
    }
    class variable {
        double v_exp(double r0, double k, double[] p1, double[] p2){
            g[0] = p1[0] - p2[0];
            g[1] = p1[1] - p2[1];
            g[2] = p1[2] - p2[2];

            double r = Math.sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);

            double exp1 = Math.exp(k * (r0 - r));
            double gg = -k * exp1 / r;

            g[0] *= gg;
            g[1] *= gg;
            g[2] *= gg;

            return exp1;

        };

        double[] g = new double[3];  // diff(value, p1 - p2)
    };
    //return sum of n-dimer energies (2-body PIP, repulsive)
    Double eval(Double[] xyz1, Double[] xyz2, Integer ndim) {

        ArrayList<Double> energies = new ArrayList<>(ndim);
        for (int i = 0; i < ndim; i++) {
            energies.add(0.0);
        }

        for (int j = 0; j < ndim; j++) {
//            double mon1[9];
//            double mon2[9];
//
//            std::copy(xyz1 + j * 9, xyz1 + (j + 1) * 9, mon1);
//            std::copy(xyz2 + j * 9, xyz2 + (j + 1) * 9, mon2);

            double[] mon1 = new double[9];
            double[] mon2 = new double[9];

            int srcPos = j * 9;
            System.arraycopy(xyz1, srcPos, mon1, 0, 9);
            System.arraycopy(xyz2, srcPos, mon2, 0, 9);
            // ##DEFINE HERE## right now it assumes 1st atom of each monomer
             double[] d12 = new double[]{mon1[0] - mon2[0], mon1[1] - mon2[1], mon1[2] - mon2[2]};

             double r12sq = d12[0] * d12[0] + d12[1] * d12[1] + d12[2] * d12[2];
             double r12 = Math.sqrt(r12sq);

            if (r12 > m_r2f) continue;

//            double[] xcrd = new double[18];  // coordinates of real sites ONLY

//            std::copy(mon1, mon1 + 9, xcrd);
//            std::copy(mon2, mon2 + 9, xcrd + 9);
//
//
//             double* A_1_a = xcrd + 0;
//             double* B_1_a = xcrd + 3;
//             double* B_2_a = xcrd + 6;
//
//             double* A_1_b = xcrd + 9;
//             double* B_1_b = xcrd + 12;
//             double* B_2_b = xcrd + 15;
            double[] A_1_a = new double[3];
            double[] B_1_a = new double[3];
            double[] B_2_a = new double[3];
            double[] A_1_b = new double[3];
            double[] B_1_b = new double[3];
            double[] B_2_b = new double[3];

            System.arraycopy(mon1, 0, A_1_a, 0, 3);
            System.arraycopy(mon1, 3, B_1_a, 0, 3);
            System.arraycopy(mon1, 6, B_2_a, 0, 3);
            System.arraycopy(mon2, 0, A_1_b, 0, 3);
            System.arraycopy(mon2, 3, B_1_b, 0, 3);
            System.arraycopy(mon2, 6, B_2_b, 0, 3);

            double[] v = new double[15];

            double sw = 0.0;
//            double gsw = 0.0;

            variable[] vr = new variable[15];

            v[0] = vr[0].v_exp(m_d_intra_AB, m_k_intra_AB, A_1_a, B_1_a);
            v[1] = vr[1].v_exp(m_d_intra_AB, m_k_intra_AB, A_1_a, B_2_a);
            v[2] = vr[2].v_exp(m_d_intra_BB, m_k_intra_BB, B_1_a, B_2_a);

            v[3] = vr[3].v_exp(m_d_intra_AB, m_k_intra_AB, A_1_b, B_1_b);
            v[4] = vr[4].v_exp(m_d_intra_AB, m_k_intra_AB, A_1_b, B_2_b);
            v[5] = vr[5].v_exp(m_d_intra_BB, m_k_intra_BB, B_1_b, B_2_b);

            v[6] = vr[6].v_exp(m_d_AA, m_k_AA, A_1_a, A_1_b);
            v[7] = vr[7].v_exp(m_d_AB, m_k_AB, A_1_a, B_1_b);
            v[8] = vr[8].v_exp(m_d_AB, m_k_AB, A_1_a, B_2_b);

            v[9] = vr[9].v_exp(m_d_AB, m_k_AB, B_1_a, A_1_b);
            v[10] = vr[10].v_exp(m_d_BB, m_k_BB, B_1_a, B_1_b);
            v[11] = vr[11].v_exp(m_d_BB, m_k_BB, B_1_a, B_2_b);

            v[12] = vr[12].v_exp(m_d_AB, m_k_AB, B_2_a, A_1_b);
            v[13] = vr[13].v_exp(m_d_BB, m_k_BB, B_2_a, B_1_b);
            v[14] = vr[14].v_exp(m_d_BB, m_k_BB, B_2_a, B_2_b);

            sw = f_switch(r12);

            energies.set(j, sw * Poly_A1B2_2b.poly_eval(coefficients, v));
        }

        double energy = 0.0;
        for (int i = 0; i < ndim; i++) {
            energy += energies.get(i);
        }


        return energy;
    }







}

