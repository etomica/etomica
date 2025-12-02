package etomica.potential.mbnrg;

import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.space.Vector;

public class P3mbnrgCO2 implements IPotentialMolecular {
    double[] coefficients;
    String mon1;
    String mon2;
    String mon3;
    double m_k_x_inter_A_A_0;
    double m_k_x_inter_A_B_0;
    double m_k_x_inter_B_B_0;
    // Inner cutoff
    double m_ri;
    // Outer cutoff
    double m_ro;

    public P3mbnrgCO2(String mon1, String mon2, String mon3) {
        this.mon1 = mon1;
        this.mon2 = mon2;
        this.mon3 = mon3;
        //initialize PIP coefficients and parameters
        A1B2_init_3b.init(this);

    }
    double f_switch(double r) {
        if (r > m_ro) {
//            g = 0.0;
            return 0.0;
        } else if (r > m_ri) {double t1 = Math.PI / (m_ro - m_ri);
        double x = (r - m_ri) * t1;
//            g = -std::sin(x) * t1 / 2.0;
            return (1.0 + Math.cos(x)) / 2.0;
        } else {
//            g = 0.0;
            return 1.0;
        }
    }

    public double energy(IMoleculeList molecules) {
            double sw;

            IAtomList atoms1 = molecules.get(0).getChildList();
            Vector vC_1 = atoms1.get(0).getPosition();
            Vector vO1_1 = atoms1.get(1).getPosition();
            Vector vO2_1 = atoms1.get(2).getPosition();

            IAtomList atoms2 = molecules.get(1).getChildList();
            Vector vC_2 = atoms2.get(0).getPosition();
            Vector vO1_2 = atoms2.get(1).getPosition();
            Vector vO2_2 = atoms2.get(2).getPosition();

            IAtomList atoms3 = molecules.get(2).getChildList();
            Vector vC_3 = atoms3.get(0).getPosition();
            Vector vO1_3 = atoms3.get(1).getPosition();
            Vector vO2_3 = atoms3.get(2).getPosition();

            double d12r = Math.sqrt(vC_1.Mv1Squared(vC_2));

            double d13r = Math.sqrt(vC_1.Mv1Squared(vC_3));

            double d23r = Math.sqrt(vC_2.Mv1Squared(vC_3));

            double[] xs = new double[27];

            xs[0] = Definitions.v_exp(m_k_x_inter_A_A_0, vC_1, vC_2);
            xs[1] = Definitions.v_exp(m_k_x_inter_A_A_0, vC_1, vC_3);
            xs[2] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_1, vO1_2);
            xs[3] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_1, vO2_2);
            xs[4] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_1, vO1_3);
            xs[5] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_1, vO2_3);
            xs[6] = Definitions.v_exp(m_k_x_inter_A_A_0, vC_2, vC_3);
            xs[7] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_2, vO1_1);
            xs[8] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_2, vO2_1);
            xs[9] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_2, vO1_3);
            xs[10] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_2, vO2_3);
            xs[11] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_3, vO1_1);
            xs[12] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_3, vO2_1);
            xs[13] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_3, vO1_2);
            xs[14] = Definitions.v_exp(m_k_x_inter_A_B_0, vC_3, vO2_2);
            xs[15] = Definitions.v_exp(m_k_x_inter_B_B_0, vO1_1, vO1_2);
            xs[16] = Definitions.v_exp(m_k_x_inter_B_B_0, vO1_1, vO2_2);
            xs[17] = Definitions.v_exp(m_k_x_inter_B_B_0, vO1_1, vO1_3);
            xs[18] = Definitions.v_exp(m_k_x_inter_B_B_0, vO1_1, vO2_3);
            xs[19] = Definitions.v_exp(m_k_x_inter_B_B_0, vO2_1, vO1_2);
            xs[20] = Definitions.v_exp(m_k_x_inter_B_B_0, vO2_1, vO2_2);
            xs[21] = Definitions.v_exp(m_k_x_inter_B_B_0, vO2_1, vO1_3);
            xs[22] = Definitions.v_exp(m_k_x_inter_B_B_0, vO2_1, vO2_3);
            xs[23] = Definitions.v_exp(m_k_x_inter_B_B_0, vO1_2, vO1_3);
            xs[24] = Definitions.v_exp(m_k_x_inter_B_B_0, vO1_2, vO2_3);
            xs[25] = Definitions.v_exp(m_k_x_inter_B_B_0, vO2_2, vO1_3);
            xs[26] = Definitions.v_exp(m_k_x_inter_B_B_0, vO2_2, vO2_3);

            double sw12 = f_switch(d12r);
            double sw13 = f_switch(d13r);
            double sw23 = f_switch(d23r);

            sw = sw12 * sw13 + sw12 * sw23 + sw13 * sw23 - 2 * sw12 * sw23 * sw13;

            return sw * Poly_A1B2_3b.poly_eval(xs, coefficients);

    }
    
}
