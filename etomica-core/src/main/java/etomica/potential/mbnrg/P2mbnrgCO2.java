package etomica.potential.mbnrg;

import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecular;
import etomica.space.Vector;


public class P2mbnrgCO2 implements IPotentialMolecular {
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
    public P2mbnrgCO2(String mon1, String mon2) {
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

    public double energy(IMoleculeList molecules) {
        
        IAtomList atoms1 = molecules.get(0).getChildList();
        Vector vC_1 = atoms1.get(0).getPosition();
        Vector vO1_1 = atoms1.get(1).getPosition();
        Vector vO2_1 = atoms1.get(2).getPosition();
        IAtomList atoms2 = molecules.get(1).getChildList();
        Vector vC_2 = atoms2.get(0).getPosition();
        Vector vO1_2 = atoms2.get(1).getPosition();
        Vector vO2_2 = atoms2.get(2).getPosition();

//         double[] d12 = new double[]{mon1[0] - mon2[0], mon1[1] - mon2[1], mon1[2] - mon2[2]};
//
//         double r12sq = d12[0] * d12[0] + d12[1] * d12[1] + d12[2] * d12[2];
//         double r12 = Math.sqrt(r12sq);
        //first atom only
         double r12 = Math.sqrt(vC_1.Mv1Squared(vC_2));

        double[] v = new double[15];

        double sw;
        
        v[0] = Definitions.v_exp(m_d_intra_AB, m_k_intra_AB, vC_1, vO1_1);
        v[1] = Definitions.v_exp(m_d_intra_AB, m_k_intra_AB, vC_1, vO2_1);
        v[2] = Definitions.v_exp(m_d_intra_BB, m_k_intra_BB, vO1_1, vO2_1);

        v[3] = Definitions.v_exp(m_d_intra_AB, m_k_intra_AB, vC_2, vO1_2);
        v[4] = Definitions.v_exp(m_d_intra_AB, m_k_intra_AB, vC_2, vO2_2);
        v[5] = Definitions.v_exp(m_d_intra_BB, m_k_intra_BB, vO1_2, vO2_2);

        v[6] = Definitions.v_exp(m_d_AA, m_k_AA, vC_1, vC_2);
        v[7] = Definitions.v_exp(m_d_AB, m_k_AB, vC_1, vO1_2);
        v[8] = Definitions.v_exp(m_d_AB, m_k_AB, vC_1, vO2_2);

        v[9] = Definitions.v_exp(m_d_AB, m_k_AB, vO1_1, vC_2);
        v[10] = Definitions.v_exp(m_d_BB, m_k_BB, vO1_1, vO1_2);
        v[11] = Definitions.v_exp(m_d_BB, m_k_BB, vO1_1, vO2_2);

        v[12] = Definitions.v_exp(m_d_AB, m_k_AB, vO2_1, vC_2);
        v[13] = Definitions.v_exp(m_d_BB, m_k_BB, vO2_1, vO1_2);
        v[14] = Definitions.v_exp(m_d_BB, m_k_BB, vO2_1, vO2_2);

        sw = f_switch(r12);

        return sw * Poly_A1B2_2b.poly_eval(coefficients, v);
    }










}

