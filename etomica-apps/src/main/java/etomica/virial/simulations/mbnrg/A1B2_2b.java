package etomica.virial.simulations.mbnrg;

public class A1B2_2b {
    Double[] coefficients;
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
    // Constructor, initialization of coefficients
    public A1B2_2b(String mon1, String mon2) {
        this.mon1 = mon1;
        this.mon2 = mon2;
        //initialize PIP coefficients and parameters
        A1B2_init_2b.init(this);

    }



}

