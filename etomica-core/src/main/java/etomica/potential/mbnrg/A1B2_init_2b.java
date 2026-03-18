package etomica.potential.mbnrg;

import java.util.Arrays;

public class A1B2_init_2b {
//    public static void main(String[] args) {
//        double[] x = Definitions.readCoefficients("co2_archive_2b");
//        System.out.println("Coefficients: " + Arrays.toString(x));
//
//    }
    public static void init(P2mbnrgCO2 t){

        if (t.mon1.equals("co2_archive") && t.mon2.equals("co2_archive")) {
            t.coefficients = Definitions.readCoefficients("co2_archive_2b");
            t.m_d_intra_AB = 2.533470783111825e+00;  // A^(-1))
            t.m_k_intra_AB = 1.772445777810116e+00;  // A^(-1))
            t.m_d_intra_BB = 1.299246633223584e+00;  // A^(-1))
            t.m_k_intra_BB = 1.546435354569379e+00;  // A^(-1))
            t.m_d_AA = 4.042621705776338e+00;        // A^(-1))
            t.m_k_AA = 1.375890559471216e+00;        // A^(-1))
            t.m_d_AB = 4.893559114105628e+00;        // A^(-1))
            t.m_k_AB = 7.248440685272545e-01;        // A^(-1))
            t.m_d_BB = 2.785825083032128e+00;        // A^(-1))
            t.m_k_BB = 2.691866638221668e+00;        // A^(-1))
            t.m_r2i = 8.000000000000000e+00;         // A
            t.m_r2f = 9.000000000000000e+00;         // A
        }

        if (t.mon1.equals("co2") && t.mon2.equals("co2")) {
            t.coefficients = Definitions.readCoefficients("co2_co2_2b");
            t.m_d_intra_AB = 3.098180719521002e+00;  // A^(-1))
            t.m_k_intra_AB = 1.267379326287500e+00;  // A^(-1))
            t.m_d_intra_BB = 2.152649569148020e+00;  // A^(-1))
            t.m_k_intra_BB = 6.904623477544016e-01;  // A^(-1))
            t.m_d_AA = 4.195717132005157e+00;        // A^(-1))
            t.m_k_AA = 6.959794865790364e-02;        // A^(-1))
            t.m_d_AB = 6.633057823695092e+00;        // A^(-1))
            t.m_k_AB = 9.311550654726656e-01;        // A^(-1))
            t.m_d_BB = 5.614642310195038e+00;        // A^(-1))
            t.m_k_BB = 1.053234040860423e+00;        // A^(-1))
            t.m_r2i = 7.000000000000000e+00;         // A
            t.m_r2f = 9.000000000000000e+00;         // A
        }

        if (t.mon1.equals("co2cm5100") && t.mon2.equals("co2cm5100")) {
            t.coefficients = Definitions.readCoefficients("co2cm5100_2b");
            t.m_d_intra_AB = 3.098180719521002e+00;  // A^(-1))
            t.m_k_intra_AB = 1.267379326287500e+00;  // A^(-1))
            t.m_d_intra_BB = 2.152649569148020e+00;  // A^(-1))
            t.m_k_intra_BB = 6.904623477544016e-01;  // A^(-1))
            t.m_d_AA = 4.195717132005157e+00;        // A^(-1))
            t.m_k_AA = 6.959794865790364e-02;        // A^(-1))
            t.m_d_AB = 6.633057823695092e+00;        // A^(-1))
            t.m_k_AB = 9.311550654726656e-01;        // A^(-1))
            t.m_d_BB = 5.614642310195038e+00;        // A^(-1))
            t.m_k_BB = 1.053234040860423e+00;        // A^(-1))
            t.m_r2i = 7.000000000000000e+00;         // A
            t.m_r2f = 9.000000000000000e+00;         // A
        }

        if (t.mon1.equals("co2cm595") && t.mon2.equals("co2cm595")) {
            t.coefficients = Definitions.readCoefficients("co2cm595_2b");

            t.m_d_intra_AB = 4.903036747511356e+00;  // A^(-1))
            t.m_k_intra_AB = 9.019913232531600e-01;  // A^(-1))
            t.m_d_intra_BB = 2.703937855283358e+00;  // A^(-1))
            t.m_k_intra_BB = 1.004834209745757e+00;  // A^(-1))
            t.m_d_AA = 4.537101521312326e+00;        // A^(-1))
            t.m_k_AA = 9.163852883742113e-02;        // A^(-1))
            t.m_d_AB = 6.300498483648834e+00;        // A^(-1))
            t.m_k_AB = 9.021417601624032e-01;        // A^(-1))
            t.m_d_BB = 5.827456239776177e+00;        // A^(-1))
            t.m_k_BB = 1.089701616045492e+00;        // A^(-1))
            t.m_r2i = 7.000000000000000e+00;         // A
            t.m_r2f = 9.000000000000000e+00;         // A
        }

        if (t.mon1.equals("co2cm590") && t.mon2.equals("co2cm590")) {
            t.coefficients = Definitions.readCoefficients("co2cm590_2b");

            t.m_d_intra_AB = 2.329138487337566e+00;  // A^(-1))
            t.m_k_intra_AB = 1.434391663074755e+00;  // A^(-1))
            t.m_d_intra_BB = 3.229816367099317e+00;  // A^(-1))
            t.m_k_intra_BB = 9.702035861774958e-01;  // A^(-1))
            t.m_d_AA = 4.006166160264772e+00;        // A^(-1))
            t.m_k_AA = 1.117857009888089e-01;        // A^(-1))
            t.m_d_AB = 5.789451971684920e+00;        // A^(-1))
            t.m_k_AB = 8.757387756335298e-01;        // A^(-1))
            t.m_d_BB = 5.798891302817732e+00;        // A^(-1))
            t.m_k_BB = 1.133531479542603e+00;        // A^(-1))
            t.m_r2i = 7.000000000000000e+00;         // A
            t.m_r2f = 9.000000000000000e+00;         // A
        }

        if (t.mon1.equals("co2cm585") && t.mon2.equals("co2cm585")) {
            t.coefficients = Definitions.readCoefficients("co2cm585_2b");

            t.m_d_intra_AB = 2.165524277685071e+00;  // A^(-1))
            t.m_k_intra_AB = 2.061892967817994e+00;  // A^(-1))
            t.m_d_intra_BB = 2.835603377312035e+00;  // A^(-1))
            t.m_k_intra_BB = 3.459793262493814e+00;  // A^(-1))
            t.m_d_AA = 4.235604053709159e+00;        // A^(-1))
            t.m_k_AA = 1.141488887995241e-01;        // A^(-1))
            t.m_d_AB = 4.667583436073858e+00;        // A^(-1))
            t.m_k_AB = 8.691887001880019e-01;        // A^(-1))
            t.m_d_BB = 5.473818228420447e+00;        // A^(-1))
            t.m_k_BB = 1.157204319628839e+00;        // A^(-1))
            t.m_r2i = 7.000000000000000e+00;         // A
            t.m_r2f = 9.000000000000000e+00;         // A
        }

        if (t.mon1.equals("co2cm5875") && t.mon2.equals("co2cm5875")) {
            t.coefficients = Definitions.readCoefficients("co2cm5875_2b");

            t.m_d_intra_AB = 3.134332315833306e+00;  // A^(-1))
            t.m_k_intra_AB = 1.282904902772805e+00;  // A^(-1))
            t.m_d_intra_BB = 4.434199799497058e+00;  // A^(-1))
            t.m_k_intra_BB = 8.407363065419096e-01;  // A^(-1))
            t.m_d_AA = 1.800437262469989e+00;        // A^(-1))
            t.m_k_AA = 1.216397961826567e-01;        // A^(-1))
            t.m_d_AB = 6.926666442517130e+00;        // A^(-1))
            t.m_k_AB = 8.810761305307210e-01;        // A^(-1))
            t.m_d_BB = 4.127797348559158e+00;        // A^(-1))
            t.m_k_BB = 1.148765431751049e+00;        // A^(-1))
            t.m_r2i = 7.000000000000000e+00;         // A
            t.m_r2f = 9.000000000000000e+00;         // A
        }

        if (t.mon1.equals("co2cm580") && t.mon2.equals("co2cm580")) {
            t.coefficients = Definitions.readCoefficients("co2cm580_2b");


            t.m_d_intra_AB = 4.158829702809511e+00;  // A^(-1))
            t.m_k_intra_AB = 8.744952741703924e-01;  // A^(-1))
            t.m_d_intra_BB = 9.999985981119849e+00;  // A^(-1))
            t.m_k_intra_BB = 4.486543492998215e-01;  // A^(-1))
            t.m_d_AA = 3.899891165573160e+00;        // A^(-1))
            t.m_k_AA = 1.074038746629509e+00;        // A^(-1))
            t.m_d_AB = 4.292889407717251e+00;        // A^(-1))
            t.m_k_AB = 6.400679425638476e-01;        // A^(-1))
            t.m_d_BB = 4.386671648943636e+00;        // A^(-1))
            t.m_k_BB = 1.326376079648352e+00;        // A^(-1))
            t.m_r2i = 7.000000000000000e+00;         // A
            t.m_r2f = 9.000000000000000e+00;         // A
        }


    }
    }
