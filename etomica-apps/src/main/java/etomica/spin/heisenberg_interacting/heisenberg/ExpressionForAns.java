package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialAtomic;
import etomica.space.Space;
import etomica.space.Vector;

import static etomica.math.SpecialFunctions.besselI;

/**
 * Anx expressions for computing v_E and v_EE in the mapping.
 *
 * @author Weisong Lin
 */

public class ExpressionForAns {
    protected Vector ei, ej;
    protected final double mu, J, bt, bJ, bmu; //TODO should I add final here
    protected double Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    protected double d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
    protected double Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
    protected double d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;

    public ExpressionForAns(Space space, double dipoleMagnitude, double interactionS, double beta) {
        ei = space.makeVector();
        ej = space.makeVector();
        J = interactionS;
        mu = dipoleMagnitude;
        bt = beta;
        bJ = bt * J;
        bmu = bt * mu;
    }

    public void doCalculation(IAtomList atoms, int nMax) {
        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());
        double t1 = Math.acos(ei.getX(0));
//        double t2 = Math.acos(ej.getX(0));
        double I0bJ = besselI(0, bJ);
        double I1bJ = besselI(1, bJ);
        double I2bJ = besselI(2, bJ);

        double test0 = 0;
        double test1 = 0;
        double test2 = 0;
        double test3 = 0;
        double test3_1 = 0;
        double test3_2 = 0;
        double test4 = 0;


        Axc0 = bmu * (I0bJ + I1bJ) * (-1 + Math.cos(t1));
        dAxc0 = -bmu * (I0bJ + I1bJ) * Math.sin(t1);
        Axs0 = 0;
        dAxs0 = 0;
        Axc1 = -0.25 * bmu * bmu * Math.sin(t1) * Math.sin(t1) * (I0bJ + 2 * I1bJ + I2bJ);
        dAxc1 = -0.25 * bmu * bmu * (I0bJ + 2 * I1bJ + I2bJ) * Math.sin(t1) * Math.sin(t1);
        Axs1 = 0;
        dAxs1 = 0;

        d2Axc0 = -bmu * Math.cos(t1) * (I0bJ + I1bJ);
        d2Axs0 = 0;
        d3Axc0 = bmu * Math.sin(t1) * (I0bJ + I1bJ);
        d3Axs0 = 0;

        test0 = d3Axs0;

        //x direction
        if (nMax > 0) {
            int n = nMax;
            int n2 = n * n;
            int n3 = n2 * n;
            int n4 = n2 * n2;
            double InbJ = besselI(n, bJ);
            double Inm1bJ = besselI(n - 1, bJ);
            double Inm2bJ = besselI(n - 2, bJ);
            double Inp1bJ = besselI(n + 1, bJ);
            Axc0 = 2 * bmu * (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
                    + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * Math.sin(t1) * Math.sin(n * t1))
                    / (bJ + 4 * bJ * n4);

            Axs0 = (-2 * bmu * n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * Math.cos(t1) * Math.cos(n * t1)
                    - 2 * bmu * (bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * Math.sin(t1) * Math.sin(n * t1)
            ) / (bJ + 4 * bJ * n4);

            dAxc0 = (2 * InbJ * (-Math.cos(n * t1) * Math.sin(t1) + (n - 2 * n3) * Math.cos(t1) * Math.sin(n * t1))
                    + Inm1bJ * (1 + n - 2 * n3) * Math.sin((n - 1) * t1)
                    + (-1 + n - 2 * n3) * Inp1bJ * Math.sin((n + 1) * t1)
            ) * bmu / (1 + 4 * n4);

            dAxs0 = ((1 + 2 * n + 2 * n2) * Inm1bJ * Math.sin((n - 1) * t1)
                    + 2 * InbJ * (-2 * n * Math.cos(n * t1) * Math.sin(t1) + (1 + 2 * n2) * Math.cos(t1) * Math.sin(n * t1))
                    + (1 - 2 * n + 2 * n2) * Inp1bJ * Math.sin((n + 1) * t1)
            ) * bmu / (1 + 4 * n4);

            Axc1 = (Inm2bJ * (Math.cos((n - 2) * t1) - Math.cosh(n * t1)) * n2 / (2 - 2 * n + n2)
                    + (-4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                    + 2 * n2 * (2 + n2 - 2 * n) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                    + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (Math.cos((n - 2) * t1) + Math.cosh(n * t1)) + 2 * Inm1bJ * (bJ * Math.cos((n - 2) * t1) + (n - 1) * Math.cosh(n * t1)))
            ) / (bJ * bJ * (4 + n4) * I0bJ)
            ) * bmu * bmu / 4 / n2;

            dAxc1 = 0.25 * bmu * bmu * (
                    -(n - 2) * Inm2bJ * Math.sin((n - 2) * t1) / (2 - 2 * n + n2)
                            +
                            (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                                    - 2 * n * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (bJ - 1 - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)

                            ) / (bJ * bJ * n * (4 + n4) * I0bJ)
            );

            Axs1 = ((n2 * Inm2bJ * Math.sin((n - 2) * t1)) / (2 - 2 * n + n2)
                    + (bJ * (bJ * n2 * (2 + 2 * n + n2) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) - (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                    + 2 * n2 * (2 - 2 * n + n2) * I0bJ * Math.sin((n + 2) * t1) * ((-bJ + bJ * bJ - n * bJ) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n2 + 2 * n) * InbJ)
            ) / (bJ * bJ * (4 + n4) * I0bJ)
            ) * bmu * bmu / 4 / n2;

            dAxs1 = 0.25 * bmu * bmu * (((n - 2) * Inm2bJ * Math.cos((n - 2) * t1)) / (2 - 2 * n + n2)
                    + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ)
                    - 4 * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
            )
                    + 2 * n * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
            ) / (bJ * bJ * n * (4 + n4) * I0bJ)
            );


            d2Axc0 = (-(1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.cos((n + 1) * t1)
                    + (-1 + n) * Inm1bJ * ((1 + n) * Math.cos((n - 1) * t1) - 2 * n3 * Math.cos((n - 1) * t1))
                    + InbJ * (-2 * (1 - n2 + 2 * n4) * Math.cos(t1) * Math.cos(n * t1) + 4 * n3 * Math.sin(t1) * Math.sin(n * t1))
            ) * bmu / (1 + 4 * n4);

            d2Axs0 = (-2 * bJ * bmu * Inm1bJ * (2 * n3 * Math.cos(n * t1) * Math.sin(t1) + (1 - n2 + 2 * n4) * Math.cos(t1) * Math.sin(n * t1))
                    + 2 * bmu * InbJ * (n * (1 - n2 - 2 * bJ * n2 + 2 * n3 + 2 * n4) * Math.cos(n * t1) * Math.sin(t1) + (n * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 2 * n4)) * Math.cos(t1) * Math.sin(n * t1))
            ) / (bJ + 4 * bJ * n4);

            d3Axc0 = (2 * InbJ * ((1 - n2 + 4 * n4) * Math.cos(n * t1) * Math.sin(t1) + n * (1 + n2 + 2 * n4) * Math.cos(t1) * Math.sin(n * t1))
                    + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.sin((n + 1) * t1)
                    - (-1 + n) * (-1 + n) * Inm1bJ * ((1 + n) * Math.sin((n - 1) * t1) + 2 * n3 * Math.sin((1 - n) * t1))
            ) * bmu / (1 + 4 * n4);

            d3Axs0 = -(bJ * Inm1bJ * (n * (1 + n2 + 2 * n4) * Math.cos(t1) * Math.cos(n * t1) + (-1 + n2 - 4 * n4) * Math.sin(t1) * Math.sin(n * t1))
                    + InbJ * (-n * Math.cos(t1) * Math.cos(n * t1) * ((1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) - bJ * (1 + n2 + 2 * n4)) + Math.sin(t1) * Math.sin(n * t1) * (n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 4 * n4)))
            ) * 2 * bmu / (bJ + 4 * bJ * n4);


        }


//        System.out.println("~~~~~~~~~~~~~~~Debug only ~~~~~~~~~~~~~~~~~~~~~");
//        System.out.println("nMax= " + nMax + ";");
//        System.out.println("bJ= " + bJ + ";");
//        System.out.println("ei={" + ei.getX(0) + "," + ei.getX(1) + "};");
//        System.out.println("ej={" + ej.getX(0) + "," + ej.getX(1) + "};");
        System.out.println("zero_order term= " + test0);

//        System.out.println("Axc0= " + Axc0);
//        System.out.println("Axs0= " + Axs0);
//        System.out.println("dAxc0= " + dAxc0);
//        System.out.println("dAxs0= " + dAxs0);
//        System.out.println("Axc1= " + Axc1);
//        System.out.println("dAxc1= " + dAxc1);
//        System.out.println("Axs1= " + Axs1);
//        System.out.println("dAxs1= " + dAxs1);
//        System.out.println("d2Axc0= " + d2Axc0);
//        System.out.println("d2Axs0= " + d2Axs0);
//        System.out.println("d3Axc0= " + d3Axc0);
        System.out.println("d3Axs0= " + d3Axs0);

        System.out.println("test1= " + test1);
        System.out.println("test2= " + test2);
        System.out.println("test3= " + test3);
//        System.out.println("test3_1= " + test3_1);
//        System.out.println("test3_2= " + test3_2);
        System.out.println("test4= " + test4);
        System.exit(2);


        Ayc0 = bmu * (I0bJ + I1bJ) * Math.sin(t1);
        dAyc0 = bmu * (I0bJ + I1bJ) * Math.cos(t1);
        Ays0 = 0;
        dAys0 = 0;
        Ayc1 = 0.25 * bmu * bmu * Math.sin(t1) * Math.sin(t1) * (I0bJ + 2 * I1bJ + I2bJ);
        dAyc1 = 0.25 * bmu * bmu * Math.sin(2 * t1) * (I0bJ + 2 * I1bJ + I2bJ);
        Ays1 = 0;
        dAys1 = 0;

        //y direction
        if (nMax > 0) {
            int n = nMax;
            int n2 = n * n;
            int n3 = n2 * n;
            int n4 = n2 * n2;
            double InbJ = besselI(n, bJ);
            double Inm1bJ = besselI(n - 1, bJ);
            double Inm2bJ = besselI(n - 2, bJ);
            double Inp1bJ = besselI(n + 1, bJ);
            Ayc0 = (2 * bmu * Math.cos(n * t1) * Math.sin(t1) * ((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * n2 + 2 * bJ * n2 - 2 * n3) * InbJ)
                    - 2 * bmu * n * Math.cos(t1) * Math.sin(n * t1) * (2 * bJ * Inm1bJ + (1 + 2 * bJ + 2 * (-1 + n) * n) * InbJ)
            ) / (bJ + 4 * bJ * n4);

            dAyc0 = 2 * bmu * ((bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
                    + n * Math.sin(t1) * Math.sin(n * t1) * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (n - bJ) * InbJ))
            ) / (bJ + 4 * bJ * n4);

            Ays0 = (bJ * (1 - 2 * n + 2 * n2) * Inp1bJ * (-Math.cos((n + 1) * t1) + Math.cosh(n * t1))
                    + bJ * Inm1bJ * ((1 + 2 * n + 2 * n2) * Math.cos((n - 1) * t1) + (-1 + 2 * n - 2 * n2) * Math.cosh(n * t1))
                    + 2 * InbJ * (2 * bJ * n * Math.cos(t1) * Math.cos(n * t1) + n * (1 - 2 * n + 2 * n2) * Math.cosh(n * t1) + bJ * (1 + 2 * n2) * Math.sin(t1) * Math.sin(n * t1))
            ) * bmu / (bJ + 4 * bJ * n4);

            dAys0 = ((1 + n - 2 * n3) * Inm1bJ * Math.sin((n - 1) * t1)
                    + 2 * InbJ * (n * (-1 + 2 * n2) * Math.cos(n * t1) * Math.sin(t1) + Math.cos(t1) * Math.sin(n * t1))
                    + (1 - n + 2 * n3) * Inp1bJ * Math.sin((n + 1) * t1)
            ) * bmu / (1 + 4 * n4);


            Ayc1 = (n2 * Inm2bJ * (-Math.cos((n - 2) * t1) + Math.cosh(n * t1)) / (2 - 2 * n + n2)
                    - (4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                    + 2 * n2 * (2 - 2 * n + n2) * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                    + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (Math.cos((n - 2) * t1) + Math.cosh(n * t1)) + 2 * Inm1bJ * (bJ * Math.cos((n - 2) * t1) + (-1 + n) * Math.cosh(n * t1)))
            ) / (bJ * bJ * (4 + n4) * I0bJ)
            ) * bmu * bmu / 4 / n2;

            dAyc1 = 0.25 * bmu * bmu * (((-2 + n) * Inm2bJ * Math.sin((n - 2) * t1)) / (2 - 2 * n + n2)
                    + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                    + 2 * n * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
            ) / (bJ * bJ * n * (4 + n4) * I0bJ)
            );

            Ays1 = (-bJ * n2 * (2 + n * (2 + n)) * I0bJ * Math.sin((n - 2) * t1) * ((-1 + bJ + n) * Inm1bJ + bJ * InbJ)
                    - 2 * bJ * (4 + n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                    - n2 * (2 - 2 * n + n2) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n * (1 + n)) * InbJ)
            ) * bmu * bmu / (2 * bJ * bJ * n2 * (4 + n4) * I0bJ);

            dAys1 = 0.25 * bmu * bmu * (-(-2 + n) * Inm2bJ * Math.cos((n - 2) * t1) / (2 - 2 * n + n2)
                    + (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ) - 4 * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                    - 2 * n * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
            ) / (bJ * bJ * n * (4 + n4) * I0bJ)
            );


        }


//        System.out.println("Ayc0= " + Ayc0);
//        System.out.println("dAyc0= " + dAyc0);
//        System.out.println("Ays0= " + Ays0);
//        System.out.println("dAys0= " + dAys0);
//        System.out.println("Ayc1= " + Ayc1);
//        System.out.println("dAyc1= " + dAyc1);
//        System.out.println("Ays1= " + Ays1);
//        System.out.println("dAys1= " + dAys1);


//        System.out.println("zero_order term= " + test0);
//        System.out.println("test1= " + test1);
//        System.out.println("test2= " + test2);
//        System.out.println("test3= " + test3);
//        System.out.println("test3_1= " + test3_1);
//        System.out.println("test3_2= " + test3_2);
//        System.out.println("test4= " + test4);
//        System.exit(2);

    }

    public void zeroSum() {
        Axc0 = 0;
    }

    public double getSum() {
        return Axc0;
    }
}
