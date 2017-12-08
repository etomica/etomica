package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

import static etomica.math.SpecialFunctions.besselI;

/**
 * Anx expressions for computing v_E and v_EE in the mapping.
 *
 * @author Weisong Lin
 */

public class PotentialCalculationHeisenberg implements PotentialCalculation {
    //public class PotentialCalculationHeisenberg {
    protected Vector ei, ej;
    protected double JEEMJEJE, UEE, VarJEMUE;
    protected final double mu, J, bt, bJ, bmu; //TODO should I add final here
    protected double[] Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    protected double[] d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
    protected double[] Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
    protected double[] d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;
    protected double psix1, psix2, psix11, psix12, psix22, psi1x1, psi1x2;
    protected double psiy1, psiy2, psiy11, psiy12, psiy22, psi1y1, psi1y2;
    protected int nMax = 0;

    public PotentialCalculationHeisenberg(Space space, double dipoleMagnitude, double interactionS, double beta, int nmax) {
        ei = space.makeVector();//TODO Do I have to do this again.
        ej = space.makeVector();
        J = interactionS;
        mu = dipoleMagnitude;
        bt = beta;
        bJ = bt * J;
        bmu = bt * mu;
        nMax = nmax;

        Axc0 = new double[nMax + 1];
        Axs0 = new double[nMax + 1];
        dAxc0 = new double[nMax + 1];
        dAxs0 = new double[nMax + 1];
        Axc1 = new double[nMax + 1];
        Axs1 = new double[nMax + 1];
        dAxc1 = new double[nMax + 1];
        dAxs1 = new double[nMax + 1];
        d2Axc0 = new double[nMax + 1];
        d2Axs0 = new double[nMax + 1];
        d3Axc0 = new double[nMax + 1];
        d3Axs0 = new double[nMax + 1];
        d2Axc1 = new double[nMax + 1];
        d2Axs1 = new double[nMax + 1];
        Ayc0 = new double[nMax + 1];
        Ays0 = new double[nMax + 1];
        dAyc0 = new double[nMax + 1];
        dAys0 = new double[nMax + 1];
        Ayc1 = new double[nMax + 1];
        Ays1 = new double[nMax + 1];
        dAyc1 = new double[nMax + 1];
        dAys1 = new double[nMax + 1];
        d2Ayc0 = new double[nMax + 1];
        d2Ays0 = new double[nMax + 1];
        d3Ayc0 = new double[nMax + 1];
        d3Ays0 = new double[nMax + 1];
        d2Ayc1 = new double[nMax + 1];
        d2Ays1 = new double[nMax + 1];

    }


    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }
        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        ei.E(atom1.getOrientation().getDirection());
        ej.E(atom2.getOrientation().getDirection());
        double t1 = Math.acos(ei.getX(0));
        double t2 = Math.acos(ej.getX(0));

        double I0bJ = besselI(0, bJ);
        double I1bJ = besselI(1, bJ);
        double I2bJ = besselI(2, bJ);

//        double test0 = 0;
//        double test1 = 0;
//        double test2 = 0;
//        double test3 = 0;
//        double test3_1 = 0;
//        double test3_2 = 0;
//        double test4 = 0;


//        System.out.println("~~~~~~~~~~~~~~~Debug only ~~~~~~~~~~~~~~~~~~~~~");
//        System.out.println("T=" + (1 / bt) + ";");
//        System.out.println("bt=" + bt + ";");
//        System.out.println("J=" + J + ";");
//        System.out.println("mu=" + mu + ";");
//        System.out.println("n= " + n + ";");
//        System.out.println("bJ= " + bJ + ";");
//        System.out.println("ei={" + ei.getX(0) + "," + ei.getX(1) + "};");
//        System.out.println("ej={" + ej.getX(0) + "," + ej.getX(1) + "};");
//        System.out.println("theta1 = ArcCos[ei[[1]]];");
//        System.out.println("theta2 = ArcCos[ej[[1]]];");
//        System.out.println("bt = bJ/J;");
//        System.out.println("bmu = bt*mu;");
        for (int n = 0; n <= nMax; n++) {
            if (n == 0) {
                Axc0[n] = bmu * (I0bJ + I1bJ) * (-1 + Math.cos(t1));
                dAxc0[n] = -bmu * (I0bJ + I1bJ) * Math.sin(t1);
                Axs0[n] = 0;
                dAxs0[n] = 0;
                Axc1[n] = -0.25 * bmu * bmu * Math.sin(t1) * Math.sin(t1) * (I0bJ + 2 * I1bJ + I2bJ);
                dAxc1[n] = -0.25 * bmu * bmu * (I0bJ + 2 * I1bJ + I2bJ) * Math.sin(2 * t1);
                Axs1[n] = 0;
                dAxs1[n] = 0;

                d2Axc0[n] = -bmu * Math.cos(t1) * (I0bJ + I1bJ);
                d2Axs0[n] = 0;
                d3Axc0[n] = bmu * Math.sin(t1) * (I0bJ + I1bJ);
                d3Axs0[n] = 0;
                d2Axc1[n] = -0.5 * bmu * bmu * Math.cos(2 * t1) * (I0bJ + 2 * I1bJ + I2bJ);
                d2Axs1[n] = 0;
            }

//        System.out.println("Axc00- " +"("+ Axc0 + ")");
//        System.out.println("dAxc00- "+"(" + dAxc0+ ")");
//        System.out.println("Axc10- "+"(" + Axc1+ ")");
//        System.out.println("dAxc10- "+"(" + dAxc1+ ")");
//        System.out.println("d2Axc00- "+"(" + d2Axc0+ ")");
//        System.out.println("d3Axc00- "+"(" + d3Axc0+ ")");
//        System.out.println("d2Axc10- "+"(" + d2Axc1+ ")");


            //x direction
            if (n > 0) {
                int n2 = n * n;
                int n3 = n2 * n;
                int n4 = n2 * n2;
                double InbJ = besselI(n, bJ);
                double Inm1bJ = besselI(n - 1, bJ);
                double Inm2bJ = besselI(n - 2, bJ);
                double Inp1bJ = besselI(n + 1, bJ);
                Axc0[n] = 2 * bmu * (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
                        + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * Math.sin(t1) * Math.sin(n * t1))
                        / (bJ + 4 * bJ * n4);

                dAxc0[n] = (2 * InbJ * (-Math.cos(n * t1) * Math.sin(t1) + (n - 2 * n3) * Math.cos(t1) * Math.sin(n * t1))
                        + Inm1bJ * (1 + n - 2 * n3) * Math.sin((n - 1) * t1)
                        + (-1 + n - 2 * n3) * Inp1bJ * Math.sin((n + 1) * t1)
                ) * bmu / (1 + 4 * n4);

                Axs0[n] = ((1 + 2 * n + 2 * n2) * Inm1bJ * Math.sin((n - 1) * t1)
                        + 2 * InbJ * (-2 * n * Math.cos(n * t1) * Math.sin(t1) + (1 + 2 * n2) * Math.cos(t1) * Math.sin(n * t1))
                        + (1 - 2 * n + 2 * n2) * Inp1bJ * Math.sin((n + 1) * t1)
                ) * bmu / (1 + 4 * n4);

                dAxs0[n] = (-2 * bmu * n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * Math.cos(t1) * Math.cos(n * t1)
                        - 2 * bmu * (bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * Math.sin(t1) * Math.sin(n * t1)
                ) / (bJ + 4 * bJ * n4);

                Axc1[n] = (Inm2bJ * (Math.cos((n - 2) * t1) - Math.cosh(n * t1)) * n2 / (2 - 2 * n + n2)
                        + (-4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 + n2 - 2 * n) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (Math.cos((n - 2) * t1) + Math.cosh(n * t1)) + 2 * Inm1bJ * (bJ * Math.cos((n - 2) * t1) + (n - 1) * Math.cosh(n * t1)))
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu * bmu / 4 / n2;

                dAxc1[n] = 0.25 * bmu * bmu * (
                        -(n - 2) * Inm2bJ * Math.sin((n - 2) * t1) / (2 - 2 * n + n2)
                                +
                                (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                                        - 2 * n * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (bJ - 1 - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)

                                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                Axs1[n] = ((n2 * Inm2bJ * Math.sin((n - 2) * t1)) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n2 * (2 + 2 * n + n2) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) - (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n2 * (2 - 2 * n + n2) * I0bJ * Math.sin((n + 2) * t1) * ((-bJ + bJ * bJ - n * bJ) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n2 + 2 * n) * InbJ)
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu * bmu / 4 / n2;

                dAxs1[n] = 0.25 * bmu * bmu * (((n - 2) * Inm2bJ * Math.cos((n - 2) * t1)) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ)
                        - 4 * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                d2Axc0[n] = (-(1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.cos((n + 1) * t1)
                        + (-1 + n) * Inm1bJ * ((1 + n) * Math.cos((n - 1) * t1) - 2 * n3 * Math.cos((n - 1) * t1))
                        + InbJ * (-2 * (1 - n2 + 2 * n4) * Math.cos(t1) * Math.cos(n * t1) + 4 * n3 * Math.sin(t1) * Math.sin(n * t1))
                ) * bmu / (1 + 4 * n4);

                d2Axs0[n] = (-2 * bJ * bmu * Inm1bJ * (2 * n3 * Math.cos(n * t1) * Math.sin(t1) + (1 - n2 + 2 * n4) * Math.cos(t1) * Math.sin(n * t1))
                        + 2 * bmu * InbJ * (n * (1 - n2 - 2 * bJ * n2 + 2 * n3 + 2 * n4) * Math.cos(n * t1) * Math.sin(t1) + (n * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 2 * n4)) * Math.cos(t1) * Math.sin(n * t1))
                ) / (bJ + 4 * bJ * n4);

                d3Axc0[n] = (2 * InbJ * ((1 - n2 + 4 * n4) * Math.cos(n * t1) * Math.sin(t1) + n * (1 + n2 + 2 * n4) * Math.cos(t1) * Math.sin(n * t1))
                        + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.sin((n + 1) * t1)
                        - (-1 + n) * (-1 + n) * Inm1bJ * ((1 + n) * Math.sin((n - 1) * t1) + 2 * n3 * Math.sin((1 - n) * t1))
                ) * bmu / (1 + 4 * n4);

                d3Axs0[n] = -(bJ * Inm1bJ * (n * (1 + n2 + 2 * n4) * Math.cos(t1) * Math.cos(n * t1) + (-1 + n2 - 4 * n4) * Math.sin(t1) * Math.sin(n * t1))
                        + InbJ * (-n * Math.cos(t1) * Math.cos(n * t1) * ((1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) - bJ * (1 + n2 + 2 * n4)) + Math.sin(t1) * Math.sin(n * t1) * (n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 4 * n4)))
                ) * 2 * bmu / (bJ + 4 * bJ * n4);

                d2Axc1[n] = 0.25 * bmu * bmu * (-(n - 2) * (n - 2) * Inm2bJ * Math.cos((n - 2) * t1) / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                d2Axs1[n] = 0.25 * bmu * bmu * (-((n - 2) * (n - 2) * Inm2bJ * Math.sin((n - 2) * t1)) / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (n + 2) * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


            }
        }


//        System.out.println("~~~~~~~~~~~~~~~Debug only ~~~~~~~~~~~~~~~~~~~~~");
//        System.out.println("nMax= " + nMax + ";");
//        System.out.println("bJ= " + bJ + ";");
//        System.out.println("ei={" + ei.getX(0) + "," + ei.getX(1) + "};");
//        System.out.println("ej={" + ej.getX(0) + "," + ej.getX(1) + "};");
//        System.out.println("zero_order term= " + test0);

//        System.out.println("Axc0- " + "(" + Axc0 + ")");
//        System.out.println("Axs0- " + "(" + Axs0 + ")");
//        System.out.println("dAxc0- " + "(" + dAxc0 + ")");
//        System.out.println("dAxs0- " + "(" + dAxs0 + ")");
//        System.out.println("Axc1- " + "(" + Axc1 + ")");
//        System.out.println("dAxc1- " + "(" + dAxc1 + ")");
//        System.out.println("Axs1- " + "(" + Axs1 + ")");
//        System.out.println("dAxs1- " + "(" + dAxs1 + ")");
//        System.out.println("d2Axc0- " + "(" + d2Axc0 + ")");
//        System.out.println("d2Axs0- " + "(" + d2Axs0 + ")");
//        System.out.println("d3Axc0- " + "(" + d3Axc0 + ")");
//        System.out.println("d3Axs0- " + "(" + d3Axs0 + ")");
//        System.out.println("d2Axc1- " + "(" + d2Axc1 + ")");
//        System.out.println("d2Axs1- " + "(" + d2Axs1 + ")");

//        System.out.println("test1= " + test1);
//        System.out.println("test2= " + test2);
//        System.out.println("test3= " + test3);
//        System.out.println("test3_1= " + test3_1);
//        System.out.println("test3_2= " + test3_2);
//        System.out.println("test4= " + test4);
//        System.exit(2);

        for (int n = 0; n <= nMax; n++) {
            if (n == 0) {
                Ayc0[n] = bmu * (I0bJ + I1bJ) * Math.sin(t1);
                dAyc0[n] = bmu * (I0bJ + I1bJ) * Math.cos(t1);
                Ays0[n] = 0;
                dAys0[n] = 0;
                Ayc1[n] = 0.25 * bmu * bmu * Math.sin(t1) * Math.sin(t1) * (I0bJ + 2 * I1bJ + I2bJ);
                dAyc1[n] = 0.25 * bmu * bmu * Math.sin(2 * t1) * (I0bJ + 2 * I1bJ + I2bJ);
                Ays1[n] = 0;
                dAys1[n] = 0;

                d2Ayc0[n] = -bmu * Math.sin(t1) * (I0bJ + I1bJ);
                d2Ays0[n] = 0;
                d3Ayc0[n] = -bmu * Math.cos(t1) * (I0bJ + I1bJ);
                d3Ays0[n] = 0;
                d2Ayc1[n] = 0.5 * bmu * bmu * (I0bJ + 2 * I1bJ + I2bJ) * Math.cos(2 * t1);
                d2Ays1[n] = 0;
            }
//        System.out.println("Ayc00- " + "(" + Ayc0 + ")");
//        System.out.println("dAyc00- " + "(" + dAyc0 + ")");
//        System.out.println("Ayc10- " + "(" + Ayc1 + ")");
//        System.out.println("dAyc10- " + "(" + dAyc1 + ")");
//        System.out.println("d2Ayc00- " + "(" + d2Ayc0 + ")");
//        System.out.println("d3Ayc00- " + "(" + d3Ayc0 + ")");
//        System.out.println("d2Ayc10- " + "(" + d2Ayc1 + ")");

            //y direction
            if (n > 0) {
                int n2 = n * n;
                int n3 = n2 * n;
                int n4 = n2 * n2;
                double InbJ = besselI(n, bJ);
                double Inm1bJ = besselI(n - 1, bJ);
                double Inm2bJ = besselI(n - 2, bJ);
                double Inp1bJ = besselI(n + 1, bJ);
                Ayc0[n] = (2 * bmu * Math.cos(n * t1) * Math.sin(t1) * ((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * n2 + 2 * bJ * n2 - 2 * n3) * InbJ)
                        - 2 * bmu * n * Math.cos(t1) * Math.sin(n * t1) * (2 * bJ * Inm1bJ + (1 + 2 * bJ + 2 * (-1 + n) * n) * InbJ)
                ) / (bJ + 4 * bJ * n4);

                dAyc0[n] = 2 * bmu * ((bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
                        + n * Math.sin(t1) * Math.sin(n * t1) * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (n - bJ) * InbJ))
                ) / (bJ + 4 * bJ * n4);

                Ays0[n] = (bJ * (1 - 2 * n + 2 * n2) * Inp1bJ * (-Math.cos((n + 1) * t1) + Math.cosh(n * t1))
                        + bJ * Inm1bJ * ((1 + 2 * n + 2 * n2) * Math.cos((n - 1) * t1) + (-1 + 2 * n - 2 * n2) * Math.cosh(n * t1))
                        + 2 * InbJ * (2 * bJ * n * Math.cos(t1) * Math.cos(n * t1) + n * (1 - 2 * n + 2 * n2) * Math.cosh(n * t1) + bJ * (1 + 2 * n2) * Math.sin(t1) * Math.sin(n * t1))
                ) * bmu / (bJ + 4 * bJ * n4);

                dAys0[n] = ((1 + n - 2 * n3) * Inm1bJ * Math.sin((n - 1) * t1)
                        + 2 * InbJ * (n * (-1 + 2 * n2) * Math.cos(n * t1) * Math.sin(t1) + Math.cos(t1) * Math.sin(n * t1))
                        + (1 - n + 2 * n3) * Inp1bJ * Math.sin((n + 1) * t1)
                ) * bmu / (1 + 4 * n4);


                Ayc1[n] = (n2 * Inm2bJ * (-Math.cos((n - 2) * t1) + Math.cosh(n * t1)) / (2 - 2 * n + n2)
                        - (4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 - 2 * n + n2) * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (Math.cos((n - 2) * t1) + Math.cosh(n * t1)) + 2 * Inm1bJ * (bJ * Math.cos((n - 2) * t1) + (-1 + n) * Math.cosh(n * t1)))
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu * bmu / 4 / n2;

                dAyc1[n] = 0.25 * bmu * bmu * (((-2 + n) * Inm2bJ * Math.sin((n - 2) * t1)) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                Ays1[n] = (-bJ * n2 * (2 + n * (2 + n)) * I0bJ * Math.sin((n - 2) * t1) * ((-1 + bJ + n) * Inm1bJ + bJ * InbJ)
                        - 2 * bJ * (4 + n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        - n2 * (2 - 2 * n + n2) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n * (1 + n)) * InbJ)
                ) * bmu * bmu / (2 * bJ * bJ * n2 * (4 + n4) * I0bJ);

                dAys1[n] = 0.25 * bmu * bmu * (-(-2 + n) * Inm2bJ * Math.cos((n - 2) * t1) / (2 - 2 * n + n2)
                        + (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ) - 4 * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                d2Ayc0[n] = -2 * bmu * (bJ * Inm1bJ * ((1 - n2 + 2 * n4) * Math.cos(n * t1) * Math.sin(t1) + 2 * n3 * Math.cos(t1) * Math.sin(n * t1))
                        + InbJ * ((bJ - n - bJ * bJ * n2 + n3 - 2 * n4 + 2 * bJ * n4 - 2 * n2 * n3) * Math.cos(n * t1) * Math.sin(t1) + n * (-1 + n2 + 2 * bJ * n2 - 2 * n3 - 2 * n4) * Math.cos(t1) * Math.sin(n * t1))
                ) / (bJ + 4 * bJ * n4);

                d2Ays0[n] = -bmu * ((n - 1) * (n - 1) * (1 + 2 * n + 2 * n2) * Inm1bJ * Math.cos((n - 1) * t1)
                        - (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.cos((n + 1) * t1)
                        + 2 * InbJ * (-2 * n3 * Math.cos(t1) * Math.cos(n * t1) + (1 - n2 + 2 * n4) * Math.sin(t1) * Math.sin(n * t1))
                ) / (1 + 4 * n4);

                d3Ayc0[n] = 2 * bmu * (bJ * Inm1bJ * (-(1 - n2 + 4 * n4) * Math.cos(t1) * Math.cos(n * t1) + n * (1 + n2 + 2 * n4) * Math.sin(t1) * Math.sin(n * t1))
                        + InbJ * (-(-n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (1 - n2 + 4 * n4)) * Math.cos(t1) * Math.cos(n * t1) + n * (-(1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (1 + n2 + 2 * n4)) * Math.sin(t1) * Math.sin(n * t1))
                ) / (bJ + 4 * bJ * n4);


                d3Ays0[n] = -bmu * (-(-1 + n) * (-1 + n) * (-1 + n) * (1 + 2 * n + 2 * n2) * Inm1bJ * Math.sin((n - 1) * t1)
                        + 2 * InbJ * (n * (1 + n2 + 2 * n4) * Math.cos(n * t1) * Math.sin(t1) + (1 - n2 + 4 * n4) * Math.cos(t1) * Math.sin(n * t1))
                        + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.sin((n + 1) * t1)
                ) / (1 + 4 * n4);

                d2Ayc1[n] = 0.25 * bmu * bmu * ((-2 + n) * (-2 + n) * Inm2bJ * Math.cos((n - 2) * t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * (-2 + n) * (-2 + n) * n * (2 + 2 * n + n2) * I0bJ * Math.cos((n - 2) * t1) * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                d2Ays1[n] = 0.25 * bmu * bmu * ((n - 2) * (n - 2) * Inm2bJ * Math.sin((n - 2) * t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * Math.sin((n - 2) * t1) * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * Math.sin(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * Math.sin((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );
            }
        }

        System.out.println("~~~~~~~~~~~~~~~Debug only ~~~~~~~~~~~~~~~~~~~~~");
        System.out.println("nMax= " + nMax + ";");
        System.out.println("bJ= " + bJ + ";");
        System.out.println("ei={" + ei.getX(0) + "," + ei.getX(1) + "};");
        System.out.println("ej={" + ej.getX(0) + "," + ej.getX(1) + "};");
//        System.exit(2);

//        System.out.println("Ayc0- " + "(" + Ayc0 + ")");
//        System.out.println("Ays0- " + "(" + Ays0 + ")");
//        System.out.println("dAyc0- " + "(" + dAyc0 + ")");
//        System.out.println("dAys0- " + "(" + dAys0 + ")");
//        System.out.println("Ayc1- " + "(" + Ayc1 + ")");
//        System.out.println("dAyc1- " + "(" + dAyc1 + ")");
//        System.out.println("Ays1- " + "(" + Ays1 + ")");
//        System.out.println("dAys1- " + "(" + dAys1 + ")");
//        System.out.println("d2Ayc0- " + "(" + d2Ayc0 + ")");
//        System.out.println("d2Ays0- " + "(" + d2Ays0 + ")");
//        System.out.println("d3Ayc0- " + "(" + d3Ayc0 + ")");
//        System.out.println("d3Ays0- " + "(" + d3Ays0 + ")");
//        System.out.println("d2Ayc1- " + "(" + d2Ayc1 + ")");
//        System.out.println("d2Ays1- " + "(" + d2Ays1 + ")");


//        System.out.println("test1= " + test1);
//        System.out.println("test2= " + test2);
//        System.out.println("test3= " + test3);
//        System.out.println("test3_1= " + test3_1);
//        System.out.println("test3_2= " + test3_2);
//        System.out.println("test4= " + test4);
//        System.exit(2);

        double p0 = Math.exp(bJ * Math.cos(t1 - t2));
        double pM1 = 1 / p0;
        double pM2 = pM1 / p0;
        double lnp1 = -bJ * Math.sin(t1 - t2);
        double lnp11 = -bJ * Math.cos(t1 - t2);
        double lnp2 = -lnp1;
        double px1 = bmu * p0 * (Math.cos(t1) + Math.cos(t2));
        double py1 = bmu * p0 * (Math.sin(t1) + Math.sin(t2));
        for (int n = 0; n <= nMax; n++) {
//             psix1, psix2,psix11,psix12,psix22,psi1x1,psi1x2;
//             psiy1, psiy2,psiy11,psiy12,psiy22,psi1y1,psi1y2;
            psix1 += dAxs0[n] * Math.sin(n * t2) + dAxc0[n] * Math.cos(n * t2);
            psiy1 += dAys0[n] * Math.sin(n * t2) + dAyc0[n] * Math.cos(n * t2);

            psix2 += n * Axs0[n] * Math.cos(n * t2) - n * Axc0[n] * Math.sin(n * t2);
            psiy2 += n * Ays0[n] * Math.cos(n * t2) - n * Ayc0[n] * Math.sin(n * t2);

            psi1x1 += dAxs1[n] * Math.sin(n * t2) + dAxc1[n] * Math.cos(n * t2);
            psi1y1 += dAys1[n] * Math.sin(n * t2) + dAyc1[n] * Math.cos(n * t2);


            psi1x2 += n * Axs1[n] * Math.cos(n * t2) - n * Axc1[n] * Math.cos(n * t2);
            psi1y2 += n * Ays1[n] * Math.cos(n * t2) - n * Ayc1[n] * Math.cos(n * t2);

            psix11 += d2Axs0[n] * Math.sin(n * t2) + d2Axc0[n] * Math.cos(n * t2);
            psiy11 += d2Ays0[n] * Math.sin(n * t2) + d2Ayc0[n] * Math.cos(n * t2);

            psix12 += n * dAxs0[n] * Math.cos(n * t2) - n * dAxc0[n] * Math.sin(n * t2);
            psiy12 += n * dAys0[n] * Math.cos(n * t2) - n * dAyc0[n] * Math.sin(n * t2);

            psix22 += -n * n * Axs0[n] * Math.sin(n * t2) - n * n * Axc0[n] * Math.cos(n * t2);
            psiy22 += -n * n * Ays0[n] * Math.sin(n * t2) - n * n * Ayc0[n] * Math.cos(n * t2);
        }


        System.out.println("psix1- " + "(" + psix1 + ")");
        System.out.println("psiy1- " + "(" + psiy1 + ")");
        System.out.println("psix2- " + "(" + psix2 + ")");
        System.out.println("psiy2- " + "(" + psix2 + ")");
        System.out.println("psi1x1- " + "(" + psi1x1 + ")");
        System.out.println("psi1y1- " + "(" + psi1y1 + ")");
        System.out.println("psi1x2- " + "(" + psi1x2 + ")");
        System.out.println("psi1y2- " + "(" + psi1y2 + ")");
        System.out.println("psix11- " + "(" + psix11 + ")");
        System.out.println("psiy11- " + "(" + psiy11 + ")");
        System.out.println("psix22- " + "(" + psix22 + ")");
        System.out.println("psiy22- " + "(" + psiy22 + ")");

        System.exit(2);


        double vEx1 = psix1 / p0;
        double vEx2 = psix2 / p0;
        double vEy1 = psiy1 / p0;
        double vEy2 = psiy2 / p0;

        double vEEx1 = psi1x1 / p0 - px1 * vEx1 / p0;
        double vEEx2 = psi1x2 / p0 - px1 * vEx2 / p0;

        double vEEy1 = psi1y1 / p0 - py1 * vEy1 / p0;
        double vEEy2 = psi1y2 / p0 - py1 * vEy2 / p0;

        //divergence of vEE
        //x
        JEEMJEJE += bmu * bmu * (1 + I1bJ / I0bJ) - (vEEx1 * lnp1 + vEEx2 * lnp2) + bmu * pM1 * (psix1 * Math.sin(t1) + psix2 * Math.sin(t2));
        //y
        JEEMJEJE += bmu * bmu * (1 + I1bJ / I0bJ) - (vEEy1 * lnp1 + vEEy2 * lnp2) - bmu * pM1 * (psiy1 * Math.cos(t1) + psiy2 * Math.cos(t2));

        //vE dot Grad vE
        //x
        JEEMJEJE += bmu * pM1 * (psix1 * Math.sin(t1) + psix2 * Math.sin(t2))
                - pM2 * lnp1 * (psix1 * psix11 - psix1 * psix12 + psix2 * psix12 - psix2 * psix22)
                - pM2 * (lnp11 - lnp1 * lnp1) * (psix1 - psix2) * (psix1 - psix2);
        //y
        JEEMJEJE += -bmu * pM1 * (psix1 * Math.cos(t1) + psix2 * Math.cos(t2))
                - pM2 * lnp1 * (psiy1 * psiy11 - psiy1 * psiy12 + psiy2 * psiy12 - psiy2 * psiy22)
                - pM2 * (lnp11 - lnp1 * lnp1) * (psiy1 - psiy2) * (psiy1 - psiy2);


        //Var<JE-UE>
        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
        Vector[][] t = potentialSecondDerivative.gradientAndTorque(atoms);//TODO
        double f1 = t[1][0].getX(0);
        double f2 = t[1][1].getX(0);
//        System.out.println("f1 = " + f1);

        double JEMUEx = bJ * Math.sin(t1 - t2) * pM1 * (psix1 - psix2) - (vEx1 * f1 + vEx2 * f2);
        double JEMUEy = bJ * Math.sin(t1 - t2) * pM1 * (psiy1 - psiy2) - (vEx1 * f1 + vEx2 * f2);

        VarJEMUE += JEMUEx * JEMUEx + JEMUEy * JEMUEy;
        //TODO how to return variance for both direction maybe I should return both x and y value?


        //PhiSum
//        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
        Tensor[] phi = potentialSecondDerivative.secondDerivative(atoms);
        double p12 = phi[0].component(0, 0);//TODO is it the right component? Do I also need to change the p2spin potential part!!
        double p11 = phi[1].component(0, 0);
        double p22 = phi[2].component(0, 0);

        double vDotGradvx1 = pM2 * (-lnp1 * psix1 * (psix1 - psix2) + psix1 * psix11 + psix2 * psix12);
        double vDotGradvy1 = pM2 * (-lnp1 * psiy1 * (psiy1 - psiy2) + psiy1 * psiy11 + psiy2 * psiy12);

        double vDotGradvx2 = pM2 * (-lnp1 * psix2 * (psix1 - psix2) + psix1 * psix12 + psix2 * psix22);
        double vDotGradvy2 = pM2 * (-lnp1 * psiy2 * (psiy1 - psiy2) + psiy1 * psiy12 + psiy2 * psiy22);

        double fxE1 = -bmu * Math.sin(t1);
        double fxE2 = -bmu * Math.sin(t2);
        double fyE1 = bmu * Math.cos(t1);
        double fyE2 = bmu * Math.cos(t2);

        UEE += vEEx1 * f1 + vEEx2 * f2 + vDotGradvx1 * f1 + vDotGradvx2 * f2
                + vEx1 * (p11 * vEx1 + p12 * vEx2) + vEx2 * (p12 * vEx1 + p22 * vEx2)
                + 2 * vEx1 * fxE1 + 2 * vEx2 * fxE2;

        UEE += vEEy1 * f1 + vEEy2 * f2 + vDotGradvy1 * f1 + vDotGradvy2 * f2
                + vEy1 * (p11 * vEy1 + p12 * vEy2) + vEy2 * (p12 * vEy1 + p22 * vEy2)
                + 2 * vEy1 * fyE1 + 2 * vEy2 * fyE2;

    }

    public void zeroSumJEEMJEJE() {
        JEEMJEJE = 0;
    }

    public double getSumJEEMJEJE() {
        return JEEMJEJE;
    }

    public void zeroSumUEE() {
        UEE = 0;
    }

    public double getSumUEE() {
        return UEE;
    }

    public void zeroSumVarJEMUE() {
        VarJEMUE = 0;
    }

    public double getSumVarJEMUE() {
        return VarJEMUE;
    }

}
