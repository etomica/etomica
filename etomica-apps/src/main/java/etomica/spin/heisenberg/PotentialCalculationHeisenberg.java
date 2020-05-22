package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialAtomicSecondDerivative;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
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
    protected double AEEJ0, JEMUExIdeal, JEMUEyIdeal, JEMUEIdealSquare, JEEMJEJE, UEE, JEMUExSquare, JEMUEySquare, JEMUEx, JEMUEy, dipolex, dipoley, JEEMJEJExtrying, UEEnow, JEMUE, dipoleconv;
    protected final double mu, J, bt, bJ, bmu; //TODO should I add final here
    protected double[] Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    protected double[] d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
    protected double[] Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
    protected double[] d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;
    protected double psix1, psix2, psix11, psix12, psix22, psi1x1, psi1x2, psi1x11, psi1x22, psix111, psix222, psix221, psix112;
    protected double psiy1, psiy2, psiy11, psiy12, psiy22, psi1y1, psi1y2, psi1y11, psi1y22, psiy111, psiy222, psiy221, psiy112;
    protected int nMax;
    protected int count = 1;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;

    public PotentialCalculationHeisenberg(Space space, double dipoleMagnitude, double interactionS, double beta, int nMax, AtomLeafAgentManager<MoleculeAgent> leafAgentManager) {
        ei = space.makeVector();//TODO Do I have to do this again.
        ej = space.makeVector();
        J = interactionS;
        mu = dipoleMagnitude;
        bt = beta;
        bJ = bt * J;
        bmu = bt * mu;
        this.nMax = nMax;
        this.leafAgentManager = leafAgentManager;

//        int nM = leafAgentManager.getBox().getLeafList().getAtomCount();
//        JEMUEx = new double[nM + 1];
//        JEMUEy = new double[nM + 1];

//        System.out.println(nM+1);
//        System.exit(2);

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
//        System.out.println(ei);
        double t1 = Math.atan2(ei.getX(1), ei.getX(0));
        double t2 = Math.atan2(ej.getX(1), ej.getX(0));


        count += 1;
//        System.out.println("(*" + count + "th term*)");
        boolean debug = (count % 150) == 0 && false;
        if (debug) {
            System.out.println("nMax= " + nMax + ";");
            System.out.println("mu= " + mu + ";");
            System.out.println("J= " + J + ";");
            System.out.println("bJ= " + bJ + ";");
            System.out.println("bt = bJ/J;");
            System.out.println("bmu = bt*mu;");
            System.out.println("t1 = " + t1 + ";");
            System.out.println("t2 = " + t2 + ";");
        }
        if (count > 2000 && debug) System.exit(2);


//        System.out.println("ei={" + ei.getX(0) + "," + ei.getX(1) + "};");
//        System.out.println("ej={" + ej.getX(0) + "," + ej.getX(1) + "};");
//        System.out.println("theta1 = " + t1 + ";");
//        System.out.println("theta2 = " + t2 + ";");
//        System.out.println("t1 = " + t1 + " c1=" + c1);
//System.out.println("Acos eix = " + Math.acos(ei.getX(0)));

        double cost1 = Math.cos(t1);
        double sint1 = Math.sin(t1);
        double sint1p2 = sint1 * sint1;
        double bJ2 = bJ * bJ;
        double bmu2 = bmu * bmu;
        double I0bJ = besselI(0, bJ);
        double I1bJ = besselI(1, bJ);
        double I2bJ = besselI(2, bJ);


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

                //passed the test
//                double Axc00 = bmu * (I0bJ + I1bJ) * (-1 + cost1);
//                System.out.println(Axc00 - Axc0[n]);
//                double dAxc00 = -(bmu * (I0bJ + I1bJ) * sint1);
//                System.out.println(dAxc00 - dAxc0[n]);
//                double Axc10 = -(bmu2 * (I0bJ + 2 * I1bJ + I2bJ) * sint1p2) / 4;
//                System.out.println(Axc10 - Axc1[n]);
//                double dAxc10 = -(bmu2 * (I0bJ + 2 * I1bJ + besselI(2, bJ)) * Math.sin(2 * t1)) / 4;
//                System.out.println(dAxc10 - dAxc1[n]);
//                double d2Axc00 = -(bmu * (I0bJ + I1bJ) * cost1);
//                System.out.println(d2Axc00 - d2Axc0[n]);
//                double d3Axc00 = bmu * (I0bJ + I1bJ) * sint1;
//                System.out.println(d3Axc00 - d3Axc0[n]);
//                double d2Axc10 = -(bmu2 * (I0bJ + 2 * I1bJ + I2bJ) * Math.cos(2 * t1)) / 2;
//                System.out.println(d2Axc10 - d2Axc1[n]);
//                System.exit(2);

            }


            //x direction
            if (n > 0) {
                int n2 = n * n;
                int n3 = n2 * n;
                int n4 = n2 * n2;
                double np1p3 = (n + 1) * (n + 1) * (n + 1);
                double np1p2 = (n + 1) * (n + 1);
                double nm2p2 = (n - 2) * (n - 2);
                double InbJ = besselI(n, bJ);
                double Inm1bJ = besselI(n - 1, bJ);//TODO
                double Inm2bJ = besselI(n - 2, bJ);
                double Inp1bJ = besselI(n + 1, bJ);
                double sinnt1 = Math.sin(n * t1);
                double cosnt1 = Math.cos(n * t2);
                double sinnm1t1 = Math.sin((n - 1) * t1);
                double sinnp1t1 = Math.sin((n + 1) * t1);
                double coshnt1 = Math.cosh(n * t1);

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
                ) * bmu * bmu / 4.0 / n2;

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
                ) * bmu * bmu / 4.0 / n2;

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
                ) * bmu / (1.0 + 4 * n4);

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


//                double Axc0test = (2 * bmu * (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * cost1 * cosnt1 + n * (2 * bJ * Inm1bJ + (1 + 2 * bJ + 2 * (-1 + n) * n) * InbJ) * sint1 * sinnt1)) / (bJ + 4 * bJ * n4);
//                System.out.println(Axc0test - Axc0[n]);
//
//
//                double dAxc0test = (bmu * (2 * InbJ * (-(cosnt1 * sint1) + n * (1 - 2 * n2) * cost1 * sinnt1) + Inm1bJ * ((1 + n) * sinnm1t1 + 2 * n3 * sinnm1t1 + (-1 + n - 2 * n3) * Inp1bJ * sinnp1t1))) / (1 + 4 * n4);
//                System.out.println(dAxc0test - dAxc0[n]);
//
//
//                double dAxc1test = (bmu2 * (-(((-2 + n) * besselI(-2 + n, bJ) * Math.sin((-2 + n) * t1)) / (2 + (-2 + n) * n)) + (bJ * (-(bJ * n * (-4 - 2 * n + n3) * I0bJ * (2 * Inm1bJ + InbJ) * Math.sin((-2 + n) * t1)) + 4 * (4 + n4) * (bJ * I1bJ * InbJ + I0bJ * (-(bJ * Inm1bJ) + n * InbJ)) * sinnt1) - 2 * n * (4 - 2 * n + n3) * I0bJ * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ2 - 2 * bJ * n + 2 * n * (1 + n)) * InbJ) * Math.sin((2 + n) * t1)) / (bJ2 * n * (4 + n4) * I0bJ))) / 4;
//                System.out.println(dAxc1test - dAxc1[n]);
//
//
//                double Axs1test = (bmu2 * ((n2 * besselI(-2 + n, bJ) * Math.sin((-2 + n) * t1)) / (2 + (-2 + n) * n) + (bJ * (bJ * n2 * (2 + n * (2 + n)) * I0bJ * (2 * Inm1bJ + InbJ) * Math.sin((-2 + n) * t1) - 4 * (4 + n4) * (bJ * I1bJ * InbJ + I0bJ * (-(bJ * Inm1bJ) + n * InbJ)) * sinnt1) + 2 * n2 * (2 + (-2 + n) * n) * I0bJ * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ2 - 2 * bJ * n + 2 * n * (1 + n)) * InbJ) * Math.sin((2 + n) * t1)) / (bJ2 * (4 + n4) * I0bJ))) / (4 * n2);
//                System.out.println(Axs1test - Axs1[n]);
//
//
//                double dAxs1test = (bmu2 * (((-2 + n) * besselI(-2 + n, bJ) * Math.cos((-2 + n) * t1)) / (2 + (-2 + n) * n) + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * (2 * Inm1bJ + InbJ) * Math.cos((-2 + n) * t1) - 4 * (4 + n4) * (bJ * I1bJ * InbJ + I0bJ * (-(bJ * Inm1bJ) + n * InbJ)) * cosnt1) + 2 * n * (4 - 2 * n + n3) * I0bJ * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ2 - 2 * bJ * n + 2 * n * (1 + n)) * InbJ) * Math.cos((2 + n) * t1)) / (bJ2 * n * (4 + n4) * I0bJ))) / 4;
//                System.out.println(dAxs1test - dAxs1[n]);
//
//
//                double d2Axc0test = (bmu * (-(np1p2 * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.cos((1 + n) * t1)) + (-1 + n) * Inm1bJ * ((1 + n) * Math.cos((-1 + n) * t1) - 2 * n3 * Math.cos(t1 - n * t1)) + InbJ * (-2 * (1 - n2 + 2 * n4) * cost1 * cosnt1 + 4 * n3 * sint1 * sinnt1))) / (1 + 4 * n4);
//                System.out.println(d2Axc0test - d2Axc0[n]);
//
//
//                double d2Axs0test = (-2 * bJ * bmu * Inm1bJ * (2 * n3 * cosnt1 * sint1 + (1 - n2 + 2 * n4) * cost1 * sinnt1) + 2 * bmu * InbJ * (n * (1 - (1 + 2 * bJ) * n2 + 2 * n3 + 2 * n4) * cosnt1 * sint1 + (n * np1p2 * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 2 * n4)) * cost1 * sinnt1)) / (bJ + 4 * bJ * n4);
//                System.out.println(d2Axs0test - d2Axs0[n]);
//
//
//                double d2Axc1test = (bmu2 * (-((nm2p2 * besselI(-2 + n, bJ) * Math.cos((-2 + n) * t1)) / (2 + (-2 + n) * n)) + (bJ * (-(bJ * nm2p2 * n * (2 + 2 * n + n2) * I0bJ * (2 * Inm1bJ + InbJ) * Math.cos((-2 + n) * t1)) + 4 * n * (4 + n4) * (bJ * I1bJ * InbJ + I0bJ * (-(bJ * Inm1bJ) + n * InbJ)) * cosnt1) - 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ2 - 2 * bJ * n + 2 * n * (1 + n)) * InbJ) * Math.cos((2 + n) * t1)) / (bJ2 * n * (4 + n4) * I0bJ))) / 4;
//                System.out.println(d2Axc1test - d2Axc1[n]);
//
//
//                double d2Axs1test = (bmu2 * (-((nm2p2 * besselI(-2 + n, bJ) * Math.sin((-2 + n) * t1)) / (2 + (-2 + n) * n)) + (bJ * (-(bJ * nm2p2 * n * (2 + 2 * n + n2) * I0bJ * (2 * Inm1bJ + InbJ) * Math.sin((-2 + n) * t1)) + 4 * n * (4 + n4) * (bJ * I1bJ * InbJ + I0bJ * (-(bJ * Inm1bJ) + n * InbJ)) * sinnt1) - 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ2 - 2 * bJ * n + 2 * n * (1 + n)) * InbJ) * Math.sin((2 + n) * t1)) / (bJ2 * n * (4 + n4) * I0bJ))) / 4;
//                System.out.println(d2Axs1test - d2Axs1[n]);
//
//                System.exit(2);
//
//


            }
        }


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

                dAyc0[n] = 2.0 * bmu * ((bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * Math.cos(t1) * Math.cos(n * t1)
                        + n * Math.sin(t1) * Math.sin(n * t1) * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (n - bJ) * InbJ))
                ) / (bJ + 4 * bJ * n4);

                Ays0[n] = (bJ * (1 - 2 * n + 2 * n2) * Inp1bJ * (-Math.cos((n + 1) * t1) + Math.cosh(n * t1))
                        + bJ * Inm1bJ * ((1 + 2 * n + 2 * n2) * Math.cos((n - 1) * t1) + (-1 + 2 * n - 2 * n2) * Math.cosh(n * t1))
                        + 2 * InbJ * (2 * bJ * n * Math.cos(t1) * Math.cos(n * t1) + n * (1 - 2 * n + 2 * n2) * Math.cosh(n * t1) + bJ * (1 + 2 * n2) * Math.sin(t1) * Math.sin(n * t1))
                ) * bmu / (bJ + 4 * bJ * n4);

                dAys0[n] = ((1 + n - 2 * n3) * Inm1bJ * Math.sin((n - 1) * t1)
                        + 2 * InbJ * (n * (-1 + 2 * n2) * Math.cos(n * t1) * Math.sin(t1) + Math.cos(t1) * Math.sin(n * t1))
                        + (1 - n + 2 * n3) * Inp1bJ * Math.sin((n + 1) * t1)
                ) * bmu / (1.0 + 4 * n4);


                Ayc1[n] = (n2 * Inm2bJ * (-Math.cos((n - 2) * t1) + Math.cosh(n * t1)) / (2 - 2 * n + n2)
                        - (4 * bJ * (4 + n4) * Math.cos(n * t1) * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 - 2 * n + n2) * I0bJ * Math.cos((n + 2) * t1) * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
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


                //     d2Ayc0[n] = -2 * bmu * (bJ * Inm1bJ * ((1 - n2 + 2 * n4) * Math.cos(n * t1) * Math.sin(t1) + 2 * n3 * Math.cos(t1) * Math.sin(n * t1))
                //      + InbJ * ((bJ - n - bJ * bJ * n2 + n3 - 2 * n4 + 2 * bJ * n4 - 2 * n2 * n3) * Math.cos(n * t1) * Math.sin(t1) + n * (-1 + n2 + 2 * bJ * n2 - 2 * n3 - 2 * n4) * Math.cos(t1) * Math.sin(n * t1))
                //
                //     ) / (bJ + 4 * bJ * n4);

                //    d2Ays0[n] = -bmu * ((n - 1) * (n - 1) * (1 + 2 * n + 2 * n2) * Inm1bJ * Math.cos((n - 1) * t1)
                //           - (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * Math.cos((n + 1) * t1)
                //           + 2 * InbJ * (-2 * n3 * Math.cos(t1) * Math.cos(n * t1) + (1 - n2 + 2 * n4) * Math.sin(t1) * Math.sin(n * t1))
                //   ) / (1 + 4 * n4);


                d2Ayc0[n] = 2 * bmu * (-(bJ * Inm1bJ + (bJ + n * (-1 + n - 2 * n3)) * InbJ) * Math.cos(n * t1) * Math.sin(t1) + n2 * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * Math.cos(n * t1) * Math.sin(t1) - n * (bJ * Inm1bJ + (bJ + n * (-1 + n - 2 * n3)) * InbJ) * Math.cos(t1) * Math.sin(n * t1) + n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * Math.cos(t1) * Math.sin(n * t1)) / (bJ + 4 * bJ * n4);

                d2Ays0[n] = (bmu * ((-1 + n) * (1 + n - 2 * n3) * Inm1bJ * Math.cos((-1 + n) * t1) + (1 + n) * (1 - n + 2 * n3) * Inp1bJ * Math.cos(t1 + n * t1) + 2 * InbJ * (n * Math.cos(t1) * Math.cos(n * t1) + n * (-1 + 2 * n2) * Math.cos(t1) * Math.cos(n * t1) - Math.sin(t1) * Math.sin(n * t1) - n2 * (-1 + 2 * n2) * Math.sin(t1) * Math.sin(n * t1)))) / (1 + 4 * n4);


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


        double p0 = Math.exp(bJ * Math.cos(t1 - t2));
        double pM1 = 1 / p0;
        double pM2 = pM1 / p0;
        double lnp1 = -bJ * Math.sin(t1 - t2);
        double lnp11 = -bJ * Math.cos(t1 - t2);
        double lnp2 = -lnp1;
        double px1 = bmu * p0 * (Math.cos(t1) + Math.cos(t2));
        double py1 = bmu * p0 * (Math.sin(t1) + Math.sin(t2));
        psix1 = 0;
        psix2 = 0;
        psix11 = 0;
        psix12 = 0;
        psix22 = 0;
        psi1x1 = 0;
        psi1x2 = 0;
        psiy1 = 0;
        psiy2 = 0;
        psiy11 = 0;
        psiy12 = 0;
        psiy22 = 0;
        psi1y1 = 0;
        psi1y2 = 0;
        psi1x11 = 0;
        psi1y11 = 0;
        psi1x22 = 0;
        psi1y22 = 0;
        psix111 = 0;
        psiy111 = 0;
        psix112 = 0;
        psiy112 = 0;
        psix221 = 0;
        psiy221 = 0;
        psix222 = 0;
        psiy222 = 0;
        psi1x11 = 0;
        psi1y11 = 0;
        psi1x22 = 0;
        psi1y22 = 0;
        for (int n = 0; n <= nMax; n++) {
            psix1 += dAxs0[n] * Math.sin(n * t2) + dAxc0[n] * Math.cos(n * t2);
            psiy1 += dAys0[n] * Math.sin(n * t2) + dAyc0[n] * Math.cos(n * t2);

            psix2 += n * Axs0[n] * Math.cos(n * t2) - n * Axc0[n] * Math.sin(n * t2);
            psiy2 += n * Ays0[n] * Math.cos(n * t2) - n * Ayc0[n] * Math.sin(n * t2);

            psi1x1 += dAxs1[n] * Math.sin(n * t2) + dAxc1[n] * Math.cos(n * t2);
            psi1y1 += dAys1[n] * Math.sin(n * t2) + dAyc1[n] * Math.cos(n * t2);


            psi1x2 += n * Axs1[n] * Math.cos(n * t2) - n * Axc1[n] * Math.sin(n * t2);
            psi1y2 += n * Ays1[n] * Math.cos(n * t2) - n * Ayc1[n] * Math.sin(n * t2);

            psix11 += d2Axs0[n] * Math.sin(n * t2) + d2Axc0[n] * Math.cos(n * t2);
            psiy11 += d2Ays0[n] * Math.sin(n * t2) + d2Ayc0[n] * Math.cos(n * t2);

            psix112 += (n * d2Axs0[n] * Math.cos(n * t2)) - (n * d2Axc0[n] * Math.sin(n * t2));
            psiy112 += (n * d2Ays0[n] * Math.cos(n * t2)) - (n * d2Ayc0[n] * Math.sin(n * t2));

            psix111 += d3Axs0[n] * Math.sin(n * t2) + d3Axc0[n] * Math.cos(n * t2);
            psiy111 += d3Ays0[n] * Math.sin(n * t2) + d3Ayc0[n] * Math.cos(n * t2);

            psi1x11 += d2Axs1[n] * Math.sin(n * t2) + d2Axc1[n] * Math.cos(n * t2);
            psi1y11 += d2Ays1[n] * Math.sin(n * t2) + d2Ayc1[n] * Math.cos(n * t2);

            psix12 += n * dAxs0[n] * Math.cos(n * t2) - n * dAxc0[n] * Math.sin(n * t2);
            psiy12 += n * dAys0[n] * Math.cos(n * t2) - n * dAyc0[n] * Math.sin(n * t2);

            psix22 += -n * n * Axs0[n] * Math.sin(n * t2) - n * n * Axc0[n] * Math.cos(n * t2);
            psiy22 += -n * n * Ays0[n] * Math.sin(n * t2) - n * n * Ayc0[n] * Math.cos(n * t2);

            psix221 += -n * n * dAxs0[n] * Math.sin(n * t2) - n * n * dAxc0[n] * Math.cos(n * t2);
            psiy221 += -n * n * dAys0[n] * Math.sin(n * t2) - n * n * dAyc0[n] * Math.cos(n * t2);

            psix222 += -n * n * n * Axs0[n] * Math.cos(n * t2) + n * n * n * Axc0[n] * Math.sin(n * t2);
            psiy222 += -n * n * n * Ays0[n] * Math.cos(n * t2) + n * n * n * Ayc0[n] * Math.sin(n * t2);

            psi1x22 += -n * n * Axs1[n] * Math.sin(n * t2) - n * n * Axc1[n] * Math.cos(n * t2);
            psi1y22 += -n * n * Ays1[n] * Math.sin(n * t2) - n * n * Ayc1[n] * Math.cos(n * t2);
        }

//        System.out.println("psix1- " + "(" + psix1 + ")");
//        System.out.println("psiy1- " + "(" + psiy1 + ")");
//        System.out.println("psix2- " + "(" + psix2 + ")");
//        System.out.println("psiy2- " + "(" + psiy2 + ")");
//        System.out.printmln("psi1x1- " + "(" + psi1x1 + ")");
//        System.out.println("psi1y1- " + "(" + psi1y1 + ")");
//        System.out.println("psi1x2- " + "(" + psi1x2 + ")");
//        System.out.println("psi1y2- " + "(" + psi1y2 + ")");
//        System.out.println("psix11- " + "(" + psix11 + ")");
//        System.out.println("psiy11- " + "(" + psiy11 + ")");
//        System.out.println("psix22- " + "(" + psix22 + ")");
//        System.out.println("psiy22- " + "(" + psiy22 + ")");
//        System.exit(2);


        double vEx1 = psix1 / p0;
        double vEx2 = psix2 / p0;
        double vEy1 = psiy1 / p0;
        double vEy2 = psiy2 / p0;

        double vEEx1 = psi1x1 / p0 - px1 * vEx1 / p0;
        double vEEx2 = psi1x2 / p0 - px1 * vEx2 / p0;

        double vEEy1 = psi1y1 / p0 - py1 * vEy1 / p0;
        double vEEy2 = psi1y2 / p0 - py1 * vEy2 / p0;

//            System.out.println("vEx1[nMax]- " + "(" + vEx1 + ")");
//            System.out.println("vEx2[nMax]- " + "(" + vEx2 + ")");
//            System.out.println("vEx2[nMax]- " + "(" + vEx2 + ")");
//            System.out.println("vEy1[nMax]- " + "(" + vEy1 + ")");
//            System.out.println("vEy2[nMax]- " + "(" + vEy2 + ")");
//            System.out.println("vEEx1[nMax]- " + "(" + vEEx1 + ")");
//            System.out.println("vEEx2[nMax]- " + "(" + vEEx2 + ")");
//            System.out.println("vEEy1[nMax]- " + "(" + vEEy1 + ")");
//            System.out.println("vEEy2[nMax]- " + "(" + vEEy2 + ")");
//        System.exit(2);


//x
        JEEMJEJE += (((psi1x11 + psi1x22) / p0) + (((psi1x1 * bJ * Math.sin(t1 - t2)) + (psi1x2 * (-bJ * Math.sin(t1 - t2)))) / p0)) + (-bmu * ((-Math.sin(t1) * vEx1 - Math.sin(t2) * vEx2) + ((Math.cos(t1) + Math.cos(t2)) * ((bJ * Math.sin(t1 - t2) * (psix1 - psix2) / p0) + ((psix11 + psix22) / p0))))) + (vEx1 * ((psix11 * bJ * Math.sin(t1 - t2) / p0) + (psix111 / p0) + (bJ * ((Math.sin(t1 - t2) * psix11) + (psix1 * Math.cos(t1 - t2)) + (psix1 * bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2))) / p0))) + (vEx1 * ((psix22 * bJ * Math.sin(t1 - t2) / p0) + (psix221 / p0) + (psix12 * (-bJ * Math.sin(t1 - t2)) / p0) + (-psix2 * bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / p0))) + (vEx2 * ((psix11 * (-bJ * Math.sin(t1 - t2)) / p0) + (psix112 / p0) + (psix12 * bJ * Math.sin(t1 - t2) / p0) + (-psix1 * bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / p0))) + (vEx2 * ((psix22 * (-bJ * Math.sin(t1 - t2)) / p0) + (psix222 / p0) + (bJ * ((-Math.sin(t1 - t2) * psix22) + (psix2 * Math.cos(t1 - t2)) + (psix2 * bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2))) / p0)));

//y
        JEEMJEJE += (((psi1y11 + psi1y22) / p0) + (((psi1y1 * bJ * Math.sin(t1 - t2)) + (psi1y2 * (-bJ * Math.sin(t1 - t2)))) / p0)) + (-bmu * ((Math.cos(t1) * vEy1 + Math.cos(t2) * vEy2) + ((Math.sin(t1) + Math.sin(t2)) * ((bJ * Math.sin(t1 - t2) * (psiy1 - psiy2) / p0) + ((psiy11 + psiy22) / p0))))) + (vEy1 * ((psiy11 * bJ * Math.sin(t1 - t2) / p0) + (psiy111 / p0) + (bJ * ((Math.sin(t1 - t2) * psiy11) + (psiy1 * Math.cos(t1 - t2)) + (psiy1 * bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2))) / p0))) + (vEy1 * ((psiy22 * bJ * Math.sin(t1 - t2) / p0) + (psiy221 / p0) + (psiy12 * (-bJ * Math.sin(t1 - t2)) / p0) + (-psiy2 * bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / p0))) + (vEy2 * ((psiy11 * (-bJ * Math.sin(t1 - t2)) / p0) + (psiy112 / p0) + (psiy12 * bJ * Math.sin(t1 - t2) / p0) + (-psiy1 * bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / p0))) + (vEy2 * ((psiy22 * (-bJ * Math.sin(t1 - t2)) / p0) + (psiy222 / p0) + (bJ * ((-Math.sin(t1 - t2) * psiy22) + (psiy2 * Math.cos(t1 - t2)) + (psiy2 * bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2))) / p0)));


        double testx = (((psi1x11 + psi1x22) / p0) + (((psi1x1 * bJ * Math.sin(t1 - t2)) + (psi1x2 * (-bJ * Math.sin(t1 - t2)))) / p0)) + (-bmu * ((-Math.sin(t1) * vEx1 - Math.sin(t2) * vEx2) + ((Math.cos(t1) + Math.cos(t2)) * ((bJ * Math.sin(t1 - t2) * (psix1 - psix2) / p0) + ((psix11 + psix22) / p0))))) + (vEx1 * ((psix11 * bJ * Math.sin(t1 - t2) / p0) + (psix111 / p0) + (bJ * ((Math.sin(t1 - t2) * psix11) + (psix1 * Math.cos(t1 - t2)) + (psix1 * bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2))) / p0))) + (vEx1 * ((psix22 * bJ * Math.sin(t1 - t2) / p0) + (psix221 / p0) + (psix12 * (-bJ * Math.sin(t1 - t2)) / p0) + (-psix2 * bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / p0))) + (vEx2 * ((psix11 * (-bJ * Math.sin(t1 - t2)) / p0) + (psix112 / p0) + (psix12 * bJ * Math.sin(t1 - t2) / p0) + (-psix1 * bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / p0))) + (vEx2 * ((psix22 * (-bJ * Math.sin(t1 - t2)) / p0) + (psix222 / p0) + (bJ * ((-Math.sin(t1 - t2) * psix22) + (psix2 * Math.cos(t1 - t2)) + (psix2 * bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2))) / p0)));
        double testy = (((psi1y11 + psi1y22) / p0) + (((psi1y1 * bJ * Math.sin(t1 - t2)) + (psi1y2 * (-bJ * Math.sin(t1 - t2)))) / p0)) + (-bmu * ((Math.cos(t1) * vEy1 + Math.cos(t2) * vEy2) + ((Math.sin(t1) + Math.sin(t2)) * ((bJ * Math.sin(t1 - t2) * (psiy1 - psiy2) / p0) + ((psiy11 + psiy22) / p0))))) + (vEy1 * ((psiy11 * bJ * Math.sin(t1 - t2) / p0) + (psiy111 / p0) + (bJ * ((Math.sin(t1 - t2) * psiy11) + (psiy1 * Math.cos(t1 - t2)) + (psiy1 * bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2))) / p0))) + (vEy1 * ((psiy22 * bJ * Math.sin(t1 - t2) / p0) + (psiy221 / p0) + (psiy12 * (-bJ * Math.sin(t1 - t2)) / p0) + (-psiy2 * bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / p0))) + (vEy2 * ((psiy11 * (-bJ * Math.sin(t1 - t2)) / p0) + (psiy112 / p0) + (psiy12 * bJ * Math.sin(t1 - t2) / p0) + (-psiy1 * bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / p0))) + (vEy2 * ((psiy22 * (-bJ * Math.sin(t1 - t2)) / p0) + (psiy222 / p0) + (bJ * ((-Math.sin(t1 - t2) * psiy22) + (psiy2 * Math.cos(t1 - t2)) + (psiy2 * bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2))) / p0)));


        if (debug && false) {
            System.out.println("JEEMJEJEx[nMax]-(" + testx + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
            System.out.println("JEEMJEJEy[nMax]-(" + testy + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
        }
        //Var<JE-UE>
        MoleculeAgent agentAtom1 = leafAgentManager.getAgent(atom1);
        MoleculeAgent agentAtom2 = leafAgentManager.getAgent(atom2);
        double f1 = bt * agentAtom1.torque.getX(0);
        double f2 = bt * agentAtom2.torque.getX(0);


        //  -bJCos(t1-t2) for //ij
        double p12 = -bJ * ei.dot(ej);
        double p21 = p12;
        double p11 = bt * agentAtom1.phi.component(0, 0);
        double p22 = bt * agentAtom2.phi.component(0, 0);

        if (debug && false) {
            System.out.println("f1=" + f1 + ";f2=" + f2 + ";p12=" + (p12) + ";p11=" + p11 + ";p21=" + p21 + ";p22=" + p22 + ";");
//            System.out.println("f1-("+ f1+ ");");
//            System.out.println("total+=%^2;" );
//            System.out.println("f2-("+ f2+ ");");
//            System.out.println("total+=%^2;" );
//            System.out.println("p11-("+ p11+ ");");
//            System.out.println("total+=%^2;" );
//            System.out.println("p12-("+ p12+ ");");
//            System.out.println("total+=%^2;" );
//            System.out.println("p21-("+ p21+ ");");
//            System.out.println("total+=%^2;" );
//            System.out.println("p22-("+ p22+ ");");
//            System.out.println("total+=%^2;" );
        }


        //x
//        JEMUE += ((bJ * Math.sin(t1 - t2) * (psix1 - psix2) / p0) + ((psix11 + psix22) / p0)) - (-bt * (Math.cos(t1) + Math.cos(t2)) * mu - (vEx1 * (f1) + vEx2 * (f2)));
        JEMUEx += ((bJ * Math.sin(t1 - t2) * (psix1 - psix2) / p0) + ((psix11 + psix22) / p0)) - (-bt * (Math.cos(t1) + Math.cos(t2)) * mu - (vEx1 * (f1) + vEx2 * (f2)));
        double JEMUExNow = ((bJ * Math.sin(t1 - t2) * (psix1 - psix2) / p0) + ((psix11 + psix22) / p0)) - (-bt * (Math.cos(t1) + Math.cos(t2)) * mu - (vEx1 * (f1) + vEx2 * (f2)));

        JEMUExSquare += JEMUExNow * JEMUExNow;
        //y
//        JEMUE += ((bJ * Math.sin(t1 - t2) * (psiy1 - psiy2) / p0) + ((psiy11 + psiy22) / p0)) - (-bt * (Math.sin(t1) + Math.sin(t2)) * mu - (vEy1 * (f1) + vEy2 * (f2)));
        JEMUEy += ((bJ * Math.sin(t1 - t2) * (psiy1 - psiy2) / p0) + ((psiy11 + psiy22) / p0)) - (-bt * (Math.sin(t1) + Math.sin(t2)) * mu - (vEy1 * (f1) + vEy2 * (f2)));
        double JEMUEyNow = ((bJ * Math.sin(t1 - t2) * (psiy1 - psiy2) / p0) + ((psiy11 + psiy22) / p0)) - (-bt * (Math.sin(t1) + Math.sin(t2)) * mu - (vEy1 * (f1) + vEy2 * (f2)));

        JEMUEySquare += JEMUEyNow * JEMUEyNow;
        if (debug && false) {
            System.out.println("JEMUEx[nMax]-(" + JEMUExNow + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%^2;");
            System.out.println("JEMUEy[nMax]-(" + JEMUEyNow + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%^2;");
        }
        double vDotGradvx1 = pM2 * (-lnp1 * psix1 * (psix1 - psix2) + psix1 * psix11 + psix2 * psix12);
        double vDotGradvy1 = pM2 * (-lnp1 * psiy1 * (psiy1 - psiy2) + psiy1 * psiy11 + psiy2 * psiy12);

        double vDotGradvx2 = pM2 * (-lnp1 * psix2 * (psix1 - psix2) + psix1 * psix12 + psix2 * psix22);
        double vDotGradvy2 = pM2 * (-lnp1 * psiy2 * (psiy1 - psiy2) + psiy1 * psiy12 + psiy2 * psiy22);

        if (debug && false) {
            System.out.println("vDGvx1[theta1,theta2,nMax]-(" + vDotGradvx1 + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
            System.out.println("vDGvx2[theta1,theta2,nMax]-(" + vDotGradvx2 + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
            System.out.println("vDGvy1[theta1,theta2,nMax]-(" + vDotGradvy1 + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
            System.out.println("vDGvy2[theta1,theta2,nMax]-(" + vDotGradvy2 + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
        }

        double fxE1 = -bmu * Math.sin(t1);
        double fxE2 = -bmu * Math.sin(t2);
        double fyE1 = bmu * Math.cos(t1);
        double fyE2 = bmu * Math.cos(t2);
        UEE += -(vEEx1 * f1 + vEEx2 * f2 + vDotGradvx1 * f1 + vDotGradvx2 * f2)
                + vEx1 * (p11 * vEx1 + p12 * vEx2) + vEx2 * (p12 * vEx1 + p22 * vEx2)
                - 2 * vEx1 * fxE1 - 2 * vEx2 * fxE2;

        UEE += -(vEEy1 * f1 + vEEy2 * f2 + vDotGradvy1 * f1 + vDotGradvy2 * f2)
                + vEy1 * (p11 * vEy1 + p12 * vEy2) + vEy2 * (p12 * vEy1 + p22 * vEy2)
                - 2 * vEy1 * fyE1 - 2 * vEy2 * fyE2;


        if (debug && false) {
            testx = -(vEEx1 * f1 + vEEx2 * f2 + vDotGradvx1 * f1 + vDotGradvx2 * f2)
                    + vEx1 * (p11 * vEx1 + p12 * vEx2) + vEx2 * (p12 * vEx1 + p22 * vEx2)
                    - 2 * vEx1 * fxE1 - 2 * vEx2 * fxE2;
            testy = -(vEEy1 * f1 + vEEy2 * f2 + vDotGradvy1 * f1 + vDotGradvy2 * f2)
                    + vEy1 * (p11 * vEy1 + p12 * vEy2) + vEy2 * (p12 * vEy1 + p22 * vEy2)
                    - 2 * vEy1 * fyE1 - 2 * vEy2 * fyE2;
            System.out.println("f1=" + f1 + ";f2=" + f2 + ";p12=" + (p12) + ";p11=" + p11 + ";p21=" + p21 + ";p22=" + p22 + ";");
            System.out.println("UEEx[nMax]-(" + testx + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
            System.out.println("UEEy[nMax]-(" + testy + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
        }


//  System.out.println("theta1= " + t1 + " theta2= " + t2 + " dipolexmap= " + dxmap + " dipolexcon= " + dxcon + " JEMUExnow= " + JEMUExnow + " JEEMJEJExtrying=" + JEEMJEJExtrying + " UEEnow=" + UEEnow);
//System.out.println("theta1= " + t1 + " theta2= " + t2 + " UEExnow= " + UEExnow + " UEEynow= " + UEEynow + " JEEMJEJExtrying=" + JEEMJEJExtrying + " p11=" + p11 + " p12=" + p12 + " p22=" + p22 + " JEMUExnow= " + JEMUExnow + " JEMUEynow= " + JEMUEynow + " psix111 " + psix111+ " psix222 " + psix222+ " psix112 " + psix112+ " psix221 " + psix221+ " d2Axs0 " + d2Axs0[5]+ " d3Axs0 " + d3Axs0[5]+ " d2Axc0 " + d2Axc0[5]+ " d3Axc0 " + d3Axc0[5]+ " dAxc0 " + dAxc0[5]+ " psix11 " + psix11+ " psix22 " + psix22+ " psix12 " + psix12);
//        System.out.println("theta1= " + t1 + " theta2= " + t2 + "JEy=" + JEy + " UEEynow= " + UEEynow + " JEEMJEJEy=" + JEEMJEJEytrying + " p11=" + p11 + " p12=" + p12 + " p22=" + p22 + " JEMUExnow= " + JEMUExnow + " JEMUEynow= " + JEMUEynow + " psiy1 " + psiy1 + " psiy111 " + psiy111 + " psiy222 " + psiy222 + " psiy112 " + psiy112 + " psiy221 " + psiy221 + " psiy11 " + psiy11 + " psiy22 " + psiy22 + " psiy12 " + psiy12);

//        JEMUExIdeal += -bmu * (f1 * Math.sin(t1) + f2 * Math.sin(t2));
//        JEMUEyIdeal += bmu * (f1 * Math.cos(t1) + f2 * Math.cos(t2));
//        testx = -bmu * (f1 * Math.sin(t1) + f2 * Math.sin(t2));
//        testy = bmu * (f1 * Math.cos(t1) + f2 * Math.cos(t2));
//        JEMUEIdealSquare = testx * testx + testy * testy;

        if (debug && false) {
            testx = -bmu * (f1 * Math.sin(t1) + f2 * Math.sin(t2));
            testy = bmu * (f1 * Math.cos(t1) + f2 * Math.cos(t2));
            System.out.println("JEMUExJ0-(" + testx + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
            System.out.println("JEMUEyJ0-(" + testy + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
        }


        AEEJ0 += -bmu2 * (2 + f1 * f1 + f2 * f2 + 2 * f1 * f2 * Math.cos(t1 - t2) - p11 - p22 -
                2 * p12 * Math.cos(t1 - t2));
//        AEEJ0 += -2 * bmu2 + 4 * bJ * Math.sin((t1 - t2) / 2) * Math.sin((t1 - t2) / 2) * (
//                bmu2 * Math.cos(t1 - t2) - bJ * bmu2 * Math.sin(t1 - t2) * Math.sin(t1 - t2)  );

//        AEEJ0 += bmu2 * (-2 + p11 + p22 + (p12 + p21) * ei.dot(ej)) - JEMUEIdealSquare;


        if (debug && false) {
            double tmp = -2 * bmu2 + 4 * bJ * Math.sin((t1 - t2) / 2) * Math.sin((t1 - t2) / 2) * (
                    bmu2 * Math.cos(t1 - t2) - bJ * bmu2 * Math.sin(t1 - t2) * Math.sin(t1 - t2));
            System.out.println("AEEIdeal-(" + tmp + ")");
        }

        if (debug && false) {
            double testAEEJ0 = bmu2 * (-2 + p11 + p22 + (p12 + p21) * Math.cos(t1 - t2));
            System.out.println("AEEJ0-(" + testAEEJ0 + ")/.{theta1->" + t1 + " ,theta2->" + t2 + "};");
            System.out.println("total+=%;");
        }
    }


    public void zeroSum() {
        JEEMJEJE = 0;
        UEE = 0;
        JEMUEx = 0;
        JEMUExIdeal = 0;
        JEMUE = 0;//not used
        JEMUEy = 0;
        JEMUEyIdeal = 0;
        JEMUExSquare = 0;
        JEMUEySquare = 0;
        AEEJ0 = 0;
        JEMUEIdealSquare = 0;
    }

    public double getSumJEEMJEJE() {
        return JEEMJEJE;
    }

    public double getSumUEE() {
        return UEE;
    }

    public double getdipolex() {
        return dipolex;
    }

    public double getdipoley() {
        return dipoley;
    }

    public double getdipoleconv() {
        return dipoleconv;
    }

    public double getSumJEMUE() {
        return JEMUE;
    }

    public double getSumJEMUEx() {
        return JEMUEx;
    }


    public double getSumJEMUExIdeal() {
        return JEMUExIdeal;
    }

    public double getSumJEMUEy() {
        return JEMUEy;
    }

    public double getSumJEMUEyIdeal() {
        return JEMUEyIdeal;
    }

    public double getSumJEMUExSquare() {
        return JEMUExSquare;
    }

    public double getSumJEMUEySquare() {
        return JEMUEySquare;
    }

    public double getAEEJ0() {
        return AEEJ0;
    }

}
