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
 * Any eypressions for computing v_E and v_EE in the mapping.
 *
 * @author Weisong Lin
 */

public class PotentialCalculationPair implements PotentialCalculation {
    //public class PotentialCalculationHeisenberg {
    protected Vector ei, ej;
    protected double AEE, JEEMJEJE, UEE, JEMUEx, JEMUEy;
    protected final double mu, J, bt, bJ, bmu;
    protected double[] Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    protected double[] d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
    protected double[] Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
    protected double[] d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;
    protected double pvx10, pvx20, dpvx10dt1, dpvx20dt1, d2pvx10dt1dt1, d2pvx10dt1dt2, d2pvx20dt1dt2, d2pvx20dt2dt2, pvx11, pvx21, dpvx11dt1, dpvx11dt2, dpvx21dt1, dpvx21dt2, dpvx20dt2, d2pvx10dt1t1, dpvx10dt2;
    protected double pvy10, pvy20, dpvy10dt1, dpvy20dt1, d2pvy10dt1dt1, d2pvy10dt1dt2, d2pvy20dt1dt2, d2pvy20dt2dt2, pvy11, pvy21, dpvy11dt1, dpvy11dt2, dpvy21dt1, dpvy21dt2, dpvy20dt2, d2pvy10dt1t1, dpvy10dt2;

    protected Vector phi;
    protected Vector vEx, vEEx, dvEx, dvEEx, d2vEx;
    protected Vector vEy, vEEy, dvEy, dvEEy, d2vEy;
    protected double fEx, fEy, force;


    protected int nMax;
    protected int count = 1;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;

    public PotentialCalculationPair(Space space, double dipoleMagnitude, double interactionS, double beta, int nMax, AtomLeafAgentManager<MoleculeAgent> leafAgentManager) {
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
        IAtomOriented atom1 = (IAtomOriented) atoms.get(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.get(1);
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
        double sint2 = Math.sin(t2);
        double cost2 = Math.cos(t2);
        double sint1Square = sint1 * sint1;
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
        double px0 = p0, py0 = p0;
        double pM2 = 1 / p0 / p0;
        double lnp1 = -bJ * Math.sin(t1 - t2);
        double lnp11 = -bJ * Math.cos(t1 - t2);
        double lnp2 = -lnp1;
        double px1 = bmu * p0 * (Math.cos(t1) + Math.cos(t2));
        double py1 = bmu * p0 * (Math.sin(t1) + Math.sin(t2));
        double dRpx0dt1 = bJ * Math.sin(t1 - t2) / px0;//Rpx0 mean reverse of px0  D[1/px0, theta1]
        double dRpx0dt2 = -dRpx0dt1; // D[1/px0, theta2]
        double dpx1Rpx0dt1 = -bmu * sint1;//D[px1/px0, theta1]
        double dpx1Rpx0dt2 = -bmu * sint2;
        double d2Rpx0dt1dt1 = bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / px0;
        double d2Rpx0dt1dt2 = -d2Rpx0dt1dt1;
        double d2Rpx0dt2dt2 = d2Rpx0dt1dt1;

        double dRpy0dt1 = bJ * Math.sin(t1 - t2) / py0;//Rpy0 mean reverse of py0  D[1/py0, theta1]
        double dRpy0dt2 = -dRpy0dt1; // D[1/py0, theta2]
        double dpy1Rpy0dt1 = -bmu * cost1;//D[py1/py0, theta1]
        double dpy1Rpy0dt2 = -bmu * cost2;
        double d2Rpy0dt1dt1 = bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / py0;
        double d2Rpy0dt1dt2 = -d2Rpy0dt1dt1;
        double d2Rpy0dt2dt2 = d2Rpy0dt1dt1;


        pvx10 = 0;
        pvx20 = 0;
        dpvx10dt1 = 0;
        dpvx20dt1 = 0;
        d2pvx10dt1dt1 = 0;
        d2pvx10dt1dt2 = 0;
        d2pvx20dt1dt2 = 0;
        d2pvx20dt2dt2 = 0;
        pvx11 = 0;
        pvx21 = 0;
        dpvx11dt1 = 0;
        dpvx11dt2 = 0;
        dpvx21dt1 = 0;
        dpvx21dt2 = 0;
        pvy10 = 0;
        pvy20 = 0;
        dpvy10dt1 = 0;
        dpvy20dt1 = 0;
        d2pvy10dt1dt1 = 0;
        d2pvy10dt1dt2 = 0;
        d2pvy20dt1dt2 = 0;
        d2pvy20dt2dt2 = 0;
        pvy11 = 0;
        pvy21 = 0;
        dpvy11dt1 = 0;
        dpvy11dt2 = 0;
        dpvy21dt1 = 0;
        dpvy21dt2 = 0;


        for (int n = 0; n <= nMax; n++) {
            double sinnt2 = Math.sin(n * t2);
            double cosnt2 = Math.cos(n * t2);

            pvx10 += dAxs0[n] * sinnt2 + dAxc0[n] * cosnt2;
            pvy10 += dAys0[n] * sinnt2 + dAyc0[n] * cosnt2;

            pvx20 += n * Axs0[n] * cosnt2 - n * Axc0[n] * sinnt2;
            pvy20 += n * Ays0[n] * cosnt2 - n * Ayc0[n] * sinnt2;

            pvx11 += dAxs1[n] * sinnt2 + dAxc1[n] * cosnt2;
            pvy11 += dAys1[n] * sinnt2 + dAyc1[n] * cosnt2;

            pvx21 += n * Axs1[n] * cosnt2 - n * Axc1[n] * sinnt2;
            pvy21 += n * Ays1[n] * cosnt2 - n * Ayc1[n] * sinnt2;

            dpvx10dt1 += d2Axs0[n] * sinnt2 + d2Axc0[n] * cosnt2;
            dpvy10dt1 += d2Ays0[n] * sinnt2 + d2Ayc0[n] * cosnt2;

            d2pvx10dt1t1 += d3Axs0[n] * sinnt2 + d3Axc0[n] * cosnt2;
            d2pvy10dt1t1 += d3Ays0[n] * sinnt2 + d3Ayc0[n] * cosnt2;

            dpvx10dt2 += n * dAxs0[n] * cosnt2 - n * dAxc0[n] * sinnt2;
            dpvy10dt2 += n * dAys0[n] * cosnt2 - n * dAyc0[n] * sinnt2;

            d2pvx10dt1dt2 += n * d2Axs0[n] * cosnt2 - n * d2Axc0[n] * sinnt2;
            d2pvy10dt1dt2 += n * d2Ays0[n] * cosnt2 - n * d2Ayc0[n] * sinnt2;

            dpvx20dt1 += n * dAxs0[n] * cosnt2 - n * dAxc0[n] * sinnt2;
            dpvy20dt1 += n * dAys0[n] * cosnt2 - n * dAyc0[n] * sinnt2;

            d2pvx20dt1dt2 += -n * n * dAxs0[n] * sinnt2 - n * n * dAxc0[n] * cosnt2;
            d2pvy20dt1dt2 += -n * n * dAys0[n] * sinnt2 - n * n * dAyc0[n] * cosnt2;

            dpvx20dt2 += -n * n * Axs0[n] * sinnt2 - n * n * Axc0[n] * cosnt2;
            dpvy20dt2 += -n * n * Ays0[n] * sinnt2 - n * n * Ayc0[n] * cosnt2;

            d2pvx20dt1dt2 += -n * n * dAxs0[n] * sinnt2 - n * n * dAxc0[n] * cosnt2;
            d2pvy20dt1dt2 += -n * n * dAys0[n] * sinnt2 - n * n * dAyc0[n] * cosnt2;

            dpvx11dt1 += d2Axs1[n] * sinnt2 + d2Axc1[n] * cosnt2;
            dpvy11dt1 += d2Ays1[n] * sinnt2 + d2Ayc1[n] * cosnt2;

            dpvx11dt2 += n * dAxs1[n] * cosnt2 - n * dAxc1[n] * sinnt2;
            dpvy11dt2 += n * dAys1[n] * cosnt2 - n * dAyc1[n] * sinnt2;

            dpvx21dt1 += n * dAxs1[n] * cosnt2 - n * dAxc1[n] * sinnt2;
            dpvy21dt1 += n * dAys1[n] * cosnt2 - n * dAyc1[n] * sinnt2;

            dpvx21dt2 += -n * n * Axs1[n] * sinnt2 - n * n * Axc1[n] * cosnt2;
            dpvy21dt2 += -n * n * Ays1[n] * sinnt2 - n * n * Ayc1[n] * cosnt2;


        }

        MoleculeAgent agentAtom1 = leafAgentManager.getAgent(atom1);
        MoleculeAgent agentAtom2 = leafAgentManager.getAgent(atom2);


        double f1 = bt * agentAtom1.torque.getX(0);
        double f2 = bt * agentAtom2.torque.getX(0);
//        double p11 = bt * agentAtom1.phi.component(0, 0);
//        double p22 = bt * agentAtom2.phi.component(0, 0);
        double p12 = -bJ * ei.dot(ej);
        double p21 = p12;
//        double fxE1 = -bmu * Math.sin(t1);
//        double fxE2 = -bmu * Math.sin(t2);
//        double fyE1 = bmu * Math.cos(t1);
//        double fyE2 = bmu * Math.cos(t2);


        // JEEMJEJE =  sum(dvEEi/dti) + sum_i sim_j (vEi*d2vEj/dtidtj )
//        double dvEEx1 = agentAtom1.dvEEx.getX(0);
//        double dvEEy1 = agentAtom1.dvEEy.getX(0);
//        double dvEEx2 = agentAtom2.dvEEx.getX(0);
//        double dvEEy2 = agentAtom2.dvEEy.getX(0);

        //sum(dvEEi/dti)
//        JEEMJEJE += dvEEx1 + dvEEy1 + dvEEx2 + dvEEy2;

        double vEx1 = agentAtom1.vEx.getX(0);
        double vEy1 = agentAtom1.vEx().getX(0);
        double d2vEx2dt1dt2 = d2pvx20dt1dt2 / px0 + dpvx20dt2 * dRpx0dt1 + dpvx20dt1 * dRpx0dt2 + pvx20 * d2Rpx0dt1dt2;
        double d2vEy2dt1dt2 = d2pvy20dt1dt2 / py0 + dpvy20dt2 * dRpy0dt1 + dpvy20dt1 * dRpy0dt2 + pvy20 * d2Rpy0dt1dt2;
        double vEx2 = agentAtom1.vEx.getX(0);
        double vEy2 = agentAtom1.vEx().getX(0);
        double d2vEx1dt1dt2 = d2pvx10dt1dt2 / px0 + dpvx10dt1 * dRpx0dt2 + dpvx10dt2 * dRpx0dt1 + pvx10 * d2Rpx0dt1dt2;
        double d2vEy1dt1dt2 = d2pvy10dt1dt2 / py0 + dpvy10dt1 * dRpy0dt2 + dpvy10dt2 * dRpy0dt1 + pvy10 * d2Rpy0dt1dt2;

        //sum_i sim_j (vEi*d2vEj/dtidtj ) i!=j case
        JEEMJEJE += vEx1 * d2vEx2dt1dt2 + vEy1 * d2vEy2dt1dt2 + vEx2 * d2vEx1dt1dt2 + vEy2 * d2vEy1dt1dt2;


        //sum_i sim_j (vEi*d2vEj/dtidtj ) i==j case
//        double d2vEx1dt1dt1 = agentAtom1.d2vEx().getX(0);
//        double d2vEy1dt1dt1 = agentAtom1.d2vEy().getX(0);
//        double d2vEx2dt2dt2 = agentAtom2.d2vEx().getX(0);
//        double d2vEy2dt2dt2 = agentAtom2.d2vEy().getX(0);
//        JEEMJEJE += vEx1 * d2vEx1dt1dt1 + vEy1 * d2vEy1dt1dt1 + vEx2 * d2vEx2dt2dt2 + vEy2 * d2vEy2dt2dt2;


        //UEE = -sum_i(vEEi*fi) - sum_i sum_j(fi*dvEi/dtj*vEj) +  sum_i sum_j(vEi phiij vEj) - 2 sum_i(vEi*fEi)
//        double vEEx1 = agentAtom1.vEEx.getX(0);
//        double vEEy1 = agentAtom1.vEEy.getX(0);
//        double vEEx2 = agentAtom2.vEEx.getX(0);
//        double vEEy2 = agentAtom2.vEEy.getX(0);
        //-sum_i(vEEi*fi)
//        UEE += -vEEx1 * f1 - vEEy1 * f1 - vEEx2 * f2 - vEEy2 * f2;


        //- sum_i sum_j(fi*dvEi/dtj*vEj) i!=j case
        double dvEx1dt2 = dpvx10dt2 / px0 + pvx10 * dRpx0dt2;
        double dvEy1dt2 = dpvy10dt2 / py0 + pvy10 * dRpy0dt2;
        double dvEx2dt1 = dpvx20dt1 / px0 + pvx20 * dRpx0dt1;
        double dvEy2dt1 = dpvy20dt1 / py0 + pvy20 * dRpy0dt1;
        UEE += -f1 * dvEx1dt2 * vEx2 - f2 * dvEx2dt1 * vEx1 - -f1 * dvEy1dt2 * vEy2 - f2 * dvEy2dt1 * vEy1;

        //- sum_i sum_j(fi*dvEi/dtj*vEj) i==j case
//        double dvEx1dt1 = agentAtom1.dvEx().getX(0);
//        double dvEy1dt1 = agentAtom1.dvEy().getX(0);
//        double dvEx2dt2 = agentAtom2.dvEx().getX(0);
//        double dvEy2dt2 = agentAtom2.dvEy().getX(0);
//        UEE += -f1 * dvEx1dt1 * vEx1 - f2 * dvEx2dt2 * vEx2 - f1 * dvEy1dt1 * vEy1 - f2 * dvEy2dt2 * vEy2;

        //sum_i sum_j(vEi phiij vEj) i!=j case
        UEE += vEx1 * p12 * vEx2 + vEx2 * p21 * vEx1 + vEy1 * p12 * vEy2 + vEy2 * p21 * vEy1;

        //sum_i sum_j(vEi phiij vEj) i==j case
//        UEE += vEx1 * p11 * vEx1 + vEx2 * p22 * vEx2 + vEy1 * p11 * vEy1 + vEy2 * p22 * vEy2;

        //- 2 sum_i(vEi*fEi)
//        UEE += -2 * (vEx1 * fxE1 + vEx2 * fxE2 + vEy1 * fyE1 + vEy2 * fyE2);

        //JEMUEx = sum_i( dvEx + bmu cos(theta) + vE*f)
//        JEMUEx += dvEx1dt1 + bmu * cost1 + vEx1 * f1 + dvEx2dt2 + bmu * cost2 + vEx2 * f2;


        //JEMUEy = sum_i( dvEy + bmu sin(theta) + vE*f)
//        JEMUEy += dvEy1dt1 + bmu * sint1 + vEy1 * f1 + dvEy2dt2 + bmu * sint2 + vEy2 * f2;

    }

    public void zeroSum() {
        JEEMJEJE = 0;
        UEE = 0;
        AEE = 0;
        JEMUEx = 0;
        JEMUEy = 0;
    }

    public double getSumJEEMJEJE() {
        return JEEMJEJE;
    }

    public double getSumUEE() {
        return UEE;
    }

    public double getSumAEE() {
        return AEE;
    }

    public double getSumJEMUEx() {
        return JEMUEx;
    }

    public double getSumJEMUEy() {
        return JEMUEy;
    }


}
