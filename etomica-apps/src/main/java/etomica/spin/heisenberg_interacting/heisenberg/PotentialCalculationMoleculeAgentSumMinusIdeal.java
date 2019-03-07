package etomica.spin.heisenberg_interacting.heisenberg;

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

public class PotentialCalculationMoleculeAgentSumMinusIdeal implements PotentialCalculation {
    //public class PotentialCalculationHeisenberg {
    protected AtomLeafAgentManager.AgentIterator leafAgentIterator;
    protected Vector ei, ej;
    protected double AEEJ0, JEMUExIdeal, JEMUEyIdeal, JEMUEIdealSquare, JEEMJEJE, UEE, JEMUExSquare, JEMUEySquare, JEMUEx, JEMUEy, dipolex, dipoley, JEEMJEJExtrying, UEEnow, JEMUE, dipoleconv;
    protected final double mu, J, bt, bJ, bmu, bmu2, I0bJ, I1bJ, I2bJ;
    protected double[] Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    protected double[] d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
    protected double[] Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
    protected double[] d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;
    protected double[] InbJArray, Inm1bJArray, Inm2bJArray, Inp1bJArray;
    protected double pvx10, pvx20, dpvx10dt1, dpvx20dt1, d2pvx10dt1dt2, d2pvx20dt1dt2, d2pvx20dt2dt2, pvx11, pvx21, dpvx11dt1, dpvx11dt2, dpvx21dt1, dpvx21dt2, dpvx20dt2, d2pvx10dt1dt1, dpvx10dt2;
    protected double pvy10, pvy20, dpvy10dt1, dpvy20dt1, d2pvy10dt1dt2, d2pvy20dt1dt2, d2pvy20dt2dt2, pvy11, pvy21, dpvy11dt1, dpvy11dt2, dpvy21dt1, dpvy21dt2, dpvy20dt2, d2pvy10dt1dt1, dpvy10dt2;

    protected int nMax;
    protected int count = 1;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;

    public PotentialCalculationMoleculeAgentSumMinusIdeal(Space space, double dipoleMagnitude, double interactionS, double beta, int nMax, AtomLeafAgentManager<MoleculeAgent> leafAgentManager) {
        ei = space.makeVector();
        ej = space.makeVector();

        J = interactionS;
        mu = dipoleMagnitude;
        bt = beta;
        bJ = bt * J;
        bmu = bt * mu;
        bmu2 = bmu * bmu;
        this.nMax = nMax;
        this.leafAgentManager = leafAgentManager;

        leafAgentIterator = leafAgentManager.makeIterator();

        I0bJ = besselI(0, bJ);
        I1bJ = besselI(1, bJ);
        I2bJ = besselI(2, bJ);
        InbJArray = new double[nMax + 1];
        Inm1bJArray = new double[nMax + 1];
        Inm2bJArray = new double[nMax + 1];
        Inp1bJArray = new double[nMax + 1];
        for (int n = 1; n <= nMax; n++) {
            InbJArray[n] = besselI(n, bJ);
            Inm1bJArray[n] = besselI(n - 1, bJ);
            Inm2bJArray[n] = besselI(n - 2, bJ);
            Inp1bJArray[n] = besselI(n + 1, bJ);
        }
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
        boolean debug = false;
        if (debug) {
//            System.out.println("nMax= " + nMax + ";");
//            System.out.println("mu= " + mu + ";");
//            System.out.println("J= " + J + ";");
//            System.out.println("bJ= " + bJ + ";");
//            System.out.println("bt = bJ/J;");
//            System.out.println("bmu = bt*mu;");
//            System.out.println("theta1 = " + t1 + ";");
//            System.out.println("theta2 = " + t2 + ";");
        }

        if (count > 20 && debug) System.exit(2);


//        System.out.println("ei={" + ei.getX(0) + "," + ei.getX(1) + "};");
//        System.out.println("ej={" + ej.getX(0) + "," + ej.getX(1) + "};");
//        System.out.println("theta1 = " + t1 + ";");
//        System.out.println("theta2 = " + t2 + ";");
//        System.out.println("t1 = " + t1 + " c1=" + c1);
//System.out.println("Acos eix = " + Math.acos(ei.getX(0)));

        double cost1 = ei.getX(0);
        double sint1 = ei.getX(1);
        double sint2 = ej.getX(1);
        double cost2 = ej.getX(0);
        double cos2t1 = 2 * cost1 * cost1 - 1;
        double cos2t2 = 2 * cost2 * cost2 - 1;
        double cost1mt2 = cost1 * cost2 + sint1 * sint2;
        double sint1mt2 = sint1 * cost2 - cost1 * sint2;
        double sin2t1 = 2 * sint1 * cost1;

        for (int n = 0; n <= nMax; n++) {
            if (n == 0) {
                Axc0[n] = bmu * (I0bJ + I1bJ) * (-1 + cost1);
                dAxc0[n] = -bmu * (I0bJ + I1bJ) * sint1;
                Axs0[n] = 0;
                dAxs0[n] = 0;
                Axc1[n] = -0.25 * bmu2 * sint1 * sint1 * (I0bJ + 2 * I1bJ + I2bJ);
                dAxc1[n] = -0.25 * bmu2 * (I0bJ + 2 * I1bJ + I2bJ) * sin2t1;
                Axs1[n] = 0;
                dAxs1[n] = 0;

                d2Axc0[n] = -bmu * cost1 * (I0bJ + I1bJ);
                d2Axs0[n] = 0;
                d3Axc0[n] = bmu * sint1 * (I0bJ + I1bJ);
                d3Axs0[n] = 0;
                d2Axc1[n] = -0.5 * bmu2 * cos2t1 * (I0bJ + 2 * I1bJ + I2bJ);
                d2Axs1[n] = 0;


                //y direction
                Ayc0[n] = bmu * (I0bJ + I1bJ) * sint1;
                dAyc0[n] = bmu * (I0bJ + I1bJ) * cost1;
                Ays0[n] = 0;
                dAys0[n] = 0;
                Ayc1[n] = 0.25 * bmu2 * sint1 * sint1 * (I0bJ + 2 * I1bJ + I2bJ);
                dAyc1[n] = 0.25 * bmu2 * sin2t1 * (I0bJ + 2 * I1bJ + I2bJ);
                Ays1[n] = 0;
                dAys1[n] = 0;

                d2Ayc0[n] = -bmu * sint1 * (I0bJ + I1bJ);
                d2Ays0[n] = 0;
                d3Ayc0[n] = -bmu * cost1 * (I0bJ + I1bJ);
                d3Ays0[n] = 0;
                d2Ayc1[n] = 0.5 * bmu2 * (I0bJ + 2 * I1bJ + I2bJ) * cos2t1;
                d2Ays1[n] = 0;
            }


            if (n > 0) {
                int n2 = n * n;
                int n3 = n2 * n;
                int n4 = n2 * n2;
                double InbJ = InbJArray[n];
                double Inm1bJ = Inm1bJArray[n];
                double Inm2bJ = Inm2bJArray[n];
                double Inp1bJ = Inp1bJArray[n];
                double sinnt1 = Math.sin(n * t1);
                double cosnt1 = Math.cos(n * t1);
                double sinnm1t1 = sinnt1 * cost1 - cosnt1 * sint1;
                double sinnp1t1 = sinnt1 * cost1 + cosnt1 * sint1;
                double coshnt1 = Math.cosh(n * t1);
                double sinnm2t1 = sinnt1 * cos2t1 - cosnt1 * sin2t1;
                double sinnp2t1 = sinnt1 * cos2t1 + cosnt1 * sin2t1;
                double cosnm2t1 = cosnt1 * cos2t1 + sinnt1 * sin2t1;
                double cosnp2t1 = cosnt1 * cos2t1 - sinnt1 * sin2t1;
                double cosnm1t1 = cosnt1 * cost1 + sinnt1 * sint1;
                double cosnp1t1 = cosnt1 * cost1 - sinnt1 * sint1;


                Axc0[n] = 2 * bmu * (((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * (1 + bJ) * n2 - 2 * n3) * InbJ) * cost1 * cosnt1
                        + (2 * bJ * Inm1bJ + (1 + 2 * bJ - 2 * n + 2 * n2) * InbJ) * n * sint1 * sinnt1)
                        / (bJ + 4 * bJ * n4);

                dAxc0[n] = (2 * InbJ * (-cosnt1 * sint1 + (n - 2 * n3) * cost1 * sinnt1)
                        + Inm1bJ * (1 + n - 2 * n3) * sinnm1t1
                        + (-1 + n - 2 * n3) * Inp1bJ * sinnp1t1
                ) * bmu / (1 + 4 * n4);

                Axs0[n] = ((1 + 2 * n + 2 * n2) * Inm1bJ * sinnm1t1
                        + 2 * InbJ * (-2 * n * cosnt1 * sint1 + (1 + 2 * n2) * cost1 * sinnt1)
                        + (1 - 2 * n + 2 * n2) * Inp1bJ * sinnp1t1
                ) * bmu / (1 + 4 * n4);

                dAxs0[n] = (-2 * bmu * n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * cost1 * cosnt1
                        - 2 * bmu * (bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * sint1 * sinnt1
                ) / (bJ + 4 * bJ * n4);

                Axc1[n] = (Inm2bJ * (cosnm2t1 - coshnt1) * n2 / (2 - 2 * n + n2)
                        + (-4 * bJ * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 + n2 - 2 * n) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (cosnm2t1 + coshnt1) + 2 * Inm1bJ * (bJ * cosnm2t1 + (n - 1) * coshnt1))
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu2 / 4.0 / n2;

                dAxc1[n] = 0.25 * bmu2 * (
                        -(n - 2) * Inm2bJ * sinnm2t1 / (2 - 2 * n + n2)
                                +
                                (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                                        - 2 * n * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (bJ - 1 - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)

                                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                Axs1[n] = ((n2 * Inm2bJ * sinnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n2 * (2 + 2 * n + n2) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) - (16 + 4 * n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n2 * (2 - 2 * n + n2) * I0bJ * sinnp2t1 * ((-bJ + bJ * bJ - n * bJ) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n2 + 2 * n) * InbJ)
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu2 / 4.0 / n2;

                dAxs1[n] = 0.25 * bmu2 * (((n - 2) * Inm2bJ * cosnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ)
                        - 4 * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                d2Axc0[n] = (-(1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * cosnp1t1
                        + (-1 + n) * Inm1bJ * ((1 + n) * cosnm1t1 - 2 * n3 * cosnm1t1)
                        + InbJ * (-2 * (1 - n2 + 2 * n4) * cost1 * cosnt1 + 4 * n3 * sint1 * sinnt1)
                ) * bmu / (1 + 4 * n4);

                d2Axs0[n] = (-2 * bJ * bmu * Inm1bJ * (2 * n3 * cosnt1 * sint1 + (1 - n2 + 2 * n4) * cost1 * sinnt1)
                        + 2 * bmu * InbJ * (n * (1 - n2 - 2 * bJ * n2 + 2 * n3 + 2 * n4) * cosnt1 * sint1 + (n * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 2 * n4)) * cost1 * sinnt1)
                ) / (bJ + 4 * bJ * n4);

                d3Axc0[n] = (2 * InbJ * ((1 - n2 + 4 * n4) * cosnt1 * sint1 + n * (1 + n2 + 2 * n4) * cost1 * sinnt1)
                        + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * sinnp1t1
                        - (-1 + n) * (-1 + n) * Inm1bJ * ((1 + n) * sinnm1t1 + 2 * n3 * sinnm1t1)
                ) * bmu / (1.0 + 4 * n4);

                d3Axs0[n] = -(bJ * Inm1bJ * (n * (1 + n2 + 2 * n4) * cost1 * cosnt1 + (-1 + n2 - 4 * n4) * sint1 * sinnt1)
                        + InbJ * (-n * cost1 * cosnt1 * ((1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) - bJ * (1 + n2 + 2 * n4)) + sint1 * sinnt1 * (n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (-1 + n2 - 4 * n4)))
                ) * 2 * bmu / (bJ + 4 * bJ * n4);

                d2Axc1[n] = 0.25 * bmu2 * (-(n - 2) * (n - 2) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                d2Axs1[n] = 0.25 * bmu2 * (-((n - 2) * (n - 2) * Inm2bJ * sinnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (n + 2) * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                Ayc0[n] = (2 * bmu * cosnt1 * sint1 * ((bJ + 2 * bJ * n2) * Inm1bJ + (bJ - n + 2 * n2 + 2 * bJ * n2 - 2 * n3) * InbJ)
                        - 2 * bmu * n * cost1 * sinnt1 * (2 * bJ * Inm1bJ + (1 + 2 * bJ + 2 * (-1 + n) * n) * InbJ)
                ) / (bJ + 4 * bJ * n4);

                dAyc0[n] = 2.0 * bmu * ((bJ * Inm1bJ + (bJ - n + n2 - 2 * n4) * InbJ) * cost1 * cosnt1
                        + n * sint1 * sinnt1 * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (n - bJ) * InbJ))
                ) / (bJ + 4 * bJ * n4);

                Ays0[n] = (bJ * (1 - 2 * n + 2 * n2) * Inp1bJ * (-cosnp1t1 + coshnt1)
                        + bJ * Inm1bJ * ((1 + 2 * n + 2 * n2) * cosnm1t1 + (-1 + 2 * n - 2 * n2) * coshnt1)
                        + 2 * InbJ * (2 * bJ * n * cost1 * cosnt1 + n * (1 - 2 * n + 2 * n2) * coshnt1 + bJ * (1 + 2 * n2) * sint1 * sinnt1)
                ) * bmu / (bJ + 4 * bJ * n4);

                dAys0[n] = ((1 + n - 2 * n3) * Inm1bJ * sinnm1t1
                        + 2 * InbJ * (n * (-1 + 2 * n2) * cosnt1 * sint1 + cost1 * sinnt1)
                        + (1 - n + 2 * n3) * Inp1bJ * sinnp1t1
                ) * bmu / (1.0 + 4 * n4);


                Ayc1[n] = (n2 * Inm2bJ * (-cosnm2t1 + coshnt1) / (2 - 2 * n + n2)
                        - (4 * bJ * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        + 2 * n2 * (2 - 2 * n + n2) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                        + bJ * n2 * (2 + 2 * n + n2) * I0bJ * (bJ * InbJ * (cosnm2t1 + coshnt1) + 2 * Inm1bJ * (bJ * cosnm2t1 + (-1 + n) * coshnt1))
                ) / (bJ * bJ * (4 + n4) * I0bJ)
                ) * bmu2 / 4 / n2;

                dAyc1[n] = 0.25 * bmu2 * (((-2 + n) * Inm2bJ * sinnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                Ays1[n] = (-bJ * n2 * (2 + n * (2 + n)) * I0bJ * sinnm2t1 * ((-1 + bJ + n) * Inm1bJ + bJ * InbJ)
                        - 2 * bJ * (4 + n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        - n2 * (2 - 2 * n + n2) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n * (1 + n)) * InbJ)
                ) * bmu2 / (2 * bJ * bJ * n2 * (4 + n4) * I0bJ);

                dAys1[n] = 0.25 * bmu2 * (-(-2 + n) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (-bJ * n * (-4 - 2 * n + n3) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ) - 4 * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                d2Ayc0[n] = 2 * bmu * (-(bJ * Inm1bJ + (bJ + n * (-1 + n - 2 * n3)) * InbJ) * cosnt1 * sint1 + n2 * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * cosnt1 * sint1 - n * (bJ * Inm1bJ + (bJ + n * (-1 + n - 2 * n3)) * InbJ) * cost1 * sinnt1 + n * (InbJ + (-1 + 2 * n2) * (-bJ * Inm1bJ + (-bJ + n) * InbJ)) * cost1 * sinnt1) / (bJ + 4 * bJ * n4);

                d2Ays0[n] = (bmu * ((-1 + n) * (1 + n - 2 * n3) * Inm1bJ * cosnm1t1 + (1 + n) * (1 - n + 2 * n3) * Inp1bJ * cosnp1t1 + 2 * InbJ * (n * cost1 * cosnt1 + n * (-1 + 2 * n2) * cost1 * cosnt1 - sint1 * sinnt1 - n2 * (-1 + 2 * n2) * sint1 * sinnt1))) / (1 + 4 * n4);


                d3Ayc0[n] = 2 * bmu * (bJ * Inm1bJ * (-(1 - n2 + 4 * n4) * cost1 * cosnt1 + n * (1 + n2 + 2 * n4) * sint1 * sinnt1)
                        + InbJ * (-(-n * (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (1 - n2 + 4 * n4)) * cost1 * cosnt1 + n * (-(1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) + bJ * (1 + n2 + 2 * n4)) * sint1 * sinnt1)
                ) / (bJ + 4 * bJ * n4);


                d3Ays0[n] = -bmu * (-(-1 + n) * (-1 + n) * (-1 + n) * (1 + 2 * n + 2 * n2) * Inm1bJ * sinnm1t1
                        + 2 * InbJ * (n * (1 + n2 + 2 * n4) * cosnt1 * sint1 + (1 - n2 + 4 * n4) * cost1 * sinnt1)
                        + (1 + n) * (1 + n) * (1 + n) * (1 - 2 * n + 2 * n2) * Inp1bJ * sinnp1t1
                ) / (1 + 4 * n4);

                d2Ayc1[n] = 0.25 * bmu2 * ((-2 + n) * (-2 + n) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (bJ * (-2 + n) * (-2 + n) * n * (2 + 2 * n + n2) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                d2Ays1[n] = 0.25 * bmu2 * ((n - 2) * (n - 2) * Inm2bJ * sinnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


            }
        }


        double p0 = Math.exp(bJ * cost1mt2);
        double px0 = p0, py0 = p0;
        double px1 = bmu * p0 * (cost1 + cost2);
        double py1 = bmu * p0 * (sint1 + sint2);
        double dRpx0dt1 = bJ * sint1mt2 / px0;//Rpx0 mean reverse of px0  D[1/px0, theta1]
        double dRpx0dt2 = -dRpx0dt1; // D[1/px0, theta2]
        double dpx1Rpx0dt1 = -bmu * sint1;//D[px1/px0, theta1]
        double dpx1Rpx0dt2 = -bmu * sint2;
        double d2Rpx0dt1dt1 = bJ * (cost1mt2 + bJ * sint1mt2 * sint1mt2) / px0;
        double d2Rpx0dt2dt2 = d2Rpx0dt1dt1;

        double dRpy0dt1 = bJ * sint1mt2 / py0;//Rpy0 mean reverse of py0  D[1/py0, theta1]
        double dRpy0dt2 = -dRpy0dt1; // D[1/py0, theta2]
        double dpy1Rpy0dt1 = bmu * cost1;//D[py1/py0, theta1]
        double dpy1Rpy0dt2 = bmu * cost2;
        double d2Rpy0dt1dt1 = bJ * (cost1mt2 + bJ * sint1mt2 * sint1mt2) / py0;
        double d2Rpy0dt2dt2 = d2Rpy0dt1dt1;


        pvx10 = 0;
        pvy10 = 0;
        pvx11 = 0;
        pvy11 = 0;
        dpvx10dt1 = 0;
        dpvy10dt1 = 0;
        dpvx10dt2 = 0;
        dpvy10dt2 = 0;
        d2pvx10dt1dt2 = 0;
        d2pvy10dt1dt2 = 0;
        dpvx11dt1 = 0;
        dpvy11dt1 = 0;
        dpvx11dt2 = 0;
        dpvy11dt2 = 0;
        pvx20 = 0;
        pvy20 = 0;
        pvx21 = 0;
        pvy21 = 0;
        dpvx20dt1 = 0;
        dpvy20dt1 = 0;
        d2pvx20dt1dt2 = 0;
        dpvx20dt2 = 0;
        dpvy20dt2 = 0;
        d2pvy20dt1dt2 = 0;
        dpvx21dt1 = 0;
        dpvy21dt1 = 0;
        dpvx21dt2 = 0;
        dpvy21dt2 = 0;
        d2pvx10dt1dt1 = 0;
        d2pvy10dt1dt1 = 0;
        d2pvx20dt2dt2 = 0;
        d2pvy20dt2dt2 = 0;


        for (int n = 0; n <= nMax; n++) {
            double sinnt2 = Math.sin(n * t2);
            double cosnt2 = Math.cos(n * t2);
            double n2 = n * n;


            pvx10 += dAxs0[n] * sinnt2 + dAxc0[n] * cosnt2;
            pvy10 += dAys0[n] * sinnt2 + dAyc0[n] * cosnt2;

            pvx11 += dAxs1[n] * sinnt2 + dAxc1[n] * cosnt2;
            pvy11 += dAys1[n] * sinnt2 + dAyc1[n] * cosnt2;


            dpvx10dt1 += d2Axs0[n] * sinnt2 + d2Axc0[n] * cosnt2;
            dpvy10dt1 += d2Ays0[n] * sinnt2 + d2Ayc0[n] * cosnt2;

            d2pvx10dt1dt1 += d3Axs0[n] * sinnt2 + d3Axc0[n] * cosnt2;
            d2pvy10dt1dt1 += d3Ays0[n] * sinnt2 + d3Ayc0[n] * cosnt2;


            dpvx10dt2 += n * dAxs0[n] * cosnt2 - n * dAxc0[n] * sinnt2;
            dpvy10dt2 += n * dAys0[n] * cosnt2 - n * dAyc0[n] * sinnt2;

            d2pvx10dt1dt2 += n * d2Axs0[n] * cosnt2 - n * d2Axc0[n] * sinnt2;
            d2pvy10dt1dt2 += n * d2Ays0[n] * cosnt2 - n * d2Ayc0[n] * sinnt2;


            dpvx11dt1 += d2Axs1[n] * sinnt2 + d2Axc1[n] * cosnt2;
            dpvy11dt1 += d2Ays1[n] * sinnt2 + d2Ayc1[n] * cosnt2;

            dpvx11dt2 += n * dAxs1[n] * cosnt2 - n * dAxc1[n] * sinnt2;
            dpvy11dt2 += n * dAys1[n] * cosnt2 - n * dAyc1[n] * sinnt2;

            if (n != 0) {
                pvx20 += n * Axs0[n] * cosnt2 - n * Axc0[n] * sinnt2;
                pvy20 += n * Ays0[n] * cosnt2 - n * Ayc0[n] * sinnt2;

                pvx21 += n * Axs1[n] * cosnt2 - n * Axc1[n] * sinnt2;
                pvy21 += n * Ays1[n] * cosnt2 - n * Ayc1[n] * sinnt2;

                dpvx20dt1 += n * dAxs0[n] * cosnt2 - n * dAxc0[n] * sinnt2;
                dpvy20dt1 += n * dAys0[n] * cosnt2 - n * dAyc0[n] * sinnt2;

                d2pvx20dt1dt2 += -n2 * dAxs0[n] * sinnt2 - n2 * dAxc0[n] * cosnt2;
                d2pvy20dt1dt2 += -n2 * dAys0[n] * sinnt2 - n2 * dAyc0[n] * cosnt2;

                dpvx20dt2 += -n2 * Axs0[n] * sinnt2 - n2 * Axc0[n] * cosnt2;
                dpvy20dt2 += -n2 * Ays0[n] * sinnt2 - n2 * Ayc0[n] * cosnt2;

                d2pvx20dt2dt2 += -n2 * n * Axs0[n] * cosnt2 + n2 * n * Axc0[n] * sinnt2;
                d2pvy20dt2dt2 += -n2 * n * Ays0[n] * cosnt2 + n2 * n * Ayc0[n] * sinnt2;

                dpvx21dt1 += n * dAxs1[n] * cosnt2 - n * dAxc1[n] * sinnt2;
                dpvy21dt1 += n * dAys1[n] * cosnt2 - n * dAyc1[n] * sinnt2;

                dpvx21dt2 += -n2 * Axs1[n] * sinnt2 - n2 * Axc1[n] * cosnt2;
                dpvy21dt2 += -n2 * Ays1[n] * sinnt2 - n2 * Ayc1[n] * cosnt2;
            }
        }


        MoleculeAgent agentAtom1 = leafAgentManager.getAgent(atom1);
        MoleculeAgent agentAtom2 = leafAgentManager.getAgent(atom2);

        double vEx1 = pvx10 / p0;
        double vEx2 = pvx20 / p0;
        double vEy1 = pvy10 / p0;
        double vEy2 = pvy20 / p0;

//        double vEx1Ideal = -bmu * sint1;
//        double vEx2Ideal = -bmu * sint2;
//        double vEy1Ideal = bmu * cost1;
//        double vEy2Ideal = bmu * cost2;

//        System.out.println("vEx1Ideal  -(" + vEx1Ideal + ")");
//        System.out.println("vEx2Ideal -(" + vEx2Ideal + ")");
//        System.out.println("vEy1Ideal -(" + vEy1Ideal + ")");
//        System.out.println("vEy2Ideal - (" + vEy2Ideal + ")");


        double vEEx1 = pvx11 / p0 - px1 * vEx1 / p0;
        double vEEx2 = pvx21 / p0 - px1 * vEx2 / p0;
        double vEEy1 = pvy11 / p0 - py1 * vEy1 / p0;
        double vEEy2 = pvy21 / p0 - py1 * vEy2 / p0;

//        double vEEx1Ideal = 0.5 * bmu2 * cost1 * sint1;
//        double vEEx2Ideal = 0.5 * bmu2 * cost2 * sint2;
//        double vEEy1Ideal = -0.5 * bmu2 * cost1 * sint1;
//        double vEEy2Ideal = -0.5 * bmu2 * cost2 * sint2;


//        System.out.println("vEEx1Ideal  -(" +vEEx1Ideal  + ")");
//        System.out.println("vEEx2Ideal  -(" +vEEx2Ideal  + ")");
//        System.out.println("vEEy1Ideal  -(" + vEEy1Ideal + ")");
//        System.out.println("vEEy2Ideal  -(" +vEEy2Ideal  + ")");


        double dvEx1dt1 = dpvx10dt1 / px0 + pvx10 * dRpx0dt1;
        double dvEy1dt1 = dpvy10dt1 / py0 + pvy10 * dRpy0dt1;

//        double dvEx1dt1Ideal = -bmu * cost1;
//        double dvEy1dt1Ideal = -bmu * sint1;


        double dvEEx1dt1 = dpvx11dt1 / px0 + pvx11 * dRpx0dt1 - dvEx1dt1 * px1 / px0 - vEx1 * dpx1Rpx0dt1;
        double dvEEy1dt1 = dpvy11dt1 / py0 + pvy11 * dRpy0dt1 - dvEy1dt1 * py1 / py0 - vEy1 * dpy1Rpy0dt1;
        double d2vEx1dt1dt1 = d2pvx10dt1dt1 / px0 + dpvx10dt1 * dRpx0dt1 + dpvx10dt1 * dRpx0dt1 + pvx10 * d2Rpx0dt1dt1;
        double d2vEy1dt1dt1 = d2pvy10dt1dt1 / py0 + dpvy10dt1 * dRpy0dt1 + dpvy10dt1 * dRpy0dt1 + pvy10 * d2Rpy0dt1dt1;

//        double dvEEx1dt1Ideal = 0.5 * bmu2 * cos2t1;
//        double dvEEy1dt1Ideal = -dvEEx1dt1Ideal;


//        double d2vEx1dt1dt1Ideal = bmu * sint1;
//        double d2vEy1dt1dt1Ideal = -bmu * cost1;

//        System.out.println("dvEEx1dt1Ideal  -(" +dvEEx1dt1Ideal  + ")");
//        System.out.println("dvEEy1dt1Ideal  -(" +dvEEy1dt1Ideal  + ")");
//        System.out.println("d2vEx1dt1dt1Ideal  -(" +d2vEx1dt1dt1Ideal  + ")");
//        System.out.println("d2vEy1dt1dt1Ideal  -(" +d2vEy1dt1dt1Ideal  + ")");
//        System.out.println("vEx1Total= " + agentAtom1.vEx().getX(0));
//        System.out.println("agentAtom1.vEy() " + agentAtom1.vEy().getX(0));

        agentAtom1.vEx().PE(vEx1);
        agentAtom1.vEy().PE(vEy1);
        agentAtom1.vEEx().PE(vEEx1);
        agentAtom1.vEEy().PE(vEEy1);
        agentAtom1.dvEx().PE(dvEx1dt1);
        agentAtom1.dvEy().PE(dvEy1dt1);
        agentAtom1.dvEEx().PE(dvEEx1dt1);
        agentAtom1.dvEEy().PE(dvEEy1dt1);
        agentAtom1.d2vEx().PE(d2vEx1dt1dt1);
        agentAtom1.d2vEy().PE(d2vEy1dt1dt1);

//        agentAtom1.vEx().PE(-vEx1Ideal);
//        agentAtom1.vEy().PE(-vEy1Ideal);
//        agentAtom1.vEEx().PE(-vEEx1Ideal);
//        agentAtom1.vEEy().PE(-vEEy1Ideal);
//        agentAtom1.dvEx().PE(-dvEx1dt1Ideal);
//        agentAtom1.dvEy().PE(-dvEy1dt1Ideal);
//        agentAtom1.dvEEx().PE(-dvEEx1dt1Ideal);
//        agentAtom1.dvEEy().PE(-dvEEy1dt1Ideal);
//        agentAtom1.d2vEx().PE(-d2vEx1dt1dt1Ideal);
//        agentAtom1.d2vEy().PE(-d2vEy1dt1dt1Ideal);


        double dvEx2dt2 = dpvx20dt2 / px0 + pvx20 * dRpx0dt2;
        double dvEy2dt2 = dpvy20dt2 / py0 + pvy20 * dRpy0dt2;


//        double dvEx2dt2Ideal = -bmu * cost2;
//        double dvEy2dt2Ideal = -bmu * sint2;

//        System.out.println("dvEx2dt2Ideal  -(" + dvEx2dt2Ideal + ")");
//        System.out.println("dvEy2dt2Ideal  -(" + dvEy2dt2Ideal + ")");

        double dvEEx2dt2 = dpvx21dt2 / px0 + pvx21 * dRpx0dt2 - (dvEx2dt2 * px1 / px0 + vEx2 * dpx1Rpx0dt2);
        double dvEEy2dt2 = dpvy21dt2 / py0 + pvy21 * dRpy0dt2 - (dvEy2dt2 * py1 / py0 + vEy2 * dpy1Rpy0dt2);
        double d2vEx2dt2dt2 = d2pvx20dt2dt2 / px0 + dpvx20dt2 * dRpx0dt2 + dpvx20dt2 * dRpx0dt2 + pvx20 * d2Rpx0dt2dt2;
        double d2vEy2dt2dt2 = d2pvy20dt2dt2 / py0 + dpvy20dt2 * dRpy0dt2 + dpvy20dt2 * dRpy0dt2 + pvy20 * d2Rpy0dt2dt2;

//        double dvEEx2dt2Ideal = 0.5 * bmu2 * cos2t2;
//        double dvEEy2dt2Ideal = -dvEEx2dt2Ideal;

//        System.out.println("dvEEx2dt2Ideal  -(" + dvEEx2dt2Ideal + ")");
//        System.out.println("dvEEy2dt2Ideal  -(" + dvEEy2dt2Ideal + ")");


//        double d2vEx2dt2dt2Ideal = bmu * sint2;
//        double d2vEy2dt2dt2Ideal = -bmu * cost2;

//        System.out.println("d2vEx2dt2dt2Ideal  -(" + d2vEx2dt2dt2Ideal + ")");
//        System.out.println("d2vEy2dt2dt2Ideal  -(" + d2vEy2dt2dt2Ideal + ")");


        agentAtom2.vEx().PE(vEx2);
        agentAtom2.vEy().PE(vEy2);
        agentAtom2.vEEx().PE(vEEx2);
        agentAtom2.vEEy().PE(vEEy2);
        agentAtom2.dvEx().PE(dvEx2dt2);
        agentAtom2.dvEy().PE(dvEy2dt2);
        agentAtom2.dvEEx().PE(dvEEx2dt2);
        agentAtom2.dvEEy().PE(dvEEy2dt2);
        agentAtom2.d2vEx().PE(d2vEx2dt2dt2);
        agentAtom2.d2vEy().PE(d2vEy2dt2dt2);


//        agentAtom2.vEx().PE(-vEx2Ideal);
//        agentAtom2.vEy().PE(-vEy2Ideal);
//        agentAtom2.vEEx().PE(-vEEx2Ideal);
//        agentAtom2.vEEy().PE(-vEEy2Ideal);
//        agentAtom2.dvEx().PE(-dvEx2dt2Ideal);
//        agentAtom2.dvEy().PE(-dvEy2dt2Ideal);
//        agentAtom2.dvEEx().PE(-dvEEx2dt2Ideal);
//        agentAtom2.dvEEy().PE(-dvEEy2dt2Ideal);
//        agentAtom2.d2vEx().PE(-d2vEx2dt2dt2Ideal);
//        agentAtom2.d2vEy().PE(-d2vEy2dt2dt2Ideal);

//        System.out.println("vEx = " + vEx1);
//        System.out.println("atom1 =" + atom1.getLeafIndex());

//        System.out.println("vEy  = " + vEy1);
//        System.out.println("vEyIdeal  = " + vEy1Ideal);
//        System.out.println("agentAtom1.vEy() " + agentAtom1.vEy().getX(0));
//        System.out.println("vEEx = " + vEEx1);
//        System.out.println("vEEy = " + vEEy1 );
//        System.out.println("dvEx = " + dvEx1dt1);
//        System.out.println("dvEy = " + dvEy1dt1);
//        System.out.println("dvEEx = " + dvEEx1dt1 );
//        System.out.println("dvEEy = " + dvEEy1dt1 );
//        System.out.println("d2vEx = " + d2vEx1dt1dt1);
//        System.out.println("d2vEy = " + d2vEy1dt1dt1);


//        double f1 = bt * agentAtom1.torque.getX(0);
//        double f2 = bt * agentAtom2.torque.getX(0);

        if (debug) {
            if (atom2.getLeafIndex() == 1) {
//                System.out.println("t2 = " + t1 );
                System.out.println("t1 = " + t1 + ";");
                System.out.println("t2 = " + t2 + ";");
//                System.out.println("f1 = " + f1);
//                System.out.println("f2 = " + f2);
//            System.out.println("vEx21 = " + vEx2) ;
//            System.out.println("vEx2Ideal-( " + vEx2Ideal +" )");
//                System.out.println("d2vEx2_12= " + d2vEx2dt2dt2);
//                System.out.println("d2vEx2_Ideal = " + d2vEx2dt2dt2Ideal);
            }
            if (atom2.getLeafIndex() == 2) {
                System.out.println("t2 = " + t1 + ";");
                System.out.println("t3 = " + t2 + ";");
//                System.out.println("f2 = " + f1);
//                System.out.println("f3 = " + f2);
//            System.out.println("vEx23 = " + vEx1) ;
//            System.out.println("vEx2Ideal-( " + vEx1Ideal +" )");
//            System.out.println(" " + (agentAtom2.vEx().getX(0) - vEx1Ideal));
//                System.out.println("d2vEx2_23= " + d2vEx1dt1dt1);
//                System.out.println("d2vEx2_Ideal = " + d2vEx1dt1dt1Ideal);
//                System.out.println("Total= " + agentAtom1.d2vEx().getX(0));
//                System.out.println("TotalPlusIdeal= " + (agentAtom1.d2vEx().getX(0) + d2vEx1dt1dt1Ideal));
            }
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

//        phi.E(0);

        if (leafAgentIterator != null) {
            leafAgentIterator.reset();
            while (leafAgentIterator.hasNext()) {
                MoleculeAgent agent = (MoleculeAgent) leafAgentIterator.next();
                agent.vEx().E(0);
                agent.vEy().E(0);
                agent.vEEx().E(0);
                agent.vEEy().E(0);
                agent.dvEx().E(0);
                agent.dvEy().E(0);
                agent.dvEEx().E(0);
                agent.dvEEy().E(0);
                agent.d2vEx().E(0);
                agent.d2vEy().E(0);
            }
        }
    }

    public double getSumJEEMJEJE() {
        return JEEMJEJE;
    }

    public double getSumUEE() {
        return UEE;
    }

}
