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

public class PotentialCalculationMoleculeAgentSumMinusIdealPair implements PotentialCalculation {
    //public class PotentialCalculationHeisenberg {
    protected AtomLeafAgentManager.AgentIterator leafAgentIterator;
    protected Vector ei, ej;
    protected double AEEJ0, JEMUExIdeal, JEMUEyIdeal, JEMUEIdealSquare, JEEMJEJE, UEE, JEMUExSquare, JEMUEySquare, JEMUEx, JEMUEy, dipolex, dipoley, JEEMJEJExtrying, UEEnow, JEMUE, dipoleconv;
    protected final double mu, J, bt, bJ, bmu, I0bJ, I1bJ, I2bJ;
    protected double[] Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    protected double[] d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
    protected double[] Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
    protected double[] d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;
    protected double[] InbJArray, Inm1bJArray, Inm2bJArray, Inp1bJArray;


    protected int nMax;
    protected int count = 1;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;

    public PotentialCalculationMoleculeAgentSumMinusIdealPair(Space space, double dipoleMagnitude, double interactionS, double beta, int nMax, AtomLeafAgentManager<MoleculeAgent> leafAgentManager) {
        ei = space.makeVector();
        ej = space.makeVector();


        J = interactionS;
        mu = dipoleMagnitude;
        bt = beta;
        bJ = bt * J;
        bmu = bt * mu;
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



//        System.out.println("(*" + count + "th term*)");
        boolean debug = false;
        if (debug) {
            count += 1;
            System.out.println("nMax= " + nMax + ";");
//            System.out.println("mu= " + mu + ";");
//            System.out.println("J= " + J + ";");
//            System.out.println("bJ= " + bJ + ";");
//            System.out.println("bt = bJ/J;");
//            System.out.println("bmu = bt*mu;");
            System.out.println("t1 = " + t1 + ";");
            System.out.println("t2 = " + t2 + ";");
        }

//        if (count > 1000 && debug) System.exit(2);


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
        double cos2t1 = 2 * cost1 - 1;
        double sin2t1 = 2 * sint1 * cost1;


        for (int n = 0; n <= nMax; n++) {
            if (n == 0) {
                Axc0[n] = bmu * (I0bJ + I1bJ) * (-1 + cost1);
                dAxc0[n] = -bmu * (I0bJ + I1bJ) * sint1;
                Axs0[n] = 0;
                dAxs0[n] = 0;
                Axc1[n] = -0.25 * bmu * bmu * sint1 * sint1 * (I0bJ + 2 * I1bJ + I2bJ);
                dAxc1[n] = -0.25 * bmu * bmu * (I0bJ + 2 * I1bJ + I2bJ) * Math.sin(2 * t1);
                Axs1[n] = 0;
                dAxs1[n] = 0;

                d2Axc0[n] = -bmu * cost1 * (I0bJ + I1bJ);
                d2Axs0[n] = 0;
                d3Axc0[n] = bmu * sint1 * (I0bJ + I1bJ);
                d3Axs0[n] = 0;
                d2Axc1[n] = -0.5 * bmu * bmu * cos2t1 * (I0bJ + 2 * I1bJ + I2bJ);
                d2Axs1[n] = 0;


                //y direction
                Ayc0[n] = bmu * (I0bJ + I1bJ) * sint1;
                dAyc0[n] = bmu * (I0bJ + I1bJ) * cost1;
                Ays0[n] = 0;
                dAys0[n] = 0;
                Ayc1[n] = 0.25 * bmu * bmu * sint1 * sint1 * (I0bJ + 2 * I1bJ + I2bJ);
                dAyc1[n] = 0.25 * bmu * bmu * Math.sin(2 * t1) * (I0bJ + 2 * I1bJ + I2bJ);
                Ays1[n] = 0;
                dAys1[n] = 0;

                d2Ayc0[n] = -bmu * sint1 * (I0bJ + I1bJ);
                d2Ays0[n] = 0;
                d3Ayc0[n] = -bmu * cost1 * (I0bJ + I1bJ);
                d3Ays0[n] = 0;
                d2Ayc1[n] = 0.5 * bmu * bmu * (I0bJ + 2 * I1bJ + I2bJ) * cos2t1;
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
                ) * bmu * bmu / 4.0 / n2;

                dAxc1[n] = 0.25 * bmu * bmu * (
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
                ) * bmu * bmu / 4.0 / n2;

                dAxs1[n] = 0.25 * bmu * bmu * (((n - 2) * Inm2bJ * cosnm2t1) / (2 - 2 * n + n2)
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

                d2Axc1[n] = 0.25 * bmu * bmu * (-(n - 2) * (n - 2) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (-bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ) + 4 * n * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        - 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                d2Axs1[n] = 0.25 * bmu * bmu * (-((n - 2) * (n - 2) * Inm2bJ * sinnm2t1) / (2 - 2 * n + n2)
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
                ) * bmu * bmu / 4 / n2;

                dAyc1[n] = 0.25 * bmu * bmu * (((-2 + n) * Inm2bJ * sinnm2t1) / (2 - 2 * n + n2)
                        + (bJ * (bJ * n * (-4 - 2 * n + n3) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ) + (16 + 4 * n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ)))
                        + 2 * n * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );

                Ays1[n] = (-bJ * n2 * (2 + n * (2 + n)) * I0bJ * sinnm2t1 * ((-1 + bJ + n) * Inm1bJ + bJ * InbJ)
                        - 2 * bJ * (4 + n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                        - n2 * (2 - 2 * n + n2) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n * (1 + n)) * InbJ)
                ) * bmu * bmu / (2 * bJ * bJ * n2 * (4 + n4) * I0bJ);

                dAys1[n] = 0.25 * bmu * bmu * (-(-2 + n) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
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

                d2Ayc1[n] = 0.25 * bmu * bmu * ((-2 + n) * (-2 + n) * Inm2bJ * cosnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (bJ * (-2 + n) * (-2 + n) * n * (2 + 2 * n + n2) * I0bJ * cosnm2t1 * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * cosnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * cosnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
                ) / (bJ * bJ * n * (4 + n4) * I0bJ)
                );


                d2Ays1[n] = 0.25 * bmu * bmu * ((n - 2) * (n - 2) * Inm2bJ * sinnm2t1 / (2 - 2 * n + n2)
                        + (bJ * (bJ * (n - 2) * (n - 2) * n * (2 + 2 * n + n2) * I0bJ * sinnm2t1 * (2 * Inm1bJ + InbJ)
                        + 4 * n * (4 + n4) * sinnt1 * (bJ * I1bJ * InbJ + I0bJ * (-bJ * Inm1bJ + n * InbJ))
                )
                        + 2 * n * (2 + n) * (4 - 2 * n + n3) * I0bJ * sinnp2t1 * (bJ * (-1 + bJ - n) * Inm1bJ + (bJ * bJ - 2 * bJ * n + 2 * n + 2 * n2) * InbJ)
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
        double dpy1Rpy0dt1 = bmu * cost1;//D[py1/py0, theta1]
        double dpy1Rpy0dt2 = bmu * cost2;
        double d2Rpy0dt1dt1 = bJ * (Math.cos(t1 - t2) + bJ * Math.sin(t1 - t2) * Math.sin(t1 - t2)) / py0;
        double d2Rpy0dt1dt2 = -d2Rpy0dt1dt1;
        double d2Rpy0dt2dt2 = d2Rpy0dt1dt1;

        double pvx10, pvx20, dpvx10dt1, dpvx20dt1, d2pvx10dt1dt2, d2pvx20dt1dt2, dpvx20dt2, dpvx10dt2;
        double pvy10, pvy20, dpvy10dt1, dpvy20dt1, d2pvy10dt1dt2, d2pvy20dt1dt2, dpvy20dt2, dpvy10dt2;


        pvx10 = 0;
        pvy10 = 0;
        dpvx10dt1 = 0;
        dpvy10dt1 = 0;
        dpvx10dt2 = 0;
        dpvy10dt2 = 0;
        d2pvx10dt1dt2 = 0;
        d2pvy10dt1dt2 = 0;
        pvx20 = 0;
        pvy20 = 0;
        dpvx20dt1 = 0;
        dpvy20dt1 = 0;
        d2pvx20dt1dt2 = 0;
        dpvx20dt2 = 0;
        dpvy20dt2 = 0;
        d2pvy20dt1dt2 = 0;

        for (int n = 0; n <= nMax; n++) {
            double sinnt2 = Math.sin(n * t2);
            double cosnt2 = Math.cos(n * t2);
            double n2 = n * n;


            pvx10 += dAxs0[n] * sinnt2 + dAxc0[n] * cosnt2;
            pvy10 += dAys0[n] * sinnt2 + dAyc0[n] * cosnt2;


            dpvx10dt1 += d2Axs0[n] * sinnt2 + d2Axc0[n] * cosnt2;
            dpvy10dt1 += d2Ays0[n] * sinnt2 + d2Ayc0[n] * cosnt2;


            dpvx10dt2 += n * dAxs0[n] * cosnt2 - n * dAxc0[n] * sinnt2;
            dpvy10dt2 += n * dAys0[n] * cosnt2 - n * dAyc0[n] * sinnt2;

            d2pvx10dt1dt2 += n * d2Axs0[n] * cosnt2 - n * d2Axc0[n] * sinnt2;
            d2pvy10dt1dt2 += n * d2Ays0[n] * cosnt2 - n * d2Ayc0[n] * sinnt2;


            if (n != 0) {
                pvx20 += n * Axs0[n] * cosnt2 - n * Axc0[n] * sinnt2;
                pvy20 += n * Ays0[n] * cosnt2 - n * Ayc0[n] * sinnt2;


                dpvx20dt1 += n * dAxs0[n] * cosnt2 - n * dAxc0[n] * sinnt2;
                dpvy20dt1 += n * dAys0[n] * cosnt2 - n * dAyc0[n] * sinnt2;

                d2pvx20dt1dt2 += -n2 * dAxs0[n] * sinnt2 - n2 * dAxc0[n] * cosnt2;
                d2pvy20dt1dt2 += -n2 * dAys0[n] * sinnt2 - n2 * dAyc0[n] * cosnt2;

                dpvx20dt2 += -n2 * Axs0[n] * sinnt2 - n2 * Axc0[n] * cosnt2;
                dpvy20dt2 += -n2 * Ays0[n] * sinnt2 - n2 * Ayc0[n] * cosnt2;
            }
        }

        MoleculeAgent agentAtom1 = leafAgentManager.getAgent(atom1);
        MoleculeAgent agentAtom2 = leafAgentManager.getAgent(atom2);

        double dvEx1dt2 = dpvx10dt2 / px0 + pvx10 * dRpx0dt2;
        double dvEy1dt2 = dpvy10dt2 / py0 + pvy10 * dRpy0dt2;


        double d2vEx1dt1dt2 = d2pvx10dt1dt2 / px0 + dpvx10dt1 * dRpx0dt2 + dpvx10dt2 * dRpx0dt1 + pvx10 * d2Rpx0dt1dt2;
        double d2vEy1dt1dt2 = d2pvy10dt1dt2 / py0 + dpvy10dt1 * dRpy0dt2 + dpvy10dt2 * dRpy0dt1 + pvy10 * d2Rpy0dt1dt2;


        double dvEx2dt1 = dpvx20dt1 / px0 + pvx20 * dRpx0dt1;
        double dvEy2dt1 = dpvy20dt1 / py0 + pvy20 * dRpy0dt1;

        double d2vEx2dt1dt2 = d2pvx20dt1dt2 / px0 + dpvx20dt2 * dRpx0dt1 + dpvx20dt1 * dRpx0dt2 + pvx20 * d2Rpx0dt1dt2;
        double d2vEy2dt1dt2 = d2pvy20dt1dt2 / py0 + dpvy20dt2 * dRpy0dt1 + dpvy20dt1 * dRpy0dt2 + pvy20 * d2Rpy0dt1dt2;


        double f1 = bt * agentAtom1.torque.getX(0);
        double f2 = bt * agentAtom2.torque.getX(0);
        double p12 = -bJ * ei.dot(ej);
        double p21 = p12;


        double vEx1 = agentAtom1.vEx().getX(0);
        double vEy1 = agentAtom1.vEy().getX(0);

        double vEx2 = agentAtom2.vEx().getX(0);
        double vEy2 = agentAtom2.vEy().getX(0);


        JEEMJEJE += vEx1 * d2vEx2dt1dt2 + vEx2 * d2vEx1dt1dt2;
        JEEMJEJE += vEy1 * d2vEy2dt1dt2 + vEy2 * d2vEy1dt1dt2;

        UEE -= f1 * vEx2 * dvEx1dt2 + f2 * vEx1 * dvEx2dt1;
        UEE -= f1 * vEy2 * dvEy1dt2 + f2 * vEy1 * dvEy2dt1;

        UEE += 2 * vEx1 * p12 * vEx2;
        UEE += 2 * vEy1 * p12 * vEy2;

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
        JEMUEIdealSquare = 0;

    }

    public double getSumJEEMJEJE() {
        return JEEMJEJE;
    }

    public double getSumUEE() {
        return UEE;
    }

}
