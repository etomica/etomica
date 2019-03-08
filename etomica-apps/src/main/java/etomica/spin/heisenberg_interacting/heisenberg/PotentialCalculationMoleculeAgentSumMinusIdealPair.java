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
//        double t1 = Math.atan2(ei.getX(1), ei.getX(0));
//        double t2 = Math.atan2(ej.getX(1), ej.getX(0));

        MoleculeAgent agentAtom1 = leafAgentManager.getAgent(atom1);
        MoleculeAgent agentAtom2 = leafAgentManager.getAgent(atom2);


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
        double cost1mt2 = cost1 * cost2 + sint1 * sint2;
        double sint1mt2 = sint1 * cost2 - cost1 * sint2;
//        double cos2t1 = 2 * cost1 - 1;
//        double sin2t1 = 2 * sint1 * cost1;



        double p0 = Math.exp(bJ * cost1mt2);
        double px0 = p0, py0 = p0;
//        double pM2 = 1 / p0 / p0;
//        double lnp1 = -bJ * sint1mt2;
//        double lnp11 = -bJ * cost1mt2;
//        double lnp2 = -lnp1;
//        double px1 = bmu * p0 * (cost1 + cost2);
//        double py1 = bmu * p0 * (sint1 + sint2);
        double dRpx0dt1 = bJ * sint1mt2 / px0;//Rpx0 mean reverse of px0  D[1/px0, theta1]
        double dRpx0dt2 = -dRpx0dt1; // D[1/px0, theta2]
//        double dpx1Rpx0dt1 = -bmu * sint1;//D[px1/px0, theta1]
//        double dpx1Rpx0dt2 = -bmu * sint2;
        double d2Rpx0dt1dt1 = bJ * (cost1mt2 + bJ * sint1mt2 * sint1mt2) / px0;
        double d2Rpx0dt1dt2 = -d2Rpx0dt1dt1;
//        double d2Rpx0dt2dt2 = d2Rpx0dt1dt1;

        double dRpy0dt1 = bJ * sint1mt2 / py0;//Rpy0 mean reverse of py0  D[1/py0, theta1]
        double dRpy0dt2 = -dRpy0dt1; // D[1/py0, theta2]
//        double dpy1Rpy0dt1 = bmu * cost1;//D[py1/py0, theta1]
//        double dpy1Rpy0dt2 = bmu * cost2;
        double d2Rpy0dt1dt1 = bJ * (cost1mt2 + bJ * sint1mt2 * sint1mt2) / py0;
        double d2Rpy0dt1dt2 = -d2Rpy0dt1dt1;
//        double d2Rpy0dt2dt2 = d2Rpy0dt1dt1;

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

        double[] sinxt2 = new double[nMax+1];
        double[] cosxt2 = new double[nMax+1];
        if (nMax>=0) {
            cosxt2[0] = 1;
            if (nMax >= 1) {
                sinxt2[1] = sint2;
                cosxt2[1] = cost2;
                if (nMax >= 2) {
                    sinxt2[2] = 2 * sint2 * cost2;
                    cosxt2[2] = 2 * cost2 * cost2 - 1;
                    if (nMax >= 3) throw new RuntimeException("fix me");
                }
            }
        }

        for (int n = 0; n <= nMax; n++) {
            double sinnt2 = sinxt2[n];
            double cosnt2 = cosxt2[n];
            double n2 = n * n;


            pvx10 += agentAtom1.dAxs0[n] * sinnt2 + agentAtom1.dAxc0[n] * cosnt2;
            pvy10 += agentAtom1.dAys0[n] * sinnt2 + agentAtom1.dAyc0[n] * cosnt2;


            dpvx10dt1 += agentAtom1.d2Axs0[n] * sinnt2 + agentAtom1.d2Axc0[n] * cosnt2;
            dpvy10dt1 += agentAtom1.d2Ays0[n] * sinnt2 + agentAtom1.d2Ayc0[n] * cosnt2;


            dpvx10dt2 += n * agentAtom1.dAxs0[n] * cosnt2 - n * agentAtom1.dAxc0[n] * sinnt2;
            dpvy10dt2 += n * agentAtom1.dAys0[n] * cosnt2 - n * agentAtom1.dAyc0[n] * sinnt2;

            d2pvx10dt1dt2 += n * agentAtom1.d2Axs0[n] * cosnt2 - n * agentAtom1.d2Axc0[n] * sinnt2;
            d2pvy10dt1dt2 += n * agentAtom1.d2Ays0[n] * cosnt2 - n * agentAtom1.d2Ayc0[n] * sinnt2;


            if (n != 0) {
                pvx20 += n * agentAtom1.Axs0[n] * cosnt2 - n * agentAtom1.Axc0[n] * sinnt2;
                pvy20 += n * agentAtom1.Ays0[n] * cosnt2 - n * agentAtom1.Ayc0[n] * sinnt2;


                dpvx20dt1 += n * agentAtom1.dAxs0[n] * cosnt2 - n * agentAtom1.dAxc0[n] * sinnt2;
                dpvy20dt1 += n * agentAtom1.dAys0[n] * cosnt2 - n * agentAtom1.dAyc0[n] * sinnt2;

                d2pvx20dt1dt2 += -n2 * agentAtom1.dAxs0[n] * sinnt2 - n2 * agentAtom1.dAxc0[n] * cosnt2;
                d2pvy20dt1dt2 += -n2 * agentAtom1.dAys0[n] * sinnt2 - n2 * agentAtom1.dAyc0[n] * cosnt2;

                dpvx20dt2 += -n2 * agentAtom1.Axs0[n] * sinnt2 - n2 * agentAtom1.Axc0[n] * cosnt2;
                dpvy20dt2 += -n2 * agentAtom1.Ays0[n] * sinnt2 - n2 * agentAtom1.Ayc0[n] * cosnt2;
            }
        }



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
//        double p21 = p12;


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
        JEMUE = 0;
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
