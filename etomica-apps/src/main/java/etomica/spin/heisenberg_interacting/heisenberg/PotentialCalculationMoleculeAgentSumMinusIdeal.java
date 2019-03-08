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
    protected final double mu, J, bt, bJ, bmu;
    //    I0bJ, I1bJ, I2bJ;
//    protected double[] Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
//    protected double[] d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
//    protected double[] Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
//    protected double[] d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;
//    protected double[] InbJArray, Inm1bJArray, Inm2bJArray, Inp1bJArray;

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
        this.nMax = nMax;
        this.leafAgentManager = leafAgentManager;
        leafAgentIterator = leafAgentManager.makeIterator();

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


        count += 1;
//        System.out.println("(*" + count + "th term*)");
        boolean debug = false;
        if (debug) {
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

//        double cos2t1 = 2 * cost1 * cost1 - 1;
//        double sin2t1 = 2 * sint1 * cost1;

        MoleculeAgent agentAtom1 = leafAgentManager.getAgent(atom1);
        MoleculeAgent agentAtom2 = leafAgentManager.getAgent(atom2);

         double pvx10, pvx20, dpvx10dt1, d2pvx20dt2dt2, pvx11, pvx21, dpvx11dt1,   dpvx21dt2, dpvx20dt2, d2pvx10dt1dt1, dpvx10dt2;
         double pvy10, pvy20, dpvy10dt1,  d2pvy20dt2dt2, pvy11, pvy21, dpvy11dt1,  dpvy21dt2, dpvy20dt2, d2pvy10dt1dt1, dpvy10dt2;


        double sint2 = ej.getX(1);
        double cost2 = ej.getX(0);
        double cost1mt2 = cost1 * cost2 + sint1 * sint2;
        double sint1mt2 = sint1 * cost2 - cost1 * sint2;
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
        dpvx11dt1 = 0;
        dpvy11dt1 = 0;
        pvx20 = 0;
        pvy20 = 0;
        pvx21 = 0;
        pvy21 = 0;
        dpvx20dt2 = 0;
        dpvy20dt2 = 0;
        dpvx21dt2 = 0;
        dpvy21dt2 = 0;
        d2pvx10dt1dt1 = 0;
        d2pvy10dt1dt1 = 0;
        d2pvx20dt2dt2 = 0;
        d2pvy20dt2dt2 = 0;


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

            pvx11 += agentAtom1.dAxs1[n] * sinnt2 + agentAtom1.dAxc1[n] * cosnt2;
            pvy11 += agentAtom1.dAys1[n] * sinnt2 + agentAtom1.dAyc1[n] * cosnt2;


            dpvx10dt1 += agentAtom1.d2Axs0[n] * sinnt2 + agentAtom1.d2Axc0[n] * cosnt2;
            dpvy10dt1 += agentAtom1.d2Ays0[n] * sinnt2 + agentAtom1.d2Ayc0[n] * cosnt2;

            d2pvx10dt1dt1 += agentAtom1.d3Axs0[n] * sinnt2 + agentAtom1.d3Axc0[n] * cosnt2;
            d2pvy10dt1dt1 += agentAtom1.d3Ays0[n] * sinnt2 + agentAtom1.d3Ayc0[n] * cosnt2;


//            dpvx10dt2 += n * agentAtom1.dAxs0[n] * cosnt2 - n * agentAtom1.dAxc0[n] * sinnt2;
//            dpvy10dt2 += n * agentAtom1.dAys0[n] * cosnt2 - n * agentAtom1.dAyc0[n] * sinnt2;

//            d2pvx10dt1dt2 += n * agentAtom1.d2Axs0[n] * cosnt2 - n * agentAtom1.d2Axc0[n] * sinnt2;
//            d2pvy10dt1dt2 += n * agentAtom1.d2Ays0[n] * cosnt2 - n * agentAtom1.d2Ayc0[n] * sinnt2;


            dpvx11dt1 += agentAtom1.d2Axs1[n] * sinnt2 + agentAtom1.d2Axc1[n] * cosnt2;
            dpvy11dt1 += agentAtom1.d2Ays1[n] * sinnt2 + agentAtom1.d2Ayc1[n] * cosnt2;

//            dpvx11dt2 += n * agentAtom1.dAxs1[n] * cosnt2 - n * agentAtom1.dAxc1[n] * sinnt2;
//            dpvy11dt2 += n * agentAtom1.dAys1[n] * cosnt2 - n * agentAtom1.dAyc1[n] * sinnt2;

            if (n != 0) {
                pvx20 += n * agentAtom1.Axs0[n] * cosnt2 - n * agentAtom1.Axc0[n] * sinnt2;
                pvy20 += n * agentAtom1.Ays0[n] * cosnt2 - n * agentAtom1.Ayc0[n] * sinnt2;

                pvx21 += n * agentAtom1.Axs1[n] * cosnt2 - n * agentAtom1.Axc1[n] * sinnt2;
                pvy21 += n * agentAtom1.Ays1[n] * cosnt2 - n * agentAtom1.Ayc1[n] * sinnt2;

//                dpvx20dt1 += n * agentAtom1.dAxs0[n] * cosnt2 - n * agentAtom1.dAxc0[n] * sinnt2;
//                dpvy20dt1 += n * agentAtom1.dAys0[n] * cosnt2 - n * agentAtom1.dAyc0[n] * sinnt2;

//                d2pvx20dt1dt2 += -n2 * agentAtom1.dAxs0[n] * sinnt2 - n2 * agentAtom1.dAxc0[n] * cosnt2;
//                d2pvy20dt1dt2 += -n2 * agentAtom1.dAys0[n] * sinnt2 - n2 * agentAtom1.dAyc0[n] * cosnt2;

                dpvx20dt2 += -n2 * agentAtom1.Axs0[n] * sinnt2 - n2 * agentAtom1.Axc0[n] * cosnt2;
                dpvy20dt2 += -n2 * agentAtom1.Ays0[n] * sinnt2 - n2 * agentAtom1.Ayc0[n] * cosnt2;

                d2pvx20dt2dt2 += -n2 * n * agentAtom1.Axs0[n] * cosnt2 + n2 * n * agentAtom1.Axc0[n] * sinnt2;
                d2pvy20dt2dt2 += -n2 * n * agentAtom1.Ays0[n] * cosnt2 + n2 * n * agentAtom1.Ayc0[n] * sinnt2;

//                dpvx21dt1 += n * agentAtom1.dAxs1[n] * cosnt2 - n * agentAtom1.dAxc1[n] * sinnt2;
//                dpvy21dt1 += n * agentAtom1.dAys1[n] * cosnt2 - n * agentAtom1.dAyc1[n] * sinnt2;

                dpvx21dt2 += -n2 * agentAtom1.Axs1[n] * sinnt2 - n2 * agentAtom1.Axc1[n] * cosnt2;
                dpvy21dt2 += -n2 * agentAtom1.Ays1[n] * sinnt2 - n2 * agentAtom1.Ayc1[n] * cosnt2;
            }
        }

        double vEx1 = pvx10 / p0;
        double vEx2 = pvx20 / p0;
        double vEy1 = pvy10 / p0;
        double vEy2 = pvy20 / p0;



        double vEEx1 = pvx11 / p0 - px1 * vEx1 / p0;
        double vEEx2 = pvx21 / p0 - px1 * vEx2 / p0;
        double vEEy1 = pvy11 / p0 - py1 * vEy1 / p0;
        double vEEy2 = pvy21 / p0 - py1 * vEy2 / p0;



        double dvEx1dt1 = dpvx10dt1 / px0 + pvx10 * dRpx0dt1;
        double dvEy1dt1 = dpvy10dt1 / py0 + pvy10 * dRpy0dt1;



        double dvEEx1dt1 = dpvx11dt1 / px0 + pvx11 * dRpx0dt1 - dvEx1dt1 * px1 / px0 - vEx1 * dpx1Rpx0dt1;
        double dvEEy1dt1 = dpvy11dt1 / py0 + pvy11 * dRpy0dt1 - dvEy1dt1 * py1 / py0 - vEy1 * dpy1Rpy0dt1;
        double d2vEx1dt1dt1 = d2pvx10dt1dt1 / px0 + dpvx10dt1 * dRpx0dt1 + dpvx10dt1 * dRpx0dt1 + pvx10 * d2Rpx0dt1dt1;
        double d2vEy1dt1dt1 = d2pvy10dt1dt1 / py0 + dpvy10dt1 * dRpy0dt1 + dpvy10dt1 * dRpy0dt1 + pvy10 * d2Rpy0dt1dt1;


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



        double dvEx2dt2 = dpvx20dt2 / px0 + pvx20 * dRpx0dt2;
        double dvEy2dt2 = dpvy20dt2 / py0 + pvy20 * dRpy0dt2;

        double dvEEx2dt2 = dpvx21dt2 / px0 + pvx21 * dRpx0dt2 - (dvEx2dt2 * px1 / px0 + vEx2 * dpx1Rpx0dt2);
        double dvEEy2dt2 = dpvy21dt2 / py0 + pvy21 * dRpy0dt2 - (dvEy2dt2 * py1 / py0 + vEy2 * dpy1Rpy0dt2);
        double d2vEx2dt2dt2 = d2pvx20dt2dt2 / px0 + dpvx20dt2 * dRpx0dt2 + dpvx20dt2 * dRpx0dt2 + pvx20 * d2Rpx0dt2dt2;
        double d2vEy2dt2dt2 = d2pvy20dt2dt2 / py0 + dpvy20dt2 * dRpy0dt2 + dpvy20dt2 * dRpy0dt2 + pvy20 * d2Rpy0dt2dt2;



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
