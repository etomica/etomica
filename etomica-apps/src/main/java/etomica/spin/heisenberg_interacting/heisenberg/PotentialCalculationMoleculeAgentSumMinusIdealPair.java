package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialAtomicSecondDerivative;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * Anx expressions for computing v_E and v_EE in the mapping.
 *
 * @author Weisong Lin
 */

public class PotentialCalculationMoleculeAgentSumMinusIdealPair implements PotentialCalculation {
    protected AtomLeafAgentManager.AgentIterator leafAgentIterator;
    protected Vector ei, ej;
    protected double JEEMJEJE, UEE, U_Map, U_Con;
    protected final double mu, J, bt, bJ, bmu;


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


        double cost1 = ei.getX(0);
        double sint1 = ei.getX(1);
        double sint2 = ej.getX(1);
        double cost2 = ej.getX(0);
        double cost1mt2 = cost1 * cost2 + sint1 * sint2;
        double sint1mt2 = sint1 * cost2 - cost1 * sint2;


        double p0 = Math.exp(bJ * cost1mt2);
        double px0 = p0, py0 = p0;
        double dRpx0dt1 = bJ * sint1mt2 / px0;//Rpx0 mean reverse of px0  D[1/px0, theta1]
        double dRpx0dt2 = -dRpx0dt1; // D[1/px0, theta2]
        double d2Rpx0dt1dt1 = bJ * (cost1mt2 + bJ * sint1mt2 * sint1mt2) / px0;
        double d2Rpx0dt1dt2 = -d2Rpx0dt1dt1;

        double dRpy0dt1 = bJ * sint1mt2 / py0;//Rpy0 mean reverse of py0  D[1/py0, theta1]
        double dRpy0dt2 = -dRpy0dt1; // D[1/py0, theta2]
        double d2Rpy0dt1dt1 = bJ * (cost1mt2 + bJ * sint1mt2 * sint1mt2) / py0;
        double d2Rpy0dt1dt2 = -d2Rpy0dt1dt1;

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
            double n2 = n * n;


            pvx10 += agentAtom1.dAxs0[n] * agentAtom2.sinntheta[n] + agentAtom1.dAxc0[n] * agentAtom2.cosntheta[n];
            pvy10 += agentAtom1.dAys0[n] * agentAtom2.sinntheta[n] + agentAtom1.dAyc0[n] * agentAtom2.cosntheta[n];


            dpvx10dt1 += agentAtom1.d2Axs0[n] * agentAtom2.sinntheta[n] + agentAtom1.d2Axc0[n] * agentAtom2.cosntheta[n];
            dpvy10dt1 += agentAtom1.d2Ays0[n] * agentAtom2.sinntheta[n] + agentAtom1.d2Ayc0[n] * agentAtom2.cosntheta[n];


            dpvx10dt2 += n * agentAtom1.dAxs0[n] * agentAtom2.cosntheta[n] - n * agentAtom1.dAxc0[n] * agentAtom2.sinntheta[n];
            dpvy10dt2 += n * agentAtom1.dAys0[n] * agentAtom2.cosntheta[n] - n * agentAtom1.dAyc0[n] * agentAtom2.sinntheta[n];

            d2pvx10dt1dt2 += n * agentAtom1.d2Axs0[n] * agentAtom2.cosntheta[n] - n * agentAtom1.d2Axc0[n] * agentAtom2.sinntheta[n];
            d2pvy10dt1dt2 += n * agentAtom1.d2Ays0[n] * agentAtom2.cosntheta[n] - n * agentAtom1.d2Ayc0[n] * agentAtom2.sinntheta[n];


            if (n != 0) {
                pvx20 += n * agentAtom1.Axs0[n] * agentAtom2.cosntheta[n] - n * agentAtom1.Axc0[n] * agentAtom2.sinntheta[n];
                pvy20 += n * agentAtom1.Ays0[n] * agentAtom2.cosntheta[n] - n * agentAtom1.Ayc0[n] * agentAtom2.sinntheta[n];


                dpvx20dt1 += n * agentAtom1.dAxs0[n] * agentAtom2.cosntheta[n] - n * agentAtom1.dAxc0[n] * agentAtom2.sinntheta[n];
                dpvy20dt1 += n * agentAtom1.dAys0[n] * agentAtom2.cosntheta[n] - n * agentAtom1.dAyc0[n] * agentAtom2.sinntheta[n];

                d2pvx20dt1dt2 += -n2 * agentAtom1.dAxs0[n] * agentAtom2.sinntheta[n] - n2 * agentAtom1.dAxc0[n] * agentAtom2.cosntheta[n];
                d2pvy20dt1dt2 += -n2 * agentAtom1.dAys0[n] * agentAtom2.sinntheta[n] - n2 * agentAtom1.dAyc0[n] * agentAtom2.cosntheta[n];

                dpvx20dt2 += -n2 * agentAtom1.Axs0[n] * agentAtom2.sinntheta[n] - n2 * agentAtom1.Axc0[n] * agentAtom2.cosntheta[n];
                dpvy20dt2 += -n2 * agentAtom1.Ays0[n] * agentAtom2.sinntheta[n] - n2 * agentAtom1.Ayc0[n] * agentAtom2.cosntheta[n];
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

        U_Map += 0.5 * J * (f1-f2) * sint1mt2;

        U_Con += -J * cost1mt2;

    }


    public void zeroSum() {
        JEEMJEJE = 0;
        UEE = 0;
        U_Map = 0;
        U_Con = 0;
    }

    public double getSumJEEMJEJE() {
        return JEEMJEJE;
    }

    public double getSumUEE() {
        return UEE;
    }

    public double getSumU_Map() {
        return U_Map;
    }

    public double getSumU_Con() {
        return U_Con;
    }

}
