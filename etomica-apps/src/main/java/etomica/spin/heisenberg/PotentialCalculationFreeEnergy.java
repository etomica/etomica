package etomica.spin.heisenberg;

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

public class PotentialCalculationFreeEnergy implements PotentialCalculation {
    protected AtomLeafAgentManager.AgentIterator leafAgentIterator;
    protected Vector ei, ej;
    protected double U_Map;
    protected final double mu, J, bt, bJ, bmu;


    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;

    public PotentialCalculationFreeEnergy(Space space, double dipoleMagnitude, double interactionS, double beta, AtomLeafAgentManager<MoleculeAgent> leafAgentManager) {
        ei = space.makeVector();
        ej = space.makeVector();
        J = interactionS;
        mu = dipoleMagnitude;
        bt = beta;
        bJ = bt * J;
        bmu = bt * mu;
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
        MoleculeAgent agentAtom1 = leafAgentManager.getAgent(atom1);
        MoleculeAgent agentAtom2 = leafAgentManager.getAgent(atom2);


        double cost1 = ei.getX(0);
        double sint1 = ei.getX(1);
        double sint2 = ej.getX(1);
        double cost2 = ej.getX(0);
        double sint1mt2 = sint1 * cost2 - cost1 * sint2;

        double f1 = bt * agentAtom1.torque.getX(0);
        double f2 = bt * agentAtom2.torque.getX(0);
        U_Map += 0.5 * J * (f1-f2) * sint1mt2;



    }


    public void zeroSum() {
        U_Map = 0;
    }


    public double getSumU_Map() {
        return U_Map;
    }


}
