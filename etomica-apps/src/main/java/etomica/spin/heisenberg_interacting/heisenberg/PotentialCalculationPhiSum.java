package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.integrator.Integrator;
import etomica.molecule.DipoleSource;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.util.numerical.BesselFunction;

public class PotentialCalculationPhiSum implements PotentialCalculation {
    protected AtomLeafAgentManager<MeterMappedAveraging.MoleculeAgent> leafAgentManager;
    protected AtomLeafAgentManager.AgentIterator leafAgentIterator;


    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }

        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
        Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);

        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);

//        System.out.println("Mapping pairs:(" + atom1 + " , " + atom2 + ")");
//        System.out.println(  atom1 +"  = " + t[1].component(0,0) );
//        System.out.println(  atom2 +"  = " + t[2].component(0,0) );
//        if (atom1.getLeafIndex() == 1) {
//            System.out.println("beforesum= " +leafAgentManager.getAgent(atom1).phi().component(0, 0));
//        }
//        if (atom2.getLeafIndex() == 1) {
//            System.out.println("beforesum= " + leafAgentManager.getAgent(atom2).phi().component(0, 0));
//        }

        leafAgentManager.getAgent(atom1).phi().PE(t[1]);
        leafAgentManager.getAgent(atom2).phi().PE(t[2]);

//        if (atom1.getLeafIndex() == 1) {
//            System.out.println("aftersum= " + leafAgentManager.getAgent(atom1).phi().component(0, 0));
//        }
//        if (atom2.getLeafIndex() == 1) {
//            System.out.println("aftersum= " + leafAgentManager.getAgent(atom2).phi().component(0, 0));
//        }


    }

    public void setAgentManager(AtomLeafAgentManager agentManager) {
        leafAgentManager = agentManager;
        leafAgentIterator = leafAgentManager.makeIterator();
    }

    public void reset() {

        leafAgentIterator.reset();
        while (leafAgentIterator.hasNext()) {
            Object agent = leafAgentIterator.next();
            ((MeterMappedAveraging.MoleculeAgent) agent).phi().E(0);
        }


    }


}
