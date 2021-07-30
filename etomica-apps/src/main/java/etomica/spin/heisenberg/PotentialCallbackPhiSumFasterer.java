package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialAtomicSecondDerivative;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Tensor;

public class PotentialCallbackPhiSumFasterer implements PotentialCallback {
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;

    @Override
    public boolean wantsHessian() {
        return true;
    }

    @Override
    public void pairComputeHessian(int i, int j, Tensor phi) {
    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof IPotentialAtomicSecondDerivative)) {
            return;
        }

        IPotentialAtomicSecondDerivative potentialSecondDerivative = (IPotentialAtomicSecondDerivative) potential;
        Tensor[] t = potentialSecondDerivative.secondDerivative(atoms);

        IAtomOriented atom1 = (IAtomOriented) atoms.get(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.get(1);

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
    }

    public void reset() {
        leafAgentManager.getAgents().values().forEach((agent) -> agent.phi().E(0));

    }


}
