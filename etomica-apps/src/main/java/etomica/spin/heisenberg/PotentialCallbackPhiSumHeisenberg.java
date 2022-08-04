package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.potential.IPotential2;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Space;
import etomica.space.Vector;

public class PotentialCallbackPhiSumHeisenberg implements PotentialCallback {
    protected Vector ei, ej;
    protected double secondDerivativeSum = 0;
    protected AtomLeafAgentManager<MoleculeAgent> leafAgentManager;


    public PotentialCallbackPhiSumHeisenberg(Space space, AtomLeafAgentManager<MoleculeAgent> leafAgentManager) {
        ei = space.makeVector();
        ej = space.makeVector();
        this.leafAgentManager = leafAgentManager;
    }

    @Override
    public boolean wantsHessian() {
        return true;
    }

    @Override
    public void pairComputeGeneral(IPotential2 pij, IAtom atom1, IAtom atom2, Vector drij, Vector fij, Vector tij, Vector tji) {
        IPotential2.Hessian h = pij.d2u(drij, atom1, atom2);
        leafAgentManager.getAgent(atom1).phi.PE(h.o1o1);
        leafAgentManager.getAgent(atom2).phi.PE(h.o2o2);

        ei.E(((IAtomOriented)atom1).getOrientation().getDirection());
        ej.E(((IAtomOriented)atom2).getOrientation().getDirection());

        secondDerivativeSum += 2 * h.o1o2.component(0, 0) * (ei.dot(ej));
        secondDerivativeSum += h.o1o1.component(0, 0);
        secondDerivativeSum += h.o2o2.component(0, 0);

//		double Cos = ei.dot(ej);
//		secondDerivativeSum += -2*Cos*Cos+2*Cos;
//		double diff = 2*t[0].component(0,0)*(ei.dot(ej))+t[1].component(0,0)+t[2].component(0,0)-1.5*( -2*Cos*Cos+2*Cos);
//		System.out.println("diff = " + diff);//check these two approach is the same or not
    }

    public void zeroSum() {
        secondDerivativeSum = 0.0;
        leafAgentManager.getAgents().values().forEach((agent) -> agent.phi().E(0));
    }

    /**
     * Returns the current value of the energy sum.
     */
    public double getSum() {
        return secondDerivativeSum;
    }


}
