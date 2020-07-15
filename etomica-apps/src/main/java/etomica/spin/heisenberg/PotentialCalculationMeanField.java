package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * Potential Calculation for mean field
 */

public class PotentialCalculationMeanField implements PotentialCalculation, AtomLeafAgentManager.AgentSource<Vector> {
    protected final double J2;
    protected final Space space;

    protected AtomLeafAgentManager<Vector> leafAgentManager;

    public PotentialCalculationMeanField(Space space, double J, Box box) {
        this.J2 = 0.5 * J;
        this.space = space;
        this.leafAgentManager = new AtomLeafAgentManager<Vector>(this, box);
    }


    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        IAtomOriented atom1 = (IAtomOriented) atoms.get(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.get(1);
        Vector ei = atom1.getOrientation().getDirection();
        Vector ej = atom2.getOrientation().getDirection();
        Vector hi = leafAgentManager.getAgent(atom1);
        Vector hj = leafAgentManager.getAgent(atom2);

        hi.PEa1Tv1(J2, ej);
        hj.PEa1Tv1(J2, ei);
    }


    public void reset() {
        leafAgentManager.getAgents().values().forEach((agent) -> agent.E(0));
    }

    public AtomLeafAgentManager<Vector> getAgentManager() {
        return leafAgentManager;
    }

    @Override
    public Vector makeAgent(IAtom a, Box agentBox) {
        return space.makeVector();
    }

    @Override
    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {

    }
}
