package etomica.spin.heisenberg_interacting.heisenberg;

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
    protected AtomLeafAgentManager.AgentIterator leafAgentIterator;
    protected final double J;
    protected final Space space;

    protected AtomLeafAgentManager<Vector> leafAgentManager;

    public PotentialCalculationMeanField(Space space, double J, Box box) {
        this.J = J;
        this.space = space;
        this.leafAgentManager = new AtomLeafAgentManager<>(this, box, Vector.class);

        leafAgentIterator = leafAgentManager.makeIterator();

    }


    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(1);
        Vector ei = atom1.getOrientation().getDirection();
        Vector ej = atom2.getOrientation().getDirection();
        Vector hi = leafAgentManager.getAgent(atom1);
        Vector hj = leafAgentManager.getAgent(atom2);

        hi.PE(ej);
        hj.PE(ei);
    }


    public void reset() {
        AtomLeafAgentManager.AgentIterator<Vector> it = leafAgentManager.makeIterator();
        for (Vector h = it.next(); it.hasNext(); h = it.next()) {
            h.E(0);
        }
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
