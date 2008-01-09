package etomica.modules.rosmosis;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.box.Box;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.IVector;
import etomica.space.Tensor;
import etomica.species.Species;

/**
 * Potential for a harmonic tether.
 * 
 * U = 0.5 eps r^2
 * 
 * @author Andrew Schultz
 */
public class P1Tether extends Potential1 implements AgentSource, PotentialSoft {

    public P1Tether(Box box, Species species) {
        super(box.getSpace());
        this.species = species;
        agentManager = new AtomAgentManager(this, box);
        work = box.getSpace().makeVector();
        gradient = new IVector[]{work};
    }
    
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    public double getEpsilon() {
        return epsilon;
    }

    public double energy(AtomSet atoms) {
        IAtomPositioned atom = (IAtomPositioned)atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME((IVector)agentManager.getAgent(atom));
        return 0.5 * epsilon * work.squared();
    }

    public IVector[] gradient(AtomSet atoms) {
        IAtomPositioned atom = (IAtomPositioned)atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME((IVector)agentManager.getAgent(atom));
        work.TE(epsilon);
        return gradient;
    }

    public IVector[] gradient(AtomSet atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(AtomSet atoms) {
        return 0;
    }

    /* AgentSource interface */
    public Class getAgentClass() {
        return IVector.class;
    }

    public Object makeAgent(IAtom a) {
        if (a.getType().getSpecies() == species && a instanceof IAtomPositioned) {
            IVector vec = space.makeVector();
            vec.E(((IAtomPositioned)a).getPosition());
            return vec;
        }
        return null;
    }

    public void releaseAgent(Object agent, IAtom atom) {
        /* do nothing */
    }

    private static final long serialVersionUID = 1L;
    protected final AtomAgentManager agentManager;
    protected final Species species;
    protected final IVector work;
    protected final IVector[] gradient;
    protected double epsilon;
}
