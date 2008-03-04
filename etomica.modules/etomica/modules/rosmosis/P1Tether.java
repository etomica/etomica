package etomica.modules.rosmosis;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVector;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;


/**
 * Potential for a harmonic tether.
 * 
 * U = 0.5 eps r^2
 * 
 * @author Andrew Schultz
 */
public class P1Tether extends Potential1 implements AgentSource, PotentialSoft {

    public P1Tether(IBox box, ISpecies species, Space _space) {
        super(_space);
        this.species = species;
        agentManager = new AtomAgentManager(this, box);
        work = _space.makeVector();
        gradient = new IVector[]{work};
    }
    
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    public double getEpsilon() {
        return epsilon;
    }

    public double energy(IAtomSet atoms) {
        IAtomPositioned atom = (IAtomPositioned)atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME((IVector)agentManager.getAgent(atom));
        return 0.5 * epsilon * work.squared();
    }

    public IVector[] gradient(IAtomSet atoms) {
        IAtomPositioned atom = (IAtomPositioned)atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME((IVector)agentManager.getAgent(atom));
        work.TE(epsilon);
        return gradient;
    }

    public IVector[] gradient(IAtomSet atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomSet atoms) {
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
    protected final ISpecies species;
    protected final IVector work;
    protected final IVector[] gradient;
    protected double epsilon;
}
