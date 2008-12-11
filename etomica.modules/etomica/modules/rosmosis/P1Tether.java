package etomica.modules.rosmosis;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;


/**
 * Potential for a harmonic tether.
 * 
 * U = 0.5 eps r^2
 * 
 * @author Andrew Schultz
 */
public class P1Tether extends Potential1 implements AgentSource, PotentialSoft {

    public P1Tether(IBox box, ISpecies species, ISpace _space) {
        super(_space);
        this.species = species;
        agentManager = new AtomLeafAgentManager(this, box);
        work = _space.makeVector();
        gradient = new IVector[]{work};
    }
    
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    public double getEpsilon() {
        return epsilon;
    }

    public double energy(IAtomList atoms) {
        IAtomPositioned atom = (IAtomPositioned)atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME((IVector)agentManager.getAgent((IAtomLeaf)atom));
        return 0.5 * epsilon * work.squared();
    }

    public IVector[] gradient(IAtomList atoms) {
        IAtomPositioned atom = (IAtomPositioned)atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME((IVector)agentManager.getAgent((IAtomLeaf)atom));
        work.TE(epsilon);
        return gradient;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    /* AgentSource interface */
    public Class getAgentClass() {
        return IVector.class;
    }

    public Object makeAgent(IAtomLeaf a) {
        if (a.getType().getSpecies() == species) {
            IVector vec = space.makeVector();
            vec.E(((IAtomPositioned)a).getPosition());
            return vec;
        }
        return null;
    }

    public void releaseAgent(Object agent, IAtomLeaf atom) {
        /* do nothing */
    }

    private static final long serialVersionUID = 1L;
    protected final AtomLeafAgentManager agentManager;
    protected final ISpecies species;
    protected final IVector work;
    protected final IVector[] gradient;
    protected double epsilon;
}
