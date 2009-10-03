package etomica.modules.rosmosis;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
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
        gradient = new IVectorMutable[]{work};
    }
    
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    public double getEpsilon() {
        return epsilon;
    }

    public double energy(IAtomList atoms) {
        IAtom atom = atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME((IVectorMutable)agentManager.getAgent(atom));
        return 0.5 * epsilon * work.squared();
    }

    public IVector[] gradient(IAtomList atoms) {
        IAtom atom = atoms.getAtom(0);
        work.E(atom.getPosition());
        work.ME((IVectorMutable)agentManager.getAgent(atom));
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
        return IVectorMutable.class;
    }

    public Object makeAgent(IAtom a) {
        if (a.getType().getSpecies() == species) {
            IVectorMutable vec = space.makeVector();
            vec.E(a.getPosition());
            return vec;
        }
        return null;
    }

    public void releaseAgent(Object agent, IAtom atom) {
        /* do nothing */
    }

    private static final long serialVersionUID = 1L;
    protected final AtomLeafAgentManager agentManager;
    protected final ISpecies species;
    protected final IVectorMutable work;
    protected final IVectorMutable[] gradient;
    protected double epsilon;
}
