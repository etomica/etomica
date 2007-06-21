package etomica.potential;

import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Length;

/**
 * Superclass for all Potential classes, which define how the atoms in the
 * system interact with each other.
 *
 * @author David Kofke
 */
 
public abstract class Potential implements java.io.Serializable, IPotential {
    
	private final int nBody;
	protected final Space space;

    /**
     * General constructor for a potential instance
     * @param nBody number of atoms to which this potential applies at a time;
     * for example with a pair potential nBody = 2; for a single-body potential,
     * nBody = 1.
     */
    public Potential(int nBody, Space space) {
        this.nBody = nBody;
        this.space = space;
    }

    public Space getSpace() {
        return space;
    }

    public abstract double getRange();
    
    public Dimension getRangeDimension() {
        return Length.DIMENSION;
    }
    
    /**
     * Returns the interaction energy between the given atoms.  There might be
     * 0, 1, 2 or more atoms in the AtomSet.
     */
    public abstract double energy(AtomSet atoms);
    
    /**
     * Informs the potential of the phase on which it acts. Typically this
     * requires at least that it update the nearestImageTransformer of its
     * coordinatePair (if it uses one), e.g.:
     * cPair.setNearestImageTransformer(phase.boundary());
     */
    public abstract void setPhase(Phase phase);
    
    /**
     * The number of atoms on which the potential depends.
     */
    public final int nBody() {return nBody;}
    

}