package etomica.potential;

import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.util.NameMaker;

/**
 * Superclass for all Potential classes, which define how the atoms in the
 * system interact with each other.
 *
 * @author David Kofke
 */
 
public abstract class Potential {
    
	private final int nBody;
	protected final Space space;
    private String name;

    /**
     * General constructor for a potential instance
     * @param nBody number of atoms to which this potential applies at a time;
     * for example with a pair potential nBody = 2; for a single-body potential,
     * nBody = 1.
     */
    public Potential(int nBody, Space space) {
        this.nBody = nBody;
        this.space = space;
        setName(NameMaker.makeName(this.getClass()));
    }

    public Space getSpace() {
        return space;
    }

    public abstract double getRange();
    
    public Dimension getRangeDimension() {
        return Length.DIMENSION;
    }
    
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
    

    /**
     * @return Returns the name.
     */
    public String getName() {
        return name;
    }
    /**
     * @param name The name to set.
     */
    public void setName(String name) {
        this.name = name;
    }
}//end of Potential