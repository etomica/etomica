package etomica;

import etomica.units.Dimension;
import etomica.utility.NameMaker;

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
        return Dimension.LENGTH;
    }
    
    public abstract double energy(AtomSet atoms);
    
    public abstract void setPhase(Phase phase);
        
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