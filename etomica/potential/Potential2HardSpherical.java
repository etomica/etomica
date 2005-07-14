package etomica.potential;

import etomica.AtomPair;
import etomica.AtomSet;
import etomica.Phase;
import etomica.Space;
import etomica.space.CoordinatePair;

/**
 * Methods for a hard (impulsive), spherically-symmetric pair potential.
 * Subclasses must provide a concrete definition for the energy (method u(double)).
 */

public abstract class Potential2HardSpherical extends Potential2 implements PotentialHard, Potential2Spherical {
   
	public Potential2HardSpherical(Space space, CoordinatePair cPair) {
	    super(space);
        this.coordinatePair = cPair;
	}
	
	/**
    * The pair energy u(r^2) with no truncation applied.
    * @param the square of the distance between the particles.
    */
    public abstract double u(double r2);

    /**
     * Energy of the pair as given by the u(double) method, with application
     * of any PotentialTruncation that may be defined for the potential.
     */
    public double energy(AtomSet pair) {
    	    coordinatePair.reset((AtomPair)pair);
    	    return u(coordinatePair.r2());
    }
    
    public void setPhase(Phase phase) {
        coordinatePair.setNearestImageTransformer(phase.boundary());
    }

    protected final CoordinatePair coordinatePair;
    
}
