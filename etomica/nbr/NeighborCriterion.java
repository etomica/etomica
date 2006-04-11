/*
 * History
 * Created on Sep 22, 2004 by kofke
 */
package etomica.nbr;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.phase.Phase;

/**
 * Atom filter used to specify whether two atoms are considered neighbors,
 * for the purpose of tabulating neighbor lists.  Neighbors are indicated 
 * if the AtomsetFilter.accept() method returns true.  Adds neighbor management
 * methods to set the phase where the criterion is being applied, and to
 * indicate if neighbor list for atom needs updating.
 */
public interface NeighborCriterion {

    public boolean accept(AtomSet pair);
    
    /**
     * Indicates whether the neighbor list for the given atom should
     * be updated, according to this criterion.  
     * @return true if the atom's list should be updated.
     */
	public boolean needUpdate(Atom atom);
	
    /**
     * Specifies the phase where the criterion is being applied.  Sometimes
     * needed if the criterion depends on features of the phase, such as the
     * boundary.
     */
	public void setPhase(Phase phase);
	
    /**
     * Indicates whether the atom has changed (e.g. moved) by an amount 
     * that might have caused its neighbor list to be invalid.  If this
     * method returns true, a neighbor list failure may have introduced
     * errors in the calculation.
     */
	public boolean unsafe();
	
    /**
     * Indicates to criterion that given atom's neighbor list has just been 
     * updated, and that properties (e.g., record of atom's position) 
     * used by needUpdate and unsafe() methods should be reset.  
     */
	public void reset(Atom atom);
    
    /**
     * Returns true if the criterion (or subcriterion) depends on the distance
     * between atoms.
     */
    public boolean isRangeDependent();
}
