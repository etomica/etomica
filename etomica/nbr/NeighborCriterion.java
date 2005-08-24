/*
 * History
 * Created on Sep 22, 2004 by kofke
 */
package etomica.nbr;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomPairFilter;
import etomica.phase.Phase;

/**
 * Atom filter used to specify whether two atoms are considered neighbors,
 * for the purpose of tabulating neighbor lists.  Neighbors are indicated 
 * if the AtomsetFilter.accept() method returns true.  Adds neighbor management
 * methods to set the phase where the criterion is being applied, and to
 * indicate if neighbor list for atom needs updating.
 */
public abstract class NeighborCriterion implements AtomPairFilter, java.io.Serializable {

    public abstract boolean accept(AtomPair pair);
    
    /**
     * Indicates whether the neighbor list for the given atom should
     * be updated, according to this criterion.  
     * @return true if the atom's list should be updated.
     */
	public abstract boolean needUpdate(Atom atom);
	
    /**
     * Specifies the phase where the criterion is being applied.  Sometimes
     * needed if the criterion depends on features of the phase, such as the
     * volume.
     */
	public abstract void setPhase(Phase phase);
	
    /**
     * Indicates whether the atom has changed (e.g. moved) by an amount 
     * that might have caused its neighbor list to be invalid.  If this
     * method returns true, a neighbor list failure may have introduced
     * errors in the calculation.
     */
	public abstract boolean unsafe();
	
    /**
     * Indicates to criterion that given atom's neighbor list has just been 
     * updated, and that properties (e.g., record of atom's position) 
     * used by needUpdate and unsafe() methods should be reset.  
     */
	public abstract void reset(Atom atom);
    
    /**
     * Returns the nominal distance within which two atoms are considered
     * neighbors.
     */
    public abstract double getNeighborRange();
    
	//TODO consider ways to ensure this is removed from nbrmanager if no longer used
    
    public static final NeighborCriterion ALL = new NeighborCriterion() {
        /**
         * Always returns false, indicating that neighbor list never needs updating.
         * This is appropriate if atoms are never added to or removed from phase,
         * because all atoms are always on neighbor list.
         */
        public boolean needUpdate(Atom atom) {return false;}
        public double getNeighborRange() {return Double.POSITIVE_INFINITY;}
        /**
         * Performs no action.
         */
        public void setPhase(Phase phase) {}
        /**
         * Always returns false, indicating that neighbor list never needs updating.
         * This is appropriate if atoms are never added to or removed from phase,
         * because all atoms are always on neighbor list.
         */
        public boolean unsafe() {return false;}
        /**
         * Performs no action.
         */
        public void reset(Atom atom) {}
        /**
         * Always returns true, indicating that all atoms pairs are neighbors.
         */
        public boolean accept(AtomPair pair) {return true;}
    };
    
    public static final NeighborCriterion NONE = new NeighborCriterion() {
        /**
         * Always returns false, indicating that neighbor list never needs updating.
         * This is appropriate if atoms are never added to or removed from phase,
         * because all atoms are always on neighbor list.
         */
        public boolean needUpdate(Atom atom) {return false;}
        public double getNeighborRange() {return 0.0;}
        /**
         * Performs no action.
         */
        public void setPhase(Phase phase) {}
        /**
         * Always returns false, indicating that neighbor list never needs updating.
         */
        public boolean unsafe() {return false;}
        /**
         * Performs no action.
         */
        public void reset(Atom atom) {}
        /**
         * Always returns false, indicating that no atoms pairs are neighbors.
         */
        public boolean accept(AtomPair pair) {return false;}
    };

}
