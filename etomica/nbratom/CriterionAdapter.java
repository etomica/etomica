/*
 * Created on Mar 2, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.nbratom;

import etomica.Atom;
import etomica.Phase;
import etomica.nbr.NeighborCriterion;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public abstract class CriterionAdapter extends NeighborCriterion {

    public CriterionAdapter(NeighborCriterion criterion) {
        subCriterion = criterion;
    }
    
    /**
     * Indicates whether the neighbor list for the given atom should
     * be updated, according to this criterion.  
     * @return true if the atom's list should be updated.
     */
    public boolean needUpdate(Atom atom) {
        return subCriterion.needUpdate(atom);
    }
    
    /**
     * Specifies the phase where the criterion is being applied.  Sometimes
     * needed if the criterion depends on features of the phase, such as the
     * volume.
     */
    public void setPhase(Phase phase) {
        subCriterion.setPhase(phase);
    }
    
    /**
     * Indicates whether the atom has changed (e.g. moved) by an amount 
     * that might have caused its neighbor list to be invalid.  If this
     * method returns true, a neighbor list failure may have introduced
     * errors in the calculation.
     */
    public boolean unsafe() {
        return subCriterion.unsafe();
    }
    
    /**
     * Indicates to criterion that given atom's neighbor list has just been 
     * updated, and that properties (e.g., record of atom's position) 
     * used by needUpdate and unsafe() methods should be reset.  
     */
    public void reset(Atom atom) {
        subCriterion.reset(atom);
    }
    
    /**
     * Returns the nominal distance within which two atoms are considered
     * neighbors.
     */
    public double getNeighborRange() {
        return subCriterion.getNeighborRange(); 
    }

    protected NeighborCriterion subCriterion;
}
