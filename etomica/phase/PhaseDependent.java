/*
 * History
 * Created on Jul 25, 2004 by kofke
 */
package etomica.phase;


/**
 * @author kofke
 *
 * Interface for an object that depends on a phase.
 */
public interface PhaseDependent {
	
    /**
     * Sets the Phase on which the object acts.
     */
    public void setPhase(Phase phase);
    
    /**
     * Accessor method for the phase on which the object performs
     * its activities.
     */
     public Phase getPhase();
        
}
