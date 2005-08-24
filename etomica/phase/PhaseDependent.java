/*
 * History
 * Created on Jul 25, 2004 by kofke
 */
package etomica.phase;

import etomica.util.IntegerRange;

/**
 * @author kofke
 *
 * Interface for an object that depends on one or more phases.
 */
public interface PhaseDependent {
	
	public Phase[] getPhase();
	
	public void setPhase(Phase[] phase);
	
	public IntegerRange phaseCountRange();
}
