/*
 * History
 * Created on Oct 27, 2004 by kofke
 */
package etomica.action;

/**
 * Toggles the state of enabling the long-range-corrections to the 
 * potential truncation in the phase.
 */
public class PhaseToggleLrc extends PhaseActionAdapter {
    
	public PhaseToggleLrc() {
		super("Toggle LRC");
	}
	
    public void actionPerformed() {
    	if(phase == null) return;
        phase.setLrcEnabled(!phase.isLrcEnabled());
    }
}