package etomica.action;

/**
 * Toggles the state of enabling the long-range-corrections to the 
 * potential truncation in the phase.
 */
public class PhaseToggleLrc extends PhaseActionAdapter {
    
    public void actionPerformed() {
    	if(phase == null) return;
        phase.setLrcEnabled(!phase.isLrcEnabled());
    }

    private static final long serialVersionUID = 1L;
}
