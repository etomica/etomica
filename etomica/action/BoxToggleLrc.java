package etomica.action;

/**
 * Toggles the state of enabling the long-range-corrections to the 
 * potential truncation in the box.
 */
public class BoxToggleLrc extends BoxActionAdapter {
    
    public void actionPerformed() {
    	if(box == null) return;
        box.setLrcEnabled(!box.isLrcEnabled());
    }

    private static final long serialVersionUID = 1L;
}
