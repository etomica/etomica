package etomica.action;
import etomica.*;

/**
 * Switches controller between paused and unpaused conditions, or starts
 * it if not yet active.
 */
public class ControllerToggle extends ControllerAction {
        
 /*   public ControllerToggle() {
        super(parentSimulation());
    }*/
        
        
    public void actionPerformed(Controller c) {
        doAction(c);
    }
        
    public static void doAction(Controller c) {
        if(!c.isActive()) c.start();
        else if(!c.isPaused()) c.pause();
        else c.unPause();
	}
}//end of ControllerToggle
    