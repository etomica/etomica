package etomica.graphics;
import etomica.*;

/**
 * Button that attaches to a controller to toggle its pause/resume state.
 * Performs toggling of button label with state.
 */
 
 /* History of changes
  * 7/03/02 (DAK/SKK) Added reset method to change label to "Start"
  */
  
public class DeviceControllerButton extends DeviceButton {
    
    private Controller controller;
    
    public DeviceControllerButton(Simulation sim) {
        super(sim);
    }
    public DeviceControllerButton(Simulation sim, Controller c) {
        this(sim);
        setController(c);
    }
    public DeviceControllerButton(Space space) {
        super(space);
    }
    /**
     * Constructor if button is to be used as part of another device.
     * Does not register with simulation.
     */
    public DeviceControllerButton(Space space, Controller c) {
        super(space);
        setController(c);
    }
    
    //final because called by contructor
    public final void setController(Controller c) {
        controller = c;
        setAction(new ActionGraphic(new Toggle(c)));
        setLabel("  Start  ");
    }
    public Controller getController() {return controller;}
    
    /**
     * Sets label of button to display "Start".
     */
    public void reset() {setLabel("  Start  ");}
    
    private class Toggle extends etomica.action.ControllerToggle {
        Toggle(Controller c) {this.setController(c);}
        public void actionPerformed(Controller c) {
            super.actionPerformed(c);
            String text;
            if(c.isPaused()) text = "Continue";
            else text = "  Pause ";
            DeviceControllerButton.this.setLabel(text);
        }
    }//end Toggle
}//end DeviceControllerButton