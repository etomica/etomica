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
    
    public DeviceControllerButton(SimulationElement parent) {
        super(parent);
    }
    public DeviceControllerButton(SimulationElement parent, Controller c) {
        this(parent);
        setController(c);
    }
    
    //final because called by constructor
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