package etomica.graphics;
import etomica.*;

/**
 * Button that attaches to a controller to toggle its pause/resume state.
 * Performs toggling of button label with state.
 */
public class DeviceControllerButton extends DeviceButton {
    
    private Controller controller;
    
    public DeviceControllerButton(Simulation sim) {
        super(sim);
    }
    public DeviceControllerButton(Simulation sim, Controller c) {
        this(sim);
        controller = c;
        setAction(new ActionGraphic(new Toggle(c)));
        setLabel("Start");
    }
    
    public void setController(Controller c) {controller = c;}
    public Controller getController() {return controller;}
    
    private class Toggle extends etomica.action.ControllerToggle {
        Toggle(Controller c) {this.setController(c);}
        public void actionPerformed(Controller c) {
            super.actionPerformed(c);
            String text;
            if(c.isPaused()) text = "Continue";
            else text = "Pause";
            DeviceControllerButton.this.setLabel(text);
        }
    }//end Toggle
}//end DeviceControllerButton