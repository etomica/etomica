package etomica.zeolite;
//This class includes a main method to demonstrate its use
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

import etomica.api.IController;

import etomica.EtomicaInfo;
import etomica.action.Action;
import etomica.graphics.Device;

/**
 * Button that causes an action to be performed.
 * 
 * @author David Kofke
 */
 
public class DeviceButtonSingle extends Device {
    
    /**
     * Constructs button that is not connected to any action.  Subsequent
     * call to setAction is needed to make the button do something.
     * @param controller
     */
    public DeviceButtonSingle(IController controller) {
        super(controller);
        button = new JButton();
    }
    
    /**
     * Constructs a button connected to the given action.  Controller
     * and action may be changed independently after construction.
     */
    public DeviceButtonSingle(IController controller, Action action) {
        this(controller);
        setAction(action);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Button that causes an elementary action to be performed");
        return info;
    }
    
    /**
     * Performs the button's action. Same effect as if the button were pressed
     * in the user interface.
     */
    public void press() {
        button.doClick();
    }
    
    /**
     * Returns the currently defined action associated with the button.
     */
    public etomica.action.Action getAction() {return targetAction;}

    /**
     * Defines the action to be performed when the button is pressed.
     */
    public void setAction(final Action newAction) {
        if(buttonAction != null) button.removeActionListener(buttonAction);
        targetAction = newAction;
        if(newAction == null) return;
        buttonAction = new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                try {
                    doAction(targetAction);
                    button.setEnabled(false);
                }
                catch (RuntimeException e) {
                    //do nothing
                }
            }
        };
        button.addActionListener(buttonAction);
    }
    
    /**
     * Returns the GUI element for display in the simulation.
     */
    public java.awt.Component graphic(Object obj) {
        return button;
    }
    
    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(String text) {button.setText(text);}
    /**
     * @return the current instance of the descriptive label.
     */
    public String getLabel() {return button.getText();}
    
    private ActionListener buttonAction;
    protected JButton button;
    protected Action targetAction;
    
    /**
     * Method to demonstrate and test the use of this class.  
     * Slider is used to control the temperature of a hard-sphere MD simulation
     */
    public static void main(String[] args) {
        final String APP_NAME = "Device Button Single";

        etomica.simulation.prototypes.HSMD2D sim = new etomica.simulation.prototypes.HSMD2D();
        etomica.graphics.SimulationGraphic graphic = new etomica.graphics.SimulationGraphic(sim, APP_NAME, sim.getSpace());
        
        //here's the part unique to this class
        etomica.action.SimulationRestart action = new etomica.action.SimulationRestart(sim, sim.getSpace());
        DeviceButtonSingle button = new DeviceButtonSingle(sim.getController(),action);
        //end of unique part
        graphic.add(button);
        graphic.makeAndDisplayFrame(APP_NAME);
    }
    
}