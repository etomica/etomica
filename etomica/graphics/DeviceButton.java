//This class includes a main method to demonstrate its use
package etomica.graphics;
import javax.swing.JButton;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;

/**
 * Button that causes an action to be performed.
 * 
 * @author David Kofke
 */
 
 /* History of changes
  * 07/03/02 (DAK/SKK) Added constructor that takes space as an argument, not registering button with a simulation.
  * 08/30/02 (DAK) Added default colors from DefaultGraphic
  * 10/07/02 (DAK) Removed default colors, now using UIManager for color defaults
  */
public class DeviceButton extends Device implements EtomicaElement {
    
    private ActionGraphic action;
    private JButton button;
    
    public DeviceButton() {
        super();
        button = new JButton();
    }
    
    /**
     * Constructs a button connected to the given action.
     */
    public DeviceButton(ActionGraphic action) {
        this();
        setAction(action);
    }
    /**
     * Connects button to given action.
     */
    public DeviceButton(etomica.Action action) {
        this(new ActionGraphic(action));
        setLabel(action.getLabel());
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Button that causes an elementary action to be performed");
        return info;
    }
    
    public etomica.Action getAction() {return action;}
    public void setAction(ActionGraphic newAction) {
        if(action != null) button.removeActionListener(action);
        action = newAction;
        button.addActionListener(action);
        button.setText(action.getLabel());
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
    
    /**
     * @return a handle to the JSlider instance used by this slider device
     */
    public JButton getButton() {return button;}
    
    /**
     * Method to demonstrate and test the use of this class.  
     * Slider is used to control the temperature of a hard-sphere MD simulation
     */
/*    public static void main(String[] args) {
        
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        
        //here's the part unique to this class
        etomica.action.SimulationRestart action = new etomica.action.SimulationRestart(sim);
        DeviceButton button = new DeviceButton(action);
        //end of unique part
 
        sim.elementCoordinator.go();
        Simulation.makeAndDisplayFrame(sim);
    }
    */
}