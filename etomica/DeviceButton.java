//This class includes a main method to demonstrate its use
package etomica;
import javax.swing.JButton;

/**
 * Button that causes an action to be performed.
 * 
 * @author David Kofke
 */
public class DeviceButton extends Device implements EtomicaElement {
    
    public String getVersion() {return "DeviceButton:01.05.25/"+Device.VERSION;}

    private etomica.Action action;
    private JButton button;
    
    public DeviceButton() {
        this(Simulation.instance);
    }
    public DeviceButton(Simulation sim) {
        super(sim);
        button = new JButton();
    }
    
    /**
     * Constructs a slider connected to the given property of the given object
     */
    public DeviceButton(etomica.Action action) {
        this();
        setAction(action);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Button that causes an elementary action to be performed");
        return info;
    }
    
    public etomica.Action getAction() {return action;}
    public void setAction(etomica.Action newAction) {
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
     * Slider is used to control the temperature of a hard-disk MD simulation
     */
    public static void main(String[] args) {
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);
                
        
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        
        //here's the part unique to this class
        etomica.action.SimulationRestart action = new etomica.action.SimulationRestart(sim);
        DeviceButton button = new DeviceButton(action);
        //end of unique part
 
        Simulation.instance.elementCoordinator.go();
        f.getContentPane().add(sim.panel());         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }
}