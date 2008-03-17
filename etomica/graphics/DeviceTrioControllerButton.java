package etomica.graphics;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.action.activity.Controller;
import etomica.api.IController;
import etomica.api.ISimulation;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.HSMD2D;
import etomica.space.Space;

/**
 * Device comprising three buttons: (1) attaches to a controller to toggle its pause/resume state; 
 * (2)causes a reset (restart) action of simulation to be performed; (3) causes reset of meter averages
 * only.
 *
 * @author Sang Kyu Kwak
 */
 
public class DeviceTrioControllerButton extends Device {
    
    private JPanel jp;
    private DeviceControllerButton startButton;
    private ISimulation simulation;
    private DeviceButton reinitButton;
    private DeviceButton resetButton;
	private double width;
	private boolean firstResized = true;
	private String shape;
	private SimulationRestart simRestart;

    /**
     * Contructs device with buttons that affect the given simulation.
     */
    public DeviceTrioControllerButton(ISimulation simulation, Space space) {
        this();
        setSimulation(simulation, space);
    }
    
    /**
     * No-argument contructor gives device that performs no action until
     * setSimulation is used to specify the target simulation.
     */
    public DeviceTrioControllerButton() {
        super();
        jp = new JPanel(new java.awt.GridLayout(1,0, 20, 20)); //default shape of panel
        jp.setBorder(new TitledBorder(null, "Control", TitledBorder.CENTER, TitledBorder.TOP));
        jp.setOpaque(false);

        startButton = new DeviceControllerButton(null);
        reinitButton = new DeviceButton(null);
        resetButton = new DeviceButton(null);
        reinitButton.setLabel("Reinitialize");
        resetButton.setLabel("Reset averages");
        
        jp.add(startButton.graphic()); 
        jp.add(reinitButton.graphic());  
        jp.add(resetButton.graphic());
                       
        setShape("VERTICAL");
    }
        
    /**
     * Sets the controller that is toggled by this device.
     */
    protected void setSimulation(ISimulation sim, Space space) {
        simulation = sim;
        simRestart = new SimulationRestart(sim, space);
        final Controller c = (Controller)sim.getController();
        setController(c);
        startButton.setController(c);
        reinitButton.setPreAction(new Action() {
            public void actionPerformed() {
                if (c.isActive()) {
                    c.halt();
                    DeviceTrioControllerButton.this.reset();
                }
                c.reset();
            }
        });
        resetButton.setController(c);
        reinitButton.setAction(simRestart);
        resetButton.setAction(simRestart.getDataResetAction());
    }

    /**
     * Returns an ArrayList of DataPumps that are headers for the data streams.
     * This streams from this ArrayList are reset when the user presses the 
     * reinitialize and reset averages buttons.
     */
    public ArrayList getDataStreamPumps() {
        return simRestart.getDataResetAction().getDataStreamPumps();
    }
    
    /**
     * Returns the simulation affected by the buttons.
     */
    public ISimulation getSimulation() {
        return simulation;
    }

    /**
     * Returns the JPanel that contains the buttons.
     */
    public java.awt.Component graphic(Object obj) {
        return jp;
    }
    
    public DeviceControllerButton getControllerButton() {return startButton;}
    public DeviceButton getReinitButton() {return reinitButton;}
    public DeviceButton getResetAveragesButton() {return resetButton;}
    
    /**
     * Sets controller toggle button to read "Start"
     */
    public void reset() {startButton.reset();}
    
    /**
     * Sets display shape of the device.  The argument must be a string having
     * one of the following values:<ul>
     * <li>"HORIZONTAL" causes buttons to be laid out in a row.
     * <li>"VERTICAL" causes buttons to be laid out in a column (default)
     * <li>"AUTOMATIC" causes buttons to adopt layout based on size of available space
     * </ul>
     * Any other value is ignored.
     */
    public void setShape(String s){
        shape = s;
        if(s=="HORIZONTAL"){jp.setLayout(new java.awt.GridLayout(1,0,4,4));jp.updateUI();}
        if(s=="VERTICAL"){jp.setLayout(new java.awt.GridLayout(0,1,4,4));jp.updateUI();}
        if(s=="AUTOMATIC"){jp.getParent().addComponentListener(new ComponentEventControllerButton());}
    }
    public String getShape() {return shape;}

    public SimulationRestart getSimRestart() {
    	return(simRestart);
    }

    /**
     * Inner class that catches action of simulation panel 
     */        
    protected class ComponentEventControllerButton implements ComponentListener {
  
        public void componentHidden(ComponentEvent e){}
        public void componentMoved(ComponentEvent e){}
        public void componentShown(ComponentEvent e){}
        public void componentResized(ComponentEvent e){
            if(firstResized) { width = jp.getParent().getWidth();firstResized =false;}
            if((double)jp.getParent().getWidth()<width){
                jp.setLayout(new java.awt.GridLayout(3,1));
                jp.updateUI();
            } else { 
                jp.setLayout(new java.awt.GridLayout(1,3)); 
                jp.updateUI();
            }
        }
    }
 
    /**
     * main method to show how to work with this class 
     */        
     public static void main(String[] args) {
        final String APP_NAME = "Device Trio Controller Button";

        etomica.space.Space sp = etomica.space2d.Space2D.getInstance();
        final HSMD2D sim = new HSMD2D(); 

        DeviceTrioControllerButton button = new DeviceTrioControllerButton(sim, sp);
            button.setShape("HORIZONTAL"); //three choices "HORIZONTAL", "AUTOMATIC"          
//        DeviceTrioControllerButton button = new DeviceTrioControllerButton(Simulation.instance, Simulation.instance.controller(0)); 
//          button.setShape("VERTICAL"); //three choices "HORIZONTAL", "AUTOMATIC"
        
        final SimulationGraphic graphic = new SimulationGraphic(sim, APP_NAME, sp);

        // Simulation Graphic will display it's own Trio button group by
        // default.  Just remove them and put ours on for this test.
        graphic.getPanel().controlPanel.removeAll();
        graphic.add(button);

        button.getReinitButton().setPostAction(graphic.getPaintAction(sim.box));

        graphic.makeAndDisplayFrame(APP_NAME);
    }
    
       
}//end DeviceTrioControllerButton
