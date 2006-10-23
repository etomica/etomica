package etomica.graphics;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

import javax.swing.JPanel;

import etomica.action.ResetAccumulators;
import etomica.action.SimulationDataAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.Controller;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.HSMD2D;

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
    private Simulation simulation;
    private DeviceButton restartButton;
    private DeviceButton resetButton;
	private double width;
	private boolean firstResized = true;
	private String shape;

    /**
     * Contructs device with buttons that affect the given simulation.
     */
    public DeviceTrioControllerButton(Simulation simulation) {
        this();
        setSimulation(simulation);
    }
    
    /**
     * No-argument contructor gives device that performs no action until
     * setSimulation is used to specify the target simulation.
     */
    public DeviceTrioControllerButton() {
        super();
        jp = new JPanel(new java.awt.GridLayout(1,0)); //default shape of panel
        jp.setBorder(new javax.swing.border.TitledBorder("Control"));
        jp.setOpaque(false);
//        jp.setBackground(DefaultGraphic.BORDER_COLOR);
/*        jp.setBorder(new javax.swing.border.TitledBorder(
                         new javax.swing.border.EtchedBorder(
                             javax.swing.border.EtchedBorder.RAISED, java.awt.Color.red, java.awt.Color.blue) 
                             ,"Control"
                             ,javax.swing.border.TitledBorder.LEFT
                             ,javax.swing.border.TitledBorder.TOP
                             ,new java.awt.Font(null,java.awt.Font.BOLD,15)
                             ,java.awt.Color.black));
                             */

        startButton = new DeviceControllerButton(null);
        restartButton = new DeviceButton(null);
        resetButton = new DeviceButton(null);
        restartButton.setLabel("Restart");
        resetButton.setLabel("Reset averages");
        
        jp.add(startButton.graphic()); 
        jp.add(restartButton.graphic());  
        jp.add(resetButton.graphic());
                       
        setShape("VERTICAL");
    }
        
    /**
     * Sets the controller that is toggled by this device.
     */
    public void setSimulation(Simulation sim) {
        simulation = sim;
        Controller c = sim.getController();
        setController(c);
        startButton.setController(c);
        restartButton.setController(c);
        resetButton.setController(c);
        restartButton.setAction(new SimulationRestart(sim));
        resetButton.setAction(new SimulationDataAction(sim,new ResetAccumulators()));
    }
    
    /**
     * Returns the simulation affected by the buttons.
     */
    public Simulation getSimulation() {
        return simulation;
    }

    /**
     * Returns the JPanel that contains the buttons.
     */
    public java.awt.Component graphic(Object obj) {
        return jp;
    }
    
    public DeviceControllerButton getControllerButton() {return startButton;}
    public DeviceButton getRestartButton() {return restartButton;}
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
        if(s=="HORIZONTAL"){jp.setLayout(new java.awt.GridLayout(1,0));jp.updateUI();}
        if(s=="VERTICAL"){jp.setLayout(new java.awt.GridLayout(0,1));jp.updateUI();}
        if(s=="AUTOMATIC"){jp.getParent().addComponentListener(new ComponentEventControllerButton());}
    }
    public String getShape() {return shape;}
    
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
        
        HSMD2D sim = new HSMD2D(); 

        DeviceTrioControllerButton button = new DeviceTrioControllerButton(sim);
            button.setShape("HORIZONTAL"); //three choices "HORIZONTAL", "AUTOMATIC"          
//        DeviceTrioControllerButton button = new DeviceTrioControllerButton(Simulation.instance, Simulation.instance.controller(0)); 
//          button.setShape("VERTICAL"); //three choices "HORIZONTAL", "AUTOMATIC"
        
        SimulationGraphic graphic = new SimulationGraphic(sim);
        graphic.add(button);
        graphic.makeAndDisplayFrame();
    }
    
       
}//end DeviceTrioControllerButton
