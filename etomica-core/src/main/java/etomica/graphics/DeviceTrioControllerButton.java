/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.ArrayList;

import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.activity.Controller;
import etomica.simulation.Simulation;
import etomica.data.DataPump;
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
    private Simulation simulation;
    private DeviceButton reinitButton;
    private DeviceButton resetButton;
	private double width;
	private boolean firstResized = true;
	private String shape;
	private SimulationRestart simRestart;

    /**
     * Contructs device with buttons that affect the given simulation.
     */
    public DeviceTrioControllerButton(Simulation simulation, Space space, Controller _controller) {
        this();
        setSimulation(simulation, space, _controller);
    }
    
    /**
     * No-argument contructor gives device that performs no action until
     * setSimulation is used to specify the target simulation.
     */
    public DeviceTrioControllerButton() {
        super();
        jp = new JPanel(new java.awt.GridLayout(1,0, 20, 20)); //default shape of panel
//        jp.setBorder(new TitledBorder(null, "Control", TitledBorder.CENTER, TitledBorder.TOP));
        jp.setOpaque(false);

        startButton = new DeviceControllerButton(null);
        reinitButton = new DeviceButton(null);
        resetButton = new DeviceButton(null);
        reinitButton.setLabel("Reinitialize");
        resetButton.setLabel("Reset averages");
        
        jp.add(startButton.graphic()); 
        jp.add(reinitButton.graphic());  
        jp.add(resetButton.graphic());
                       
        setShape("HORIZONTAL");
    }
        
    /**
     * Sets the controller that is toggled by this device.
     */
    protected void setSimulation(Simulation sim, Space space, Controller controller) {
        simulation = sim;
        simRestart = new SimulationRestart(sim, space, controller);
        final Controller c = controller;
        setController(c);
        startButton.setController(c);
        reinitButton.setPreAction(new IAction() {
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
    public ArrayList<DataPump> getDataStreamPumps() {
        return simRestart.getDataResetAction().getDataStreamPumps();
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

        DeviceTrioControllerButton button = new DeviceTrioControllerButton(sim, sp, sim.getController());
            button.setShape("HORIZONTAL"); //three choices "HORIZONTAL", "AUTOMATIC"          
//        DeviceTrioControllerButton button = new DeviceTrioControllerButton(Simulation.instance, Simulation.instance.controller(0)); 
//          button.setShape("VERTICAL"); //three choices "HORIZONTAL", "AUTOMATIC"
        
        final SimulationGraphic graphic = new SimulationGraphic(sim, APP_NAME, sp, sim.getController());

        // Simulation Graphic will display it's own Trio button group by
        // default.  Just remove them and put ours on for this test.
        graphic.getPanel().controlPanel.removeAll();
        graphic.add(button);

        button.getReinitButton().setPostAction(graphic.getPaintAction(sim.box));

        graphic.makeAndDisplayFrame(APP_NAME);
    }
    
       
}//end DeviceTrioControllerButton
