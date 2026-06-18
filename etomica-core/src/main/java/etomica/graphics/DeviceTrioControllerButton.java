/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.IAction;
import etomica.action.SimulationRestart;
import etomica.action.controller.Controller;
import etomica.data.DataPump;
import etomica.simulation.Simulation;
import etomica.space.Space;

import javax.swing.*;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.ArrayList;

/**
 * Device comprising three buttons: (1) attaches to a controller to toggle its pause/resume state; 
 * (2)causes a reset (restart) action of simulation to be performed; (3) causes reset of meter averages
 * only.
 *
 * @author Sang Kyu Kwak
 */
 
public class DeviceTrioControllerButton extends Device {
    
    private final JPanel jp;
    private final DeviceRunControls runControls;
    private final Simulation simulation;
    private final DeviceButton reinitButton;
    private final DeviceButton resetButton;
	private double width;
	private boolean firstResized = true;
	private String shape;
	private final SimulationRestart simRestart;

    /**
     * Contructs device with buttons that affect the given simulation.
     */
    public DeviceTrioControllerButton(Simulation simulation, Space space, Controller controller) {
        super(controller);
        this.simulation = simulation;
        this.simRestart = new SimulationRestart(simulation);
        this.controller = controller;
        resetButton = new DeviceButton(controller);
        reinitButton = new DeviceButton(controller);
        this.runControls = new DeviceRunControls(controller);

        reinitButton.setLabel("Reinitialize");
        resetButton.setLabel("Reset averages");
        reinitButton.setPreAction(new IAction() {
            public void actionPerformed() {
                controller.pause().whenComplete((res, ex) -> {
                    runControls.reset();
                });
            }
        });
        reinitButton.setAction(simRestart);
        resetButton.setAction(simRestart.getDataResetAction());
        jp = new JPanel(new java.awt.GridLayout(1,0, 20, 20)); //default shape of panel
//        jp.setBorder(new TitledBorder(null, "Control", TitledBorder.CENTER, TitledBorder.TOP));
        jp.setOpaque(false);


//        jp.add(startButton.graphic());
        jp.add(runControls.graphic());
        jp.add(reinitButton.graphic());
        jp.add(resetButton.graphic());

        setShape("HORIZONTAL");
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
    public java.awt.Component graphic() {
        return jp;
    }
    
    public DeviceRunControls getRunControls() {return runControls;}
    public DeviceButton getReinitButton() {return reinitButton;}
    public DeviceButton getResetAveragesButton() {return resetButton;}
    
    /**
     * Sets controller toggle button to read "Start"
     */
    public void reset() {runControls.reset();}
    
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
}//end DeviceTrioControllerButton
