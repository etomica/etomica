/**
 * ControllerActions
 *
 * The ControllerActions class is responsible for creating static action listeners to the controller 
 * etomica.Controller class.
 *
 * @author Bryan C. Mihalick
 * 1/17/01
 */

package etomica.gui;

import etomica.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.Component;
import javax.swing.JEditorPane;
import javax.swing.JInternalFrame;
import javax.swing.JScrollPane;
import java.io.IOException;

public class ControllerActions {
    /**
     * Static action listener that starts the simulation
     */
    public static final ActionListener START = new StartAction();
    
    /**
     * Static action listener that stops the simulation
     */
    public static final ActionListener STOP = new StopAction();
    
    /**
     * Static action listener that pauses the simulation
     */
    public static final ActionListener PAUSE = new PauseAction();

    /**
     * Handles the start event and calls the start method of the Controller class
     */
    private static class StartAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            ((Controller)Simulation.instance.controllerList().get(0)).start();
            ((Controller.Button)((Controller)Simulation.instance.controllerList().get(0)).getButton()).clickForUnpause();
        }// end of actionPerformed
    }// end of StartAction

    /**
     * Handles the stop event and calls the stop method of the Controller class
     */
    private static class StopAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            ((Controller)Simulation.instance.controllerList().get(0)).halt();
            ((Controller.Button)((Controller)Simulation.instance.controllerList().get(0)).getButton()).reset();
        }// end of actionPerformed
    }// end of StopAction

    /**
     * Handles the pause event and calls the pause method of the Controller class
     */
    private static class PauseAction implements ActionListener {
        
        //need to revise so that attempt to pause cannot cause a start
        public void actionPerformed(ActionEvent event) {
            ((Controller.Button)((Controller)Simulation.instance.controllerList().get(0)).getButton()).button.doClick();
        }// end of actionPerformed
    }// end of PauseAction
}// end of ControllerActions class
    