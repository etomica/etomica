/**
 * SimulateActions
 *
 * The SimulateActions class creates static action listeners to the simulate drop-down menu of the
 * EtomicaMenuBar.
 *
 * @author Bryan C. Mihalick
 * 8/16/00
 */

package etomica.gui;

import java.awt.Component;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDesktopPane;
import javax.swing.JFrame;
import javax.swing.JRadioButton;
import javax.swing.JToolBar;
import javax.swing.JInternalFrame;
import java.beans.PropertyVetoException;
import etomica.*;

public class SimulateActions {
        
    /**
     * Static action listener for opening the editor window corresponding to the current simulation
     */
    public static final ActionListener EDITSIMULATION = new EditSimulationAction();
    
    /**
     * Static action listener for performing the necessary actions when a new simulation is selected.
     */
    public static final ActionListener SELECTSIMULATIONACTION = new SelectSimulationAction(Simulation.instance);
    
    /**
     * Static class that handles the edit simulation action
     */
    private static class EditSimulationAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            EtomicaMenuBar.editSimulationItem.setEnabled(false);
            Etomica.simulationEditorFrame.setSimulationEditor(new SimulationEditor(Simulation.instance));
            try {
                Etomica.simulationEditorFrame.setClosed(false);
                Etomica.DesktopFrame.desktop.add(Etomica.simulationEditorFrame);
                Etomica.simulationEditorFrame.setSelected(true);
            }
            catch(java.beans.PropertyVetoException pve){}
        }// end of actionPerformed
    }// end of EditSimulationAction class
    
    /**
     * Class describing actions to be performed if a different simulation is selected.
     * Each instantiated simulation has its own instance of this action, which is called
     * if it becomes the current simulation.
     */
    public static class SelectSimulationAction implements ActionListener {
        Simulation simulation;
        /**
         * Constructor takes as argument the simulation associated with this action.
         */
        public SelectSimulationAction(Simulation sim) {
            simulation = sim;
        }
        public void actionPerformed(ActionEvent event) {
            Etomica.selectSimulation(simulation);
        }
    }
}// end of SimulateActions class
    