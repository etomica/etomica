/**
 * SimulationEditor
 *
 * The SimulationEditor class is responsible for creating handles to the splitpanes contained 
 * in the tabs, as well as adding them to the tabs.  
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */
 
package etomica.gui;

import etomica.*;
import etomica.SimulationElement;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

//Java2 imports
//import java.util.Iterator;

import etomica.utility.Iterator;

public class SimulationEditor extends javax.swing.JTabbedPane {
    /**
     * Static splitpane that displays all the potential 1 classes on the left (as radiobuttons) and 
     * all the potentials added to the simulation.instance on the right (in a JList).
     */
    public final Potential1EditorPane potential1Editor;
    
    /**
     * Static splitpane that displays all the potential 2 classes on the left (as radiobuttons) and 
     * all the potentials added to the simulation.instance on the right (in a JList).
     */
    public final Potential2EditorPane potential2Editor;
    
    /**
     * Static splitpane that displays all the species classes on the left (as radiobuttons) and 
     * all the species added to the simulation.instance on the right (in a JList).
     */
    public final SimulationEditorPane speciesEditor;
    
    /**
     * Static splitpane that displays all the integrator classes on the left (as radiobuttons) and 
     * all the integrator added to the simulation.instance on the right (in a JList).
     */
    public final SimulationEditorPane integratorEditor;
    
    /**
     * Static splitpane that displays all the phase classes on the left (as radiobuttons) and 
     * all the phase added to the simulation.instance on the right (in a JList).
     */
    public final SimulationEditorPane phaseEditor;
    
    /**
     * Static splitpane that displays all the controller classes on the left (as radiobuttons) and 
     * all the controller added to the simulation.instance on the right (in a JList).
     */
    public final SimulationEditorPane controllerEditor;
    
    /**
     * Static splitpane that displays all the display classes on the left (as radiobuttons) and 
     * all the display added to the simulation.instance on the right (in a JList).
     */
    public final SimulationEditorPane displayEditor;
    
    /**
     * Static splitpane that displays all the meters classes on the left (as radiobuttons) and 
     * all the meter added to the simulation.instance on the right (in a JList).
     */
    public final SimulationEditorPane meterEditor;
    
    /**
     * Static splitpane that displays all the device classes on the left (as radiobuttons) and 
     * all the device added to the simulation.instance on the right (in a JList).
     */
    public final SimulationEditorPane deviceEditor;
    
    /**
     * Instance of the simulation being edited by this tab pane.
     */
    private Simulation simulation;
    
    private java.util.HashMap editorPanes = new java.util.HashMap(16);

    /**
     * Constructor that adds all the splitpanes to the tabmenu of the SimulationEditorFrame
     */
    public SimulationEditor(Simulation sim){
        simulation = sim;
        potential1Editor = new Potential1EditorPane(this);
        potential2Editor = new Potential2EditorPane(this);
        speciesEditor = new SimulationEditorPane(this,"Species");
        integratorEditor = new SimulationEditorPane(this,"Integrator");
        phaseEditor = new SimulationEditorPane(this,"Phase");
        controllerEditor = new SimulationEditorPane(this,"Controller");
        displayEditor = new SimulationEditorPane(this,"Display");
        meterEditor = new SimulationEditorPane(this,"Meter");
        deviceEditor = new SimulationEditorPane(this,"Device");
        editorPanes.put(Potential1.class, potential1Editor);
        editorPanes.put(Potential2.class, potential2Editor);
        editorPanes.put(Species.class, speciesEditor);
        editorPanes.put(Integrator.class, integratorEditor);
        editorPanes.put(Phase.class, phaseEditor);
        editorPanes.put(Controller.class, controllerEditor);
        editorPanes.put(Display.class, displayEditor);
        editorPanes.put(MeterAbstract.class, meterEditor);
        editorPanes.put(Device.class, deviceEditor);
        addTab("Species", speciesEditor);
        addTab("Potential1", potential1Editor);
        addTab("Potential2", potential2Editor);
        addTab("Integrator", integratorEditor);
        addTab("Phase", phaseEditor);
        addTab("Controller", controllerEditor);
        addTab("Display", displayEditor);
        addTab("Meter", meterEditor);
        addTab("Device", deviceEditor);
        updateElements();
        checkSimFeasibility();
    }// end of SimEditorTabMenu constructor
    
    public Simulation getSimulation() {return simulation;}
    
    /** Check if a sufficient number of components are added to allow a working simulation.
     *  If so, enable the start button.
     */
    public final void checkSimFeasibility(){
	    if (allRemoveEnabled()) setAllStart(true);
	    else setAllStart(false);
    }
    public void updateElements() {
        resetAllComponentLists();
        for(Iterator iter=simulation.allElements().iterator(); iter.hasNext(); ) {
            SimulationElement element = (SimulationElement)iter.next();
            EditorPane editorPane = (EditorPane)editorPanes.get(element.baseClass());
            editorPane.getComponentList().addElement(element);
            editorPane.setEnabled(true);
        }
        speciesEditor.accountForNewSpecies();
    }
    
    public void setAllRemove(boolean r) { 
	    potential1Editor.remove.setEnabled(r);
	    potential2Editor.remove.setEnabled(r);
	    speciesEditor.remove.setEnabled(r);
	    integratorEditor.remove.setEnabled(r);
	    phaseEditor.remove.setEnabled(r);
	    controllerEditor.remove.setEnabled(r);
	    displayEditor.remove.setEnabled(r);
	    meterEditor.remove.setEnabled(r);
	    deviceEditor.remove.setEnabled(r);
    }
    
    public boolean allRemoveEnabled() {
        return (potential2Editor.remove.isEnabled() && speciesEditor.remove.isEnabled() && 
	        integratorEditor.remove.isEnabled() && phaseEditor.remove.isEnabled() && 
	        controllerEditor.remove.isEnabled() && displayEditor.remove.isEnabled());
    }
    
    public void setAllStart(boolean d) { 
	    potential1Editor.start.setEnabled(d);
	    potential2Editor.start.setEnabled(d);
	    speciesEditor.start.setEnabled(d);
	    integratorEditor.start.setEnabled(d);
	    phaseEditor.start.setEnabled(d);
	    controllerEditor.start.setEnabled(d);
	    displayEditor.start.setEnabled(d);
	    meterEditor.start.setEnabled(d);
	    deviceEditor.start.setEnabled(d);
    }

    public boolean allStartEnabled() {
        if (potential2Editor.start.isEnabled() && speciesEditor.start.isEnabled() && 
	        integratorEditor.start.isEnabled() && phaseEditor.start.isEnabled() && 
	        controllerEditor.start.isEnabled() && displayEditor.start.isEnabled())
	        return true;
        else return false;
    }
    
    public void resetAllComponentLists(){
        for(java.util.Iterator iter=editorPanes.values().iterator(); iter.hasNext(); ) {
            ((EditorPane)iter.next()).getComponentList().removeAllElements();
        }
    }
    
    public SimulationEditorPane getSpeciesEditor(){ return speciesEditor; }
    public Potential1EditorPane getPotential1Editor(){ return potential1Editor; }
    public Potential2EditorPane getPotential2Editor(){ return potential2Editor; }
    public SimulationEditorPane getIntegratorEditor(){ return integratorEditor; }
    public SimulationEditorPane getPhaseEditor(){ return phaseEditor; }
    public SimulationEditorPane getControllerEditor(){ return controllerEditor; }
    public SimulationEditorPane getDisplayEditor(){ return displayEditor; }
    public SimulationEditorPane getMeterEditor(){ return meterEditor; }
    public SimulationEditorPane getDeviceEditor(){ return deviceEditor; }
}// end of SimEditorTabMenu class