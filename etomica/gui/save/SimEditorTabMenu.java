/**
 * SimEditorTabMenu
 *
 * The SimEditorTabMenu class is responsible for creating static handles to the splitpanes contained 
 * in the tabs, as well as adding them to the tabs.  
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */
 
package simulate.gui;

import simulate.Simulation;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class SimEditorTabMenu extends javax.swing.JTabbedPane {
    /**
     * Static splitpane that displays all the potential 1 classes on the left (as radiobuttons) and 
     * all the potentials added to the simulation.instance on the right (in a JList).
     */
    public final Potential1EditorPane potential1Editor = new Potential1EditorPane();
    
    /**
     * Static splitpane that displays all the potential 2 classes on the left (as radiobuttons) and 
     * all the potentials added to the simulation.instance on the right (in a JList).
     */
    public final Potential2EditorPane potential2Editor = new Potential2EditorPane();
    
    /**
     * Static splitpane that displays all the species classes on the left (as radiobuttons) and 
     * all the species added to the simulation.instance on the right (in a JList).
     */
    public final SpeciesEditorPane speciesEditor = new SpeciesEditorPane();
    
    /**
     * Static splitpane that displays all the integrator classes on the left (as radiobuttons) and 
     * all the integrator added to the simulation.instance on the right (in a JList).
     */
    public final IntegratorEditorPane integratorEditor = new IntegratorEditorPane();
    
    /**
     * Static splitpane that displays all the phase classes on the left (as radiobuttons) and 
     * all the phase added to the simulation.instance on the right (in a JList).
     */
    public final PhaseEditorPane phaseEditor = new PhaseEditorPane();
    
    /**
     * Static splitpane that displays all the controller classes on the left (as radiobuttons) and 
     * all the controller added to the simulation.instance on the right (in a JList).
     */
    public final ControllerEditorPane controllerEditor = new ControllerEditorPane();
    
    /**
     * Static splitpane that displays all the display classes on the left (as radiobuttons) and 
     * all the display added to the simulation.instance on the right (in a JList).
     */
    public final DisplayEditorPane displayEditor = new DisplayEditorPane();
    
    /**
     * Static splitpane that displays all the meters classes on the left (as radiobuttons) and 
     * all the meter added to the simulation.instance on the right (in a JList).
     */
    public final MeterEditorPane meterEditor = new MeterEditorPane();
    
    /**
     * Static splitpane that displays all the device classes on the left (as radiobuttons) and 
     * all the device added to the simulation.instance on the right (in a JList).
     */
    public final DeviceEditorPane deviceEditor = new DeviceEditorPane();
    
    /**
     * Instance of the simulation being edited by this tab pane.
     */
    public Simulation simulation;

    /**
     * Constructor that add all the splitpane to the tabmenu of the SimulationEditorFrame
     */
    public SimEditorTabMenu(Simulation sim){
        addTab("Species", speciesEditor);
        addTab("Potential1", potential1Editor);
        addTab("Potential2", potential2Editor);
        addTab("Integrator", integratorEditor);
        addTab("Phase", phaseEditor);
        addTab("Controller", controllerEditor);
        addTab("Display", displayEditor);
        addTab("Meter", meterEditor);
        addTab("Device", deviceEditor);
        sim.setEditor(this);
        simulation = sim;
    }// end of SimEditorTabMenu constructor
    
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
        if (potential2Editor.remove.isEnabled() && speciesEditor.remove.isEnabled() && 
	        integratorEditor.remove.isEnabled() && phaseEditor.remove.isEnabled() && 
	        controllerEditor.remove.isEnabled() && displayEditor.remove.isEnabled())
	        return true;
        else return false;
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
        potential1Editor.componentList.removeAllElements();
        potential2Editor.componentList.removeAllElements();
        speciesEditor.componentList.removeAllElements();
        integratorEditor.componentList.removeAllElements();
        phaseEditor.componentList.removeAllElements();
        controllerEditor.componentList.removeAllElements();
        displayEditor.componentList.removeAllElements();
        meterEditor.componentList.removeAllElements();
        deviceEditor.componentList.removeAllElements();
    }
    
    public SpeciesEditorPane getSpeciesEditor(){ return speciesEditor; }
    public Potential1EditorPane getPotential1Editor(){ return potential1Editor; }
    public Potential2EditorPane getPotential2Editor(){ return potential2Editor; }
    public IntegratorEditorPane getIntegratorEditor(){ return integratorEditor; }
    public PhaseEditorPane getPhaseEditor(){ return phaseEditor; }
    public ControllerEditorPane getControllerEditor(){ return controllerEditor; }
    public DisplayEditorPane getDisplayEditor(){ return displayEditor; }
    public MeterEditorPane getMeterEditor(){ return meterEditor; }
    public DeviceEditorPane getDeviceEditor(){ return deviceEditor; }
}// end of SimEditorTabMenu class