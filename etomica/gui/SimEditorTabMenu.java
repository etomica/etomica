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

import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class SimEditorTabMenu extends javax.swing.JTabbedPane {
    /**
     * Static splitpane that displays all the potential 1 classes on the left (as radiobuttons) and 
     * all the potentials added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ Potential1EditorPane potential1Editor = new Potential1EditorPane();
    
    /**
     * Static splitpane that displays all the potential 2 classes on the left (as radiobuttons) and 
     * all the potentials added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ Potential2EditorPane potential2Editor = new Potential2EditorPane();
    
    /**
     * Static splitpane that displays all the species classes on the left (as radiobuttons) and 
     * all the species added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ SpeciesEditorPane speciesEditor = new SpeciesEditorPane();
    
    /**
     * Static splitpane that displays all the integrator classes on the left (as radiobuttons) and 
     * all the integrator added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ IntegratorEditorPane integratorEditor = new IntegratorEditorPane();
    
    /**
     * Static splitpane that displays all the phase classes on the left (as radiobuttons) and 
     * all the phase added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ PhaseEditorPane phaseEditor = new PhaseEditorPane();
    
    /**
     * Static splitpane that displays all the controller classes on the left (as radiobuttons) and 
     * all the controller added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ ControllerEditorPane controllerEditor = new ControllerEditorPane();
    
    /**
     * Static splitpane that displays all the display classes on the left (as radiobuttons) and 
     * all the display added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ DisplayEditorPane displayEditor = new DisplayEditorPane();
    
    /**
     * Static splitpane that displays all the meters classes on the left (as radiobuttons) and 
     * all the meter added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ MeterEditorPane meterEditor = new MeterEditorPane();
    
    /**
     * Static splitpane that displays all the device classes on the left (as radiobuttons) and 
     * all the device added to the simulation.instance on the right (in a JList).
     */
    static /*final*/ DeviceEditorPane deviceEditor = new DeviceEditorPane();

    /**
     * Constructor that add all the splitpane to the tabmenu of the SimulationEditorFrame
     */
    public SimEditorTabMenu(){
        addTab("Species", speciesEditor);
        addTab("Potential1", potential1Editor);
        addTab("Potential2", potential2Editor);
        addTab("Integrator", integratorEditor);
        addTab("Phase", phaseEditor);
        addTab("Controller", controllerEditor);
        addTab("Display", displayEditor);
        addTab("Meter", meterEditor);
        addTab("Device", deviceEditor);        
    }// end of SimEditorTabMenu constructor
    
    private void writeObject(ObjectOutputStream out) throws java.io.IOException {
        out.defaultWriteObject();
        out.writeObject(speciesEditor);
        out.writeObject(potential1Editor);
        out.writeObject(potential2Editor);
        out.writeObject(integratorEditor);
        out.writeObject(phaseEditor);
        out.writeObject(controllerEditor);
        out.writeObject(displayEditor);
        out.writeObject(meterEditor);
        out.writeObject(deviceEditor);
    }
    
    private void readObject(ObjectInputStream in) throws java.io.IOException, java.lang.ClassNotFoundException {
        in.defaultReadObject();
        speciesEditor = (SpeciesEditorPane)in.readObject();
        potential1Editor = (Potential1EditorPane)in.readObject();
        potential2Editor = (Potential2EditorPane)in.readObject();
        integratorEditor = (IntegratorEditorPane)in.readObject();
        phaseEditor = (PhaseEditorPane)in.readObject();
        controllerEditor = (ControllerEditorPane)in.readObject();
        displayEditor = (DisplayEditorPane)in.readObject();
        meterEditor = (MeterEditorPane)in.readObject();
        deviceEditor = (DeviceEditorPane)in.readObject();
    }    
    
    public static void setAllRemove(boolean r) { 
	    SimEditorTabMenu.potential1Editor.remove.setEnabled(r);
	    SimEditorTabMenu.potential2Editor.remove.setEnabled(r);
	    SimEditorTabMenu.speciesEditor.remove.setEnabled(r);
	    SimEditorTabMenu.integratorEditor.remove.setEnabled(r);
	    SimEditorTabMenu.phaseEditor.remove.setEnabled(r);
	    SimEditorTabMenu.controllerEditor.remove.setEnabled(r);
	    SimEditorTabMenu.displayEditor.remove.setEnabled(r);
	    SimEditorTabMenu.meterEditor.remove.setEnabled(r);
	    SimEditorTabMenu.deviceEditor.remove.setEnabled(r);
    }
    
    public static boolean allRemoveEnabled() {
        if (potential2Editor.remove.isEnabled() && speciesEditor.remove.isEnabled() && 
	        integratorEditor.remove.isEnabled() && phaseEditor.remove.isEnabled() && 
	        controllerEditor.remove.isEnabled() && displayEditor.remove.isEnabled())
	        return true;
        else return false;
    }
    
    public static void setAllStart(boolean d) { 
	    SimEditorTabMenu.potential1Editor.start.setEnabled(d);
	    SimEditorTabMenu.potential2Editor.start.setEnabled(d);
	    SimEditorTabMenu.speciesEditor.start.setEnabled(d);
	    SimEditorTabMenu.integratorEditor.start.setEnabled(d);
	    SimEditorTabMenu.phaseEditor.start.setEnabled(d);
	    SimEditorTabMenu.controllerEditor.start.setEnabled(d);
	    SimEditorTabMenu.displayEditor.start.setEnabled(d);
	    SimEditorTabMenu.meterEditor.start.setEnabled(d);
	    SimEditorTabMenu.deviceEditor.start.setEnabled(d);
    }

    public static boolean allStartEnabled() {
        if (potential2Editor.start.isEnabled() && speciesEditor.start.isEnabled() && 
	        integratorEditor.start.isEnabled() && phaseEditor.start.isEnabled() && 
	        controllerEditor.start.isEnabled() && displayEditor.start.isEnabled())
	        return true;
        else return false;
    }
    
    public static void resetAllComponentLists(){
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
    
    public static SpeciesEditorPane getSpeciesEditor(){ return speciesEditor; }
    public static Potential1EditorPane getPotential1Editor(){ return potential1Editor; }
    public static Potential2EditorPane getPotential2Editor(){ return potential2Editor; }
    public static IntegratorEditorPane getIntegratorEditor(){ return integratorEditor; }
    public static PhaseEditorPane getPhaseEditor(){ return phaseEditor; }
    public static ControllerEditorPane getControllerEditor(){ return controllerEditor; }
    public static DisplayEditorPane getDisplayEditor(){ return displayEditor; }
    public static MeterEditorPane getMeterEditor(){ return meterEditor; }
    public static DeviceEditorPane getDeviceEditor(){ return deviceEditor; }
}// end of SimEditorTabMenu class