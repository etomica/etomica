//This class includes a main method to demonstrate its use
package etomica.graphics;
import javax.swing.JCheckBox;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.modifier.ModifierBoolean;

/**
 * Button that toggles a boolean value via a check box.  This device can connect to any object
 * capable of switching between two states.  The device operates through a 
 * ModifierBoolean instance that must be connected to the state of the
 * controlled object.
 * 
 * @author David Kofke
 */
public class DeviceCheckBox extends Device implements EtomicaElement {
    
    private ModifierBoolean modifier;
    private JCheckBox box;
    private boolean currentValue = false;
    
    public DeviceCheckBox(ModifierBoolean modifier) {
        this("Select", modifier);
    }
    public DeviceCheckBox(String label, ModifierBoolean modifier) {
        super();
        init(label, modifier);
    }
    
    private void init(String label, ModifierBoolean modifier) {
        this.modifier = modifier;
        currentValue = modifier.getBoolean();
        box = new JCheckBox(label,currentValue);
        box.addActionListener( new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                toggle();
            }
        });
    }
        
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Check box that toggles a boolean value");
        return info;
    }
    
    /**
     * Sets the button to the given state.
     */
    public void setState(boolean b) {
        if(b != currentValue) toggle();
    }
    /**
     * Returns the currents true/false state of the button.
     */
    public boolean getState() {return currentValue;}
    
    /**
     * Toggles the state of the button (true to false, false to true).
     */
    public void toggle() {
        currentValue = !currentValue;
        modifier.setBoolean(currentValue);
        box.setSelected(currentValue);
    }
    
    /**
     * Specifies the boolean modifier that is set between true and false by the button.
     */
    public ModifierBoolean getModifier() {return modifier;}
    /**
     * Returns the boolean modifier set by this button.
     */
    public void setModifier(ModifierBoolean newModifier) {
        modifier = newModifier;
        modifier.setBoolean(currentValue);
    }
    
    /**
     * Returns the GUI box element for display in the simulation.
     */
    public java.awt.Component graphic(Object obj) {
        return box;
    }
    
    
    /**
     * @return a handle to the JCheckBox instance used by this device
     */
    public JCheckBox getCheckBox() {return box;}
        
    /**
     * Method to demonstrate and test the use of this class.  
     */
//    public static void main(String[] args) {
//        
//        final etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
//        Simulation.instance = sim;
//        
//        //here's the part unique to this class
//        //sets up box to toggle atoms between red and blue
//        final ColorSchemeByType colorScheme = new ColorSchemeByType();
//        sim.display.setColorScheme(colorScheme);
//        ModifierBoolean modifier = new ModifierBoolean() {
//            public void setBoolean(boolean b) {
//                if(b) ColorSchemeByType.setColor((SpeciesSpheresMono)sim.species, java.awt.Color.blue);
//                else ColorSchemeByType.setColor((SpeciesSpheresMono)sim.species, java.awt.Color.red);
////                if(b) sim.species.allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.setColor(java.awt.Color.red);}});
////                else  sim.species.allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.setColor(java.awt.Color.blue);}});
//                sim.panel().repaint();
//            }
//            public boolean getBoolean() {return colorScheme.atomColor(sim.phase.firstAtom()) == java.awt.Color.blue;}
//        };
//        modifier.setBoolean(true);
//        DeviceCheckBox button = new DeviceCheckBox(sim, "Blue", modifier);
//        //end of unique part
// 
//        sim.elementCoordinator.go();
//        SimulationGraphic.makeAndDisplayFrame(sim);
//    }
//    */
}
