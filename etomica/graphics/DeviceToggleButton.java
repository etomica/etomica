//This class includes a main method to demonstrate its use
package etomica.graphics;
import javax.swing.JButton;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.modifier.ModifierBoolean;

/**
 * Button that toggles a boolean value.  This device can connect to any object
 * capable of switching between two states.  The device operates through a 
 * ModifierBoolean instance that must be connected to the state of the
 * controlled object.
 * 
 * @author David Kofke
 */
public class DeviceToggleButton extends Device implements EtomicaElement {
    
    private ModifierBoolean modifier;
    private JButton button;
    private String trueLabel = "True";
    private String falseLabel = "False";
    private boolean currentValue = false;
    
    public DeviceToggleButton(ModifierBoolean modifier) {
        this(modifier, "True", "False");
    }
    public DeviceToggleButton(ModifierBoolean modifier, 
                                String trueText, String falseText) {
        super();
        button = new JButton();
        button.addActionListener( new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                toggle();
            }
        });
        currentValue = modifier.getBoolean();
        setModifier(modifier);
        setTrueLabel(trueText);
        setFalseLabel(falseText);//sets label to correct state
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Button that toggles a boolean value");
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
        if(currentValue == true) button.setText(trueLabel);
        else button.setText(falseLabel);
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
     * Returns the GUI button element for display in the simulation.
     */
    public java.awt.Component graphic(Object obj) {
        return button;
    }
    
    /**
     * Specifies the button's label when the toggle is set to true.
     */
    public void setTrueLabel(String text) {
        trueLabel = text;
        if(currentValue == true) button.setText(trueLabel);
    }
    /**
     * @return the label displayed when the toggle is set to true.
     */
    public String getTrueLabel() {return trueLabel;}
    /**
     * Specifies the button's label when the toggle is set to false.
     */
    public void setFalseLabel(String text) {
        falseLabel = text;
        if(currentValue == false) button.setText(falseLabel);
    }
    /**
     * @return the label displayed when the toggle is set to false.
     */
    public String getFalseLabel() {return falseLabel;}
    
    /**
     * @return a handle to the JButton instance used by this device
     */
    public JButton getButton() {return button;}
        
    /**
     * Method to demonstrate and test the use of this class.  
     * Slider is used to control the temperature of a hard-sphere MD simulation
     */
/*    public static void main(String[] args) {
        
    //see DeviceCheckBox for an updated version of this method
    
        final etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        
        //here's the part unique to this class
        //sets up button to toggle atoms between red and blue
        sim.display.setColorScheme(new ColorSchemeNull());
        ModifierBoolean modifier = new ModifierBoolean() {
            public void setBoolean(boolean b) {
                if(b) sim.species.allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.setColor(java.awt.Color.red);}});
                else  sim.species.allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.setColor(java.awt.Color.blue);}});
                sim.panel().repaint();
            }
            public boolean getBoolean() {return sim.phase.firstAtom().getColor() == java.awt.Color.red;}
        };
        DeviceToggleButton button = new DeviceToggleButton(sim, modifier);
        button.setTrueLabel("Red");
        button.setFalseLabel("Blue");
        //end of unique part
 
        sim.elementCoordinator.go();
        Simulation.makeAndDisplayFrame(sim);
    }
    */
}
