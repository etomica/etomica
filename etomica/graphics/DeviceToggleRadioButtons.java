//This class includes a main method to demonstrate its use
package etomica.graphics;
import etomica.*;
import javax.swing.*;

/**
 * Button that toggles a boolean value using a pair of radio buttons.  
 * This device can connect to any object
 * capable of switching between two states.  The device operates through a 
 * ModulatorBoolean instance that must be connected to the state of the
 * controlled object.
 * 
 * @author David Kofke
 */
public class DeviceToggleRadioButtons extends Device implements EtomicaElement {
    
    private ModulatorBoolean modulator;
    private JPanel panel;
    private JRadioButton trueButton, falseButton;
    
    /**
     * Constructor with default labels of a blank title and "True" and "False" for
     * the true/false labels.
     */
    public DeviceToggleRadioButtons(ModulatorBoolean modulator) {
        this(modulator, "", "True", "False");
    }
    
    /**
     * @param sim       the parent simulation of this device
     * @param modulator the boolean modulator controlled by this device
     * @param title     a descriptive string.  If empty ("") provides plain border; if null, provides no border.
     * @param trueText  text associated with "true" state of modulator
     * @param falseText text associated with "false" state of modulator
     */
    public DeviceToggleRadioButtons(final ModulatorBoolean modulator, 
                                String title, String trueText, String falseText) {

        java.awt.event.ActionListener al = new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent e) {
                modulator.setBoolean(getState());
            }
        };
  
        trueButton = new JRadioButton(trueText, true);
        falseButton = new JRadioButton(falseText, false);
        trueButton.addActionListener(al);
        falseButton.addActionListener(al);
        ButtonGroup g = new ButtonGroup();
        g.add(trueButton);
        g.add(falseButton);
          
        panel = new JPanel();
        panel.add(trueButton);
        panel.add(falseButton);

        falseButton.setSelected(!modulator.getBoolean());
        setModulator(modulator);
        
        if(title != null /*&& !title.equals("")*/) setTitle(title);
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
        trueButton.setSelected(b);
        falseButton.setSelected(!b);
    }
    /**
     * Returns the currents true/false state of the button.
     */
    public boolean getState() {return trueButton.isSelected();}
        
    /**
     * Specifies the boolean modulator that is set between true and false by the button.
     */
    public ModulatorBoolean getModulator() {return modulator;}
    /**
     * Returns the boolean modulator set by this button.
     */
    public void setModulator(ModulatorBoolean newModulator) {
        modulator = newModulator;
        modulator.setBoolean(getState());
    }
    
    /**
     * Returns the GUI button element for display in the simulation.
     */
    public java.awt.Component graphic(Object obj) {
        return panel;
    }
    
    public void setTitle(String text) {
        panel.setBorder(new javax.swing.border.TitledBorder(text));
    }
    public String getTitle() {
        return ((javax.swing.border.TitledBorder)panel.getBorder()).getTitle();
    }
    
    
    /**
     * Specifies the button's label when the toggle is set to true.
     */
    public void setTrueLabel(String text) {
        trueButton.setText(text);
    }
    /**
     * @return the label displayed when the toggle is set to true.
     */
    public String getTrueLabel() {return trueButton.getText();}
    /**
     * Specifies the button's label when the toggle is set to false.
     */
    public void setFalseLabel(String text) {
        falseButton.setText(text);
    }
    /**
     * @return the label displayed when the toggle is set to false.
     */
    public String getFalseLabel() {return falseButton.getText();}
    
    /**
     * @return the instance of the button selecting the "true" state.
     * Useful to change it color or other graphics attributes.
     */
    public JRadioButton trueButton() {return trueButton;}
    /**
     * @return the instance of the button selecting the "false" state.
     * Useful to change it color or other graphics attributes.
     */
    public JRadioButton falseButton() {return falseButton;}
            
    /**
     * Method to demonstrate and test the use of this class.  
     * Buttons are used to toggle the color of atoms in a simulation.
     */
/*    public static void main(String[] args) {
        
        final etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        
        //here's the part unique to this class
        //sets up button to toggle atoms between red and blue
        sim.display.setColorScheme(new ColorSchemeNull());
        ModulatorBoolean modulator = new ModulatorBoolean() {
            public void setBoolean(boolean b) {
                if(b) sim.species.allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.setColor(java.awt.Color.red);}});
                else  sim.species.allAtoms(new AtomAction() {public void actionPerformed(Atom a) {a.setColor(java.awt.Color.blue);}});
                sim.panel().repaint();
            }
            public boolean getBoolean() {return sim.phase.firstAtom().getColor() == java.awt.Color.red;}
        };
        DeviceToggleRadioButtons selector = new DeviceToggleRadioButtons(sim, modulator);
        selector.setTrueLabel("Red");
        selector.setFalseLabel("Blue");
//        selector.setTitle("Color");
        //end of unique part
 
        sim.elementCoordinator.go();
        Simulation.makeAndDisplayFrame(sim);
    }
    */
}