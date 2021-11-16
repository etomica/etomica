/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

//This class includes a main method to demonstrate its use
package etomica.graphics;

import etomica.action.controller.Controller;
import etomica.modifier.ModifierBoolean;

import javax.swing.*;

/**
 * Button that toggles a boolean value using a pair of radio buttons.
 * This device can connect to any object
 * capable of switching between two states.  The device operates through a
 * ModifierBoolean instance that must be connected to the state of the
 * controlled object.
 *
 * @author David Kofke
 */
public class DeviceToggleRadioButtons extends Device {
    
    private ModifierBoolean modifier;
    private JPanel panel;
    private JRadioButton trueButton, falseButton;

    /**
     * @param modifier the boolean modifier controlled by this device
     * @param title     a descriptive string.  If empty ("") provides plain border; if null, provides no border.
     * @param trueText  text associated with "true" state of modifier
     * @param falseText text associated with "false" state of modifier
     */
    public DeviceToggleRadioButtons(Controller controller, final ModifierBoolean modifier,
                                    String title, String trueText, String falseText) {
        super(controller);
        java.awt.event.ActionListener al = new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent e) {
                modifier.setBoolean(getState());
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

        falseButton.setSelected(!modifier.getBoolean());
        setModifier(modifier);
        
        if(title != null /*&& !title.equals("")*/) setTitle(title);
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
     * Specifies the boolean modifier that is set between true and false by the button.
     */
    public ModifierBoolean getModifier() {return modifier;}
    /**
     * Returns the boolean modifier set by this button.
     */
    public void setModifier(ModifierBoolean newModifier) {
        modifier = newModifier;
        modifier.setBoolean(getState());
    }
    
    /**
     * Returns the GUI button element for display in the simulation.
     */
    public java.awt.Component graphic() {
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
        
        final etomica.simulation.prototypes.HSMD2D sim = new etomica.simulation.prototypes.HSMD2D();
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
            public boolean getBoolean() {return sim.box.firstAtom().getColor() == java.awt.Color.red;}
        };
        DeviceToggleRadioButtons selector = new DeviceToggleRadioButtons(sim, modifier);
        selector.setTrueLabel("Red");
        selector.setFalseLabel("Blue");
//        selector.setTitle("Color");
        //end of unique part
 
        sim.elementCoordinator.go();
        Simulation.makeAndDisplayFrame(sim);
    }
    */
}
